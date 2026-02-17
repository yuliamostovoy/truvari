"""Remap VCF alleles' sequence to the reference genome to annotate REMAP.

Classification:
    novel        : Allele has no hits in reference
    tandem       : Allele's closest hit is within len(allele) bp of the SV's position and same orientation
    tandem_inverted : Tandem-proximal hit but inverted orientation
    interspersed : Closest hit is not tandem
    partial      : Only partial hit(s) passing min score but < threshold coverage
    failed       : Alignment infrastructure failure (aligner error). Indicates retry needed.
    over_max_size: Skipped because ALT sequence exceeded --max-length cutoff

Parameters:
    -m/--min-length    Minimum SV length to be remapped (50)
    -M/--max-length    Maximum ALT sequence length (bp) to attempt remap (10000000)
    --cov-threshold    Minimum fraction of allele aligned (.8). For >100kb: use absolute threshold * 100kb
    --aligner          bwa (default), minimap2, or blastn
    --mm2-preset       minimap2 preset (asm20 default; options: asm5, asm10, asm20, map-ont, sr, etc.)
    --threads          Number of threads
    --mm2-threshold    Insertions > this size (bp) use minimap2 (5000)
    --save-alignments-prefix  Path prefix for writing <prefix>small.sam and <prefix>large.sam

REQUIREMENTS:
    * minimap2, bwa, or blastn binary in PATH (choose via --aligner)
    * Reference FASTA (indexed for chosen aligner)
    * When using blastn, provide a pre-built BLAST database (via --blast-db)
"""
import sys
import logging
import argparse
import subprocess
import tempfile
import os
import shlex
import hashlib
import json
from shutil import which

import truvari
from truvari.annotations.grm import cigmatch


def log_subprocess_failure(tool_name, proc):
    """Log exit details and stderr content from a failed subprocess."""
    rc = proc.returncode
    if rc < 0:
        logging.error(f"{tool_name} terminated by signal {abs(rc)}")
    else:
        logging.error(f"{tool_name} exited with code {rc}")
    stderr = (proc.stderr or "").strip()
    if stderr:
        logging.error(f"{tool_name} stderr:\n{stderr}")

def infer_query_size(qname, seq=None, fallback=0):
    """Infer query length from qname (q<chrom>_<pos>_<len>_<hash>) or fall back to sequence."""
    if qname:
        try:
            _, qlen_str, _ = qname.rsplit('_', 2)
            return int(qlen_str)
        except Exception:
            pass
    if seq and seq != '*':
        return len(seq)
    return fallback
def parse_sam_lines(lines, results, source="sam"):
    """Parse SAM lines into AlignmentHit objects and append into results dict keyed by qName.

    - Computes match_bases from CIGAR (M, =, X)
    - Computes end from ref-advancing ops (M, =, X, D, N)
    - Derives query size from qName pattern q<chrom>_<pos>_<len>_<hash>, else falls back to SEQ length or match_bases
    """
    try:
        for line in lines:
            if not line or line.startswith("@"):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 11:
                continue
            qName = parts[0]
            try:
                flag = int(parts[1])
            except ValueError:
                continue
            rName = parts[2]
            pos_str = parts[3]
            cigar = parts[5]
            seq = parts[9]
            # Skip unmapped
            if (flag & 0x4) != 0 or rName == '*':
                continue
            orientation = '-' if (flag & 0x10) != 0 else '+'
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            # match_bases from CIGAR
            try:
                match_bases = cigmatch(cigar)
            except Exception:
                match_bases = 0
                num = ''
                for ch in cigar:
                    if ch.isdigit():
                        num += ch
                    else:
                        if ch in ('M', '=', 'X') and num:
                            match_bases += int(num)
                        num = ''
            qSize = infer_query_size(qName, seq, fallback=match_bases)
            # end from ref-advancing ops
            ref_advance = 0
            num = ''
            for ch in cigar:
                if ch.isdigit():
                    num += ch
                else:
                    val = int(num) if num else 0
                    if ch in ('M', 'D', 'N', '=', 'X'):
                        ref_advance += val
                    num = ''
            end = pos + max(0, ref_advance) - 1
            score = match_bases
            as_score = None
            for opt in parts[11:]:
                if opt.startswith("AS:i:"):
                    try:
                        as_score = int(opt[5:])
                    except ValueError:
                        pass
                    break
            align_score = as_score if as_score is not None else score
            if qName not in results:
                continue
            aln = AlignmentHit(rName, pos, end, score, qSize, match_bases, orientation, align_score)
            results[qName].append(aln)
    except Exception as e:
        logging.error(f"SAM parse error ({source}): {e}")


class AlignmentHit:
    """Generic alignment hit container."""
    def __init__(self, rname, pos, end, score, query_size, match_size, orientation='+', align_score=None):
        self.rname = rname
        self.pos = pos
        self.end = end
        self.score = score
        self.query_size = query_size
        self.match_size = match_size
        self.orientation = orientation
        self.align_score = align_score if align_score is not None else score


class Minimap2BatchAligner:
    """Align queries with minimap2 producing SAM (-a) for consistency with bwa."""
    def __init__(self, reference, preset="asm20", threads=1, save_output_path=None):
        self.reference = reference
        self.preset = preset
        self.threads = max(1, int(threads))
        self.binary = None
        self.save_output_path = save_output_path
        self._check()

    def _check(self):
        self.binary = which("minimap2")
        if self.binary is None:
            raise FileNotFoundError("minimap2 binary not found in PATH")

    def align_batch(self, queries):
        results = {q: [] for q, _ in queries}
        if not queries:
            return results
        fd_q, q_fa = tempfile.mkstemp(prefix="remap_mm2_", suffix=".fa")
        with os.fdopen(fd_q, 'w') as out:
            for name, seq in queries:
                out.write(f">{name}\n{seq}\n")
        try:
            cmd = [self.binary, "-x", self.preset, "-a", "--max-chain-skip", "50000", "-N", "500","-k","15","-w","5","-t", str(self.threads), self.reference, q_fa]
            logging.info(f"minimap2 path: {self.binary}")
            logging.info(f"minimap2 cmd: {shlex.join(cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                log_subprocess_failure("minimap2", proc)
                raise RuntimeError(f"minimap2 exited {proc.returncode}")
            if self.save_output_path:
                try:
                    with open(self.save_output_path, 'w') as out_raw:
                        out_raw.write(proc.stdout)
                    logging.info(f"Saved minimap2 output to {self.save_output_path}")
                except Exception as e:
                    logging.warning(f"Failed to save minimap2 output to {self.save_output_path}: {e}")
            # Parse SAM from stdout using shared parser
            parse_sam_lines(proc.stdout.splitlines(), results, source="minimap2")
        finally:
            try:
                os.unlink(q_fa)
            except OSError:
                pass
        return results

    # Removed per-aligner SAM parser; using shared parse_sam_lines

class BwaBatchAligner:
    """Align queries with bwa (mem); parse SAM to produce alignment hits.

    Optionally save raw SAM output to a file for debugging.
    """
    def __init__(self, reference, threads=1, save_output_path=None):
        self.reference = reference
        self.threads = max(1, int(threads))
        self.save_output_path = save_output_path
        self.binary = None
        self._check()

    def _check(self):
        self.binary = which("bwa")
        if self.binary is None:
            raise FileNotFoundError("bwa binary not found in PATH")

    def align_batch(self, queries):
        results = {q: [] for q, _ in queries}
        if not queries:
            return results
        fd_q, q_fa = tempfile.mkstemp(prefix="remap_bwa_", suffix=".fa")
        with os.fdopen(fd_q, 'w') as out:
            for name, seq in queries:
                out.write(f">{name}\n{seq}\n")
        try:
            cmd = [self.binary, "mem", "-a", "-t", str(self.threads), self.reference, q_fa]
            logging.info(f"bwa path: {self.binary}")
            logging.info(f"bwa cmd: {shlex.join(cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                log_subprocess_failure("bwa", proc)
                raise RuntimeError(f"bwa exited {proc.returncode}")
            # Save raw SAM output if requested
            if self.save_output_path:
                try:
                    with open(self.save_output_path, 'w') as sam_out:
                        sam_out.write(proc.stdout)
                    logging.info(f"Saved bwa SAM to {self.save_output_path}")
                except Exception as e:
                    logging.warning(f"Failed to save bwa SAM to {self.save_output_path}: {e}")
            parse_sam_lines(proc.stdout.splitlines(), results, source="bwa")
        finally:
            try:
                os.unlink(q_fa)
            except OSError:
                pass
        return results


class BlastBatchAligner:
    """Align queries with blastn using a pre-built database and optimized settings."""

    def __init__(self, blast_db, threads=1, save_output_path=None,
                 task="megablast", max_target_seqs=50, evalue=1e-10):
        if not blast_db:
            raise ValueError("blastn aligner requires --blast-db pointing to a BLAST database prefix")
        self.blast_db = blast_db
        self.threads = max(1, int(threads))
        self.save_output_path = save_output_path
        self.task = task
        self.max_target_seqs = max(1, int(max_target_seqs))
        self.evalue = evalue
        self.binary = None
        self._check()

    def _check(self):
        self.binary = which("blastn")
        if self.binary is None:
            raise FileNotFoundError("blastn binary not found in PATH")

    def align_batch(self, queries):
        results = {q: [] for q, _ in queries}
        if not queries:
            return results
        fd_q, q_fa = tempfile.mkstemp(prefix="remap_blast_", suffix=".fa")
        with os.fdopen(fd_q, 'w') as out:
            for name, seq in queries:
                out.write(f">{name}\n{seq}\n")
        try:
            outfmt = "6 qseqid sseqid sstart send qlen length bitscore"
            cmd = [
                self.binary,
                "-task", self.task,
                "-db", self.blast_db,
                "-query", q_fa,
                "-outfmt", outfmt,
                "-max_target_seqs", str(self.max_target_seqs),
                "-culling_limit", "5",
                "-dust", "no",
                "-soft_masking", "false",
                "-word_size", "11",
                "-num_threads", str(self.threads),
                "-evalue", f"{self.evalue}",
            ]
            logging.info(f"blastn path: {self.binary}")
            logging.info(f"blastn cmd: {shlex.join(cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                log_subprocess_failure("blastn", proc)
                raise RuntimeError(f"blastn exited {proc.returncode}")
            if self.save_output_path:
                try:
                    with open(self.save_output_path, 'w') as blast_out:
                        blast_out.write(proc.stdout)
                    logging.info(f"Saved blastn output to {self.save_output_path}")
                except Exception as e:
                    logging.warning(f"Failed to save blastn output to {self.save_output_path}: {e}")
            for line in proc.stdout.splitlines():
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 7:
                    continue
                qname = parts[0]
                if qname not in results:
                    continue
                rname = parts[1]
                try:
                    sstart = int(float(parts[2]))
                    send = int(float(parts[3]))
                    qlen = int(float(parts[4]))
                    match_size = int(float(parts[5]))
                    bitscore = float(parts[6])
                except ValueError:
                    continue
                orientation = '+' if send >= sstart else '-'
                pos = min(sstart, send)
                end = max(sstart, send)
                qsize = qlen if qlen else infer_query_size(qname, fallback=match_size)
                int_bitscore = int(round(bitscore))
                aln = AlignmentHit(rname, pos, end, int_bitscore, qsize, match_size, orientation, int_bitscore)
                results[qname].append(aln)
        finally:
            try:
                os.unlink(q_fa)
            except OSError:
                pass
        return results

class Remap:
    """ Class for remapping annotation """

    def __init__(self, in_vcf, reference, out_vcf="/dev/stdout", min_length=50,
                 cov_threshold=0.8, aligner="bwa",
                 mm2_preset="asm20", threads=4,
                 mm2_threshold=5000, max_length=10_000_000,
                 save_alignments_prefix=None, blast_db=None,
                 param_state=None):
        self.in_vcf = in_vcf
        self.reference = reference
        self.out_vcf = out_vcf
        self.min_length = min_length
        self.cov_threshold = cov_threshold
        self.aligner_choice = aligner
        self.mm2_preset = mm2_preset
        self.threads = threads
        self.n_header = None
        self.mm2_threshold = int(mm2_threshold)
        self.max_length = int(max_length)
        self.save_alignments_prefix = save_alignments_prefix
        self._small_out_path = None
        self._large_out_path = None
        self.blast_db = blast_db
        self.param_state = param_state
        if self.save_alignments_prefix:
            prefix = self.save_alignments_prefix
            self._small_out_path = f"{prefix}small.sam"
            self._large_out_path = f"{prefix}large.sam"
            for path in (self._small_out_path, self._large_out_path):
                directory = os.path.dirname(path)
                if not directory:
                    continue
                try:
                    os.makedirs(directory, exist_ok=True)
                except Exception as e:
                    logging.warning(f"Failed to create directory for {path}: {e}")

        # Always have a minimap2 aligner available for large insertions
        self.mm2_aligner = Minimap2BatchAligner(self.reference, preset=self.mm2_preset, threads=self.threads, save_output_path=self._large_out_path)
        # Aligner used for small insertions
        if self.aligner_choice == "minimap2":
            self.small_aligner = Minimap2BatchAligner(self.reference, preset=self.mm2_preset, threads=self.threads, save_output_path=self._small_out_path)
        elif self.aligner_choice == "bwa":
            self.small_aligner = BwaBatchAligner(self.reference, threads=self.threads, save_output_path=self._small_out_path)
        elif self.aligner_choice == "blastn":
            self.small_aligner = BlastBatchAligner(blast_db=self.blast_db, threads=self.threads, save_output_path=self._small_out_path)
        else:
            raise ValueError("Unsupported aligner. Choose 'minimap2', 'bwa', or 'blastn'.")

        self._batch_results = {}
        self._failed_queries = set()  # infrastructure failures only
        self._over_max_queries = set()

    def edit_header(self, header=None):
        if header is None:
            fh = truvari.VariantFile(self.in_vcf, 'r')
            header = fh.header.copy()
        header.add_line(('##INFO=<ID=remap_classification,Number=1,Type=String,'
                         'Description="Classification of alt-seq remapping (tandem/interspersed/etc.)">'))
        header.add_line(('##INFO=<ID=remap_coords,Number=1,Type=String,'
                         'Description="Best remap hit location as chr:start-end">'))
        header.add_line(('##INFO=<ID=remap_ori,Number=1,Type=String,'
                         'Description="Orientation of best remap hit relative to reference (+/-)">'))
        header.add_line(('##INFO=<ID=remap_perc,Number=1,Type=Float,'
                         'Description="Percent of allele aligned in best remap hit">'))
        if self.param_state:
            header.add_line(f"##truvari_remapParams={self.param_state}")
        self.n_header = header

    def get_end_and_coverage(self, aln):  # pylint: disable=no-self-use
        end = aln.end
        soft_bases = aln.query_size - aln.match_size
        return end, soft_bases

    def make_qname(self, entry, seq):
        """Create a unique query name for an entry based on chrom, pos, length, and sequence hash."""
        try:
            s = str(seq)
        except Exception:
            s = seq
        digest = hashlib.sha1(s.encode('utf-8')).hexdigest()[:12]
        return f"q{entry.chrom}_{entry.pos}_{len(s)}_{digest}"

    def _align_simple(self, seq, chrom=None, pos=None, qname=None):
        """Return list of AlignmentHit for given query name after batch alignment."""
        if seq is None:
            return []
        seq = str(seq)
        if qname is None:
            return []
        return self._batch_results.get(qname, [])

    def remap_entry(self, entry, cov_threshold=None):
        # Only insertions are supported for remapping/annotation
        seq = entry.alts[0]
        qname = self.make_qname(entry, seq)
        if cov_threshold is None:
            cov_threshold = self.cov_threshold
        if qname in self._over_max_queries:
            return "over_max_size", []
        if qname in self._failed_queries:
            return "failed", []

        hit_records = []
        num_hits = 0
        partial_hits = 0
        close_dist = None
        inv_close_dist = None
        best_partial_hit = None

        alignments = self._align_simple(seq, entry.chrom, entry.pos, qname=qname)
        if not alignments:
            return "novel", hit_records

        for aln in alignments:
            end, soft = self.get_end_and_coverage(aln)
            dist = min(abs(aln.pos - entry.pos), abs(end - entry.pos))
            seq_len = len(seq)
            aligned_bases = seq_len - soft
            pct_query = aligned_bases / seq_len
            passes_threshold = pct_query >= cov_threshold
            hit_coords = f"{aln.rname}:{aln.pos}-{end}"
            chrom_dist = float('inf') if aln.rname != entry.chrom else dist
            hit_meta = {
                "coords": hit_coords,
                "orientation": aln.orientation,
                "perc": pct_query * 100.0,
                "chrom_dist": chrom_dist,
                "score": aln.align_score,
                "summary": f"{hit_coords}|{int(pct_query*100)}"
            }
            if not passes_threshold:
                partial_hits += 1
                if (best_partial_hit is None or
                        hit_meta["score"] > best_partial_hit["score"] or
                        (hit_meta["score"] == best_partial_hit["score"] and hit_meta["perc"] > best_partial_hit["perc"]) or
                        (hit_meta["score"] == best_partial_hit["score"] and hit_meta["perc"] == best_partial_hit["perc"] and chrom_dist < best_partial_hit["chrom_dist"])):
                    best_partial_hit = hit_meta
                continue
            hit_records.append(hit_meta)
            num_hits += 1
            if aln.rname == entry.chrom:
                if aln.orientation == '+':
                    if close_dist is None or dist < close_dist:
                        close_dist = dist
                else:
                    if inv_close_dist is None or dist < inv_close_dist:
                        inv_close_dist = dist

        if not hit_records and best_partial_hit:
            hit_records = [best_partial_hit]

        if num_hits == 0 and partial_hits == 0:
            classification = "novel"
        elif close_dist is not None and close_dist <= len(seq):
            classification = "tandem"
        elif inv_close_dist is not None and inv_close_dist <= len(seq):
            classification = "tandem_inverted"
        elif num_hits == 0 and partial_hits != 0:
            classification = "partial"
        else:
            classification = "interspersed"

        if hit_records:
            if classification in ("tandem", "tandem_inverted"):
                hit_records.sort(key=lambda h: (h["chrom_dist"], -h["perc"], -h["score"]))
            else:
                hit_records.sort(key=lambda h: (-h["score"], -h["perc"], h["chrom_dist"]))

        return classification, hit_records

    def annotate_entry(self, entry):
        # Only annotate insertions; all other SVs are written unchanged
        if entry.var_type() == truvari.SV.INS and entry.var_size() >= self.min_length:
            entry.translate(self.n_header)
            remap, hits = self.remap_entry(entry)
            entry.info["remap_classification"] = remap
            if hits:
                best_hit = hits[0]
                entry.info["remap_coords"] = best_hit["coords"]
                entry.info["remap_ori"] = best_hit["orientation"]
                entry.info["remap_perc"] = round(best_hit["perc"], 2)
            else:
                for key in ("remap_coords", "remap_ori", "remap_perc"):
                    try:
                        entry.info.pop(key, None)
                    except AttributeError:
                        if key in entry.info:
                            try:
                                del entry.info[key]
                            except Exception:  # pragma: no cover - best effort cleanup
                                pass
        return entry

    def annotate_vcf(self):
        fh = truvari.VariantFile(self.in_vcf)
        self.edit_header(fh.header.copy())
        out = truvari.VariantFile(self.out_vcf, 'w', header=self.n_header)
        queries = []
        small_queries = []
        mm2_queries = []
        entries = []
        logging.info(f"Parsing VCF {self.in_vcf} for insertion sequences to remap")
        for entry in fh:
            entries.append(entry)
            # Only build alignment queries for insertions
            if entry.var_type() != truvari.SV.INS:
                continue
            if entry.var_size() < self.min_length:
                continue
            seq = entry.alts[0]
            if seq.startswith("<") and seq.endswith(">"):
                continue
            seq_len = len(seq)
            if seq_len < self.min_length:
                continue
            qname = self.make_qname(entry, seq)
            if seq_len > self.max_length:
                logging.warning(
                    "Skipping %s:%d insertion (%d bp) exceeding max_length=%d",
                    entry.chrom, entry.pos, seq_len, self.max_length)
                self._over_max_queries.add(qname)
                continue
            queries.append((qname, seq))
            # Insertions larger than threshold go to minimap2, others to user's choice
            if seq_len > self.mm2_threshold:
                mm2_queries.append((qname, seq))
            else:
                small_queries.append((qname, seq))
        logging.info(f"Total sequences: {len(queries)} | small INS: {len(small_queries)} | large INS (mm2): {len(mm2_queries)}")
        self._batch_results = {}
        if small_queries:
            logging.info(f"Running {self.aligner_choice} on {len(small_queries)} sequences")
            try:
                small_res = self.small_aligner.align_batch(small_queries)
            except Exception as e:
                logging.error(f"{self.aligner_choice} aligner failure: {e}")
                for qname, _ in small_queries:
                    self._failed_queries.add(qname)
            else:
                self._batch_results.update(small_res)
        if mm2_queries:
            logging.info(f"Running minimap2 on {len(mm2_queries)} large insertions (>{self.mm2_threshold} bp)")
            try:
                mm2_res = self.mm2_aligner.align_batch(mm2_queries)
            except Exception as e:
                logging.error(f"minimap2 aligner failure: {e}")
                for qname, _ in mm2_queries:
                    self._failed_queries.add(qname)
            else:
                # Merge, mm2 results take precedence for those qnames
                self._batch_results.update(mm2_res)
        for entry in entries:
            entry = self.annotate_entry(entry)
            out.write(entry)


def parse_args(args):
    parser = argparse.ArgumentParser(prog="remap", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", type=str, default="/dev/stdin",
                        help="Input VCF (%(default)s)")
    parser.add_argument("-r", "--reference", required=True,
                        help="Reference FASTA file")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output VCF (%(default)s)")
    parser.add_argument("-m", "--min-length", default=50, type=truvari.restricted_int,
                        help="Smallest length of allele to remap (%(default)s)")
    parser.add_argument("-M", "--max-length", type=truvari.restricted_int, default=10_000_000,
                        help="Largest ALT sequence length (bp) to attempt remap (%(default)s)")
    parser.add_argument("--mm2-threshold", type=truvari.restricted_int, default=5000,
                        help="Insertions larger than this (bp) are aligned with minimap2 (%(default)s)")
    parser.add_argument("--cov-threshold", type=truvari.restricted_float, default=.8,
                        help="Threshold for pct of allele covered to consider hit (%(default)s)")
    parser.add_argument("--aligner", choices=["minimap2", "bwa", "blastn"], default="bwa",
                        help="Aligner choice (%(default)s)")
    parser.add_argument("--mm2-preset", default="asm20",
                        help="minimap2 preset (asm20 default; asm5/asm10/map-ont/sr etc.) (%(default)s)")
    parser.add_argument("--threads", type=truvari.restricted_int, default=1,
                        help="Threads for the aligner (%(default)s)")
    parser.add_argument("--blast-db", default=None,
                        help="Pre-built BLAST database prefix (required when --aligner blastn)")
    parser.add_argument("--save-alignments-prefix", default=None,
                        help="Path prefix; files named <prefix>small.sam and <prefix>large.sam will be written")
    a = parser.parse_args(args)
    truvari.setup_logging(True, show_version=True)
    return a

def remap_main(cmdargs):
    args = parse_args(cmdargs)
    try:
        param_state = json.dumps(vars(args), sort_keys=True)
    except TypeError:
        # Fallback to string representation if JSON serialization fails
        param_state = str(vars(args))
    try:
        anno = Remap(in_vcf=args.input,
             reference=args.reference,
             out_vcf=args.output,
             min_length=args.min_length,
             cov_threshold=args.cov_threshold,
             aligner=args.aligner,
             mm2_preset=args.mm2_preset,
             threads=args.threads,
             mm2_threshold=args.mm2_threshold,
             max_length=args.max_length,
             save_alignments_prefix=args.save_alignments_prefix,
             blast_db=args.blast_db,
             param_state=param_state)
        anno.annotate_vcf()
        logging.info("Finished remap")
    except Exception as e:
        logging.error(f"Error during remap: {e}")
        sys.exit(1)
