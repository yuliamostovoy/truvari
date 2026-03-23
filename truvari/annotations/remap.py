"""Remap VCF alleles' sequence to the reference genome to annotate REMAP.

Classification:
    novel        : Allele has no hits in reference
    tandem       : Allele's closest hit is within len(allele) bp of the SV's position and same orientation
    tandem_inverted : Tandem-proximal hit but inverted orientation
    tandem_complex: Local remap can be explained by multiple stitched segments near the source locus
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
        logging.error(f"{tool_name} Killed (signal {abs(rc)})")
    elif rc == 137:
        logging.error(f"{tool_name} Killed (exit {rc})")
    else:
        logging.error(f"{tool_name} exited with code {rc}")
    stderr = (proc.stderr or "").strip()
    if stderr:
        logging.error(f"{tool_name} stderr:\n{stderr}")
        lower = stderr.lower()
        if "out of memory" in lower or "oom" in lower:
            logging.error(f"{tool_name} OutOfMemory detected")

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


def parse_cigar_metrics(cigar):
    """Return alignment metrics parsed from a CIGAR string."""
    match_bases = 0
    ref_advance = 0
    query_advance = 0
    left_soft = 0
    right_soft = 0
    ops = []
    num = ''
    for ch in cigar:
        if ch.isdigit():
            num += ch
            continue
        val = int(num) if num else 0
        ops.append((val, ch))
        if ch in ('M', '=', 'X'):
            match_bases += val
        if ch in ('M', 'D', 'N', '=', 'X'):
            ref_advance += val
        if ch in ('M', 'I', 'S', '=', 'X'):
            query_advance += val
        num = ''
    if ops:
        if ops[0][1] == 'S':
            left_soft = ops[0][0]
        if ops[-1][1] == 'S':
            right_soft = ops[-1][0]
    return {
        "match_bases": match_bases,
        "ref_advance": ref_advance,
        "query_advance": query_advance,
        "left_soft": left_soft,
        "right_soft": right_soft,
    }


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
            metrics = parse_cigar_metrics(cigar)
            match_bases = metrics["match_bases"]
            if match_bases == 0:
                try:
                    match_bases = cigmatch(cigar)
                except Exception:
                    pass
            qSize = infer_query_size(qName, seq, fallback=match_bases)
            ref_advance = metrics["ref_advance"]
            end = pos + max(0, ref_advance) - 1
            query_advance = metrics["query_advance"]
            left_soft = metrics["left_soft"]
            right_soft = metrics["right_soft"]
            if orientation == '+':
                query_start = left_soft + 1
                query_end = max(left_soft, query_advance - right_soft)
            else:
                query_start = right_soft + 1
                query_end = max(right_soft, query_advance - left_soft)
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
            aln = AlignmentHit(rName, pos, end, score, qSize, match_bases, orientation, align_score,
                               query_start=query_start, query_end=query_end, cigar=cigar)
            results[qName].append(aln)
    except Exception as e:
        logging.error(f"SAM parse error ({source}): {e}")


class AlignmentHit:
    """Generic alignment hit container."""
    def __init__(self, rname, pos, end, score, query_size, match_size, orientation='+', align_score=None,
                 query_start=None, query_end=None, cigar=None):
        self.rname = rname
        self.pos = pos
        self.end = end
        self.score = score
        self.query_size = query_size
        self.match_size = match_size
        self.orientation = orientation
        self.align_score = align_score if align_score is not None else score
        self.query_start = query_start if query_start is not None else 1
        self.query_end = query_end if query_end is not None else match_size
        self.cigar = cigar


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
            cmd = [self.binary, "-x", self.preset, "-a", "-Y", "--max-chain-skip", "50000", "-N", "500","-k","15","-w","5","-t", str(self.threads), self.reference, q_fa]
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
            outfmt = "6 qseqid sseqid sstart send qstart qend qlen length bitscore"
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
                if len(parts) < 9:
                    continue
                qname = parts[0]
                if qname not in results:
                    continue
                rname = parts[1]
                try:
                    sstart = int(float(parts[2]))
                    send = int(float(parts[3]))
                    qstart = int(float(parts[4]))
                    qend = int(float(parts[5]))
                    qlen = int(float(parts[6]))
                    match_size = int(float(parts[7]))
                    bitscore = float(parts[8])
                except ValueError:
                    continue
                orientation = '+' if send >= sstart else '-'
                pos = min(sstart, send)
                end = max(sstart, send)
                qsize = qlen if qlen else infer_query_size(qname, fallback=match_size)
                int_bitscore = int(round(bitscore))
                aln = AlignmentHit(rname, pos, end, int_bitscore, qsize, match_size, orientation, int_bitscore,
                                   query_start=min(qstart, qend), query_end=max(qstart, qend))
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
                 local_rescue=True, local_window_mult=2,
                 local_window_min=2000, local_window_max=50000,
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
        self._alignments_dir = None
        self._small_out_path = None
        self._large_out_path = None
        self._local_debug_dir = None
        self.blast_db = blast_db
        self.local_rescue = bool(local_rescue)
        self.local_window_mult = max(1, int(local_window_mult))
        self.local_window_min = max(100, int(local_window_min))
        self.local_window_max = max(self.local_window_min, int(local_window_max))
        self.param_state = param_state
        if self.save_alignments_prefix:
            prefix = self.save_alignments_prefix
            self._alignments_dir = prefix if prefix.endswith(os.sep) else f"{prefix}alignments"
            self._small_out_path = os.path.join(self._alignments_dir, "small.sam")
            self._large_out_path = os.path.join(self._alignments_dir, "large.sam")
            self._local_debug_dir = os.path.join(self._alignments_dir, "local")
            for directory in (self._alignments_dir, self._local_debug_dir):
                try:
                    os.makedirs(directory, exist_ok=True)
                except Exception as e:
                    logging.warning(f"Failed to create directory {directory}: {e}")

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

        self._samtools = None
        self._contig_lengths = {}
        if self.local_rescue:
            self._samtools = which("samtools")
            fai = f"{self.reference}.fai"
            if self._samtools is None:
                logging.warning("Local rescue disabled because samtools is not in PATH")
                self.local_rescue = False
            elif not os.path.exists(fai):
                logging.warning("Local rescue disabled because FASTA index is missing: %s", fai)
                self.local_rescue = False
            else:
                with open(fai, 'r', encoding='utf-8') as fin:
                    for line in fin:
                        if not line.strip():
                            continue
                        fields = line.rstrip().split('\t')
                        if len(fields) < 2:
                            continue
                        try:
                            self._contig_lengths[fields[0]] = int(fields[1])
                        except ValueError:
                            continue
                if not self._contig_lengths:
                    logging.warning("Local rescue disabled because no contig lengths were parsed from %s", fai)
                    self.local_rescue = False

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
        header.add_line(('##INFO=<ID=remap_segments,Number=1,Type=String,'
                         'Description="Detailed stitched remap segments as chr:start-end(ori):qstart-qend joined by |">'))
        header.add_line(('##INFO=<ID=remap_local_rescue,Number=0,Type=Flag,'
                         'Description="Best remap hit or stitched local interpretation was recovered by local rescue">'))
        if self.param_state:
            header.add_line(f"##truvari_remapParams={self.param_state}")
        self.n_header = header

    def get_end_and_coverage(self, aln):  # pylint: disable=no-self-use
        end = aln.end
        aligned_bases = max(0, aln.query_end - aln.query_start + 1)
        soft_bases = max(0, aln.query_size - aligned_bases)
        return end, soft_bases

    def make_qname(self, entry, seq):
        """Create a unique query name for an entry based on chrom, pos, length, and sequence hash."""
        try:
            s = str(seq)
        except Exception:
            s = seq
        digest = hashlib.sha1(s.encode('utf-8')).hexdigest()[:12]
        return f"q{entry.chrom}_{entry.pos}_{len(s)}_{digest}"

    def _aligned_query_bases(self, aln):
        return max(0, aln.query_end - aln.query_start + 1)

    def _complement_base(self, base):
        return {
            'A': 'T',
            'C': 'G',
            'G': 'C',
            'T': 'A',
            'N': 'N',
        }.get(base.upper(), 'N')

    def _build_alignment_profile(self, aln, query_seq, ref_seq, ref_offset=0):
        profile = {}
        cigar = getattr(aln, "cigar", None)
        if not cigar:
            return profile

        ref_pos = aln.pos
        query_pos = aln.query_start if aln.orientation == '+' else aln.query_end
        query_step = 1 if aln.orientation == '+' else -1
        num = ''
        for ch in cigar:
            if ch.isdigit():
                num += ch
                continue
            length = int(num) if num else 0
            num = ''
            if ch in ('M', '=', 'X'):
                for _ in range(length):
                    query_base = query_seq[query_pos - 1].upper() if 1 <= query_pos <= len(query_seq) else 'N'
                    if aln.orientation == '-':
                        query_base = self._complement_base(query_base)
                    ref_base = ref_seq[ref_pos - 1].upper() if 1 <= ref_pos <= len(ref_seq) else 'N'
                    if ch == '=':
                        mismatch = 0
                    elif ch == 'X':
                        mismatch = 1
                    else:
                        mismatch = int(query_base != ref_base)
                    profile[query_pos] = {
                        "mismatch": mismatch,
                        "ref_pos": ref_offset + ref_pos - 1,
                    }
                    query_pos += query_step
                    ref_pos += 1
            elif ch == 'I':
                for _ in range(length):
                    profile[query_pos] = {
                        "mismatch": 1,
                        "ref_pos": None,
                    }
                    query_pos += query_step
            elif ch in ('D', 'N'):
                ref_pos += length
            elif ch in ('S', 'H', 'P'):
                # `query_start`/`query_end` already describe the aligned interval on the
                # original query coordinate system. Advancing `query_pos` across clipped
                # bases here double-counts clipping and breaks reverse-strand profiles.
                continue
        return profile

    def _build_hit_meta(self, aln, entry, seq_len, hit_pos=None, hit_end=None, local_rescue=False):
        pos = aln.pos if hit_pos is None else hit_pos
        end = aln.end if hit_end is None else hit_end
        dist = min(abs(pos - entry.pos), abs(end - entry.pos))
        aligned_bases = self._aligned_query_bases(aln)
        pct_query = (aligned_bases / seq_len) if seq_len else 0.0
        hit_coords = f"{aln.rname}:{pos}-{end}"
        chrom_dist = float('inf') if aln.rname != entry.chrom else dist
        hit_meta = {
            "coords": hit_coords,
            "rname": aln.rname,
            "start": pos,
            "end": end,
            "orientation": aln.orientation,
            "perc": pct_query * 100.0,
            "chrom_dist": chrom_dist,
            "score": aln.align_score,
            "summary": f"{hit_coords}|{int(pct_query * 100)}",
            "query_start": aln.query_start,
            "query_end": aln.query_end,
            "query_size": seq_len,
            "aligned_bases": aligned_bases,
            "local_rescue": local_rescue,
        }
        return hit_meta

    def _segment_to_string(self, seg):
        return f"{seg['coords']}({seg['orientation']}):q{seg['query_start']}-{seg['query_end']}"

    def _debug_log(self, debug_lines, message):
        if debug_lines is not None:
            debug_lines.append(str(message))

    def _write_local_complex_debug(self, qname, debug_lines):
        if not self._local_debug_dir or not debug_lines:
            return
        base_name = self._safe_debug_name(qname)
        debug_path = os.path.join(self._local_debug_dir, f"{base_name}.complex_debug.txt")
        try:
            with open(debug_path, 'w', encoding='utf-8') as fout:
                fout.write("\n".join(debug_lines))
                fout.write("\n")
        except Exception as exc:  # pragma: no cover - debug best effort
            logging.warning("Failed to save local complex debug output to %s: %s", debug_path, exc)

    def _trim_segment_interval(self, seg, keep_start, keep_end, entry, debug_lines=None, label=None):
        if keep_start > keep_end:
            self._debug_log(debug_lines, f"trim {label or '.'}: dropped because keep_start {keep_start} > keep_end {keep_end} for {self._segment_to_string(seg)}")
            return None
        profile = seg.get("profile", {})
        kept_query_positions = [pos for pos in sorted(profile) if keep_start <= pos <= keep_end]
        if not kept_query_positions:
            self._debug_log(debug_lines, f"trim {label or '.'}: dropped because no kept query positions in {keep_start}-{keep_end} for {self._segment_to_string(seg)}")
            return None
        ref_positions = [profile[pos]["ref_pos"] for pos in kept_query_positions if profile[pos].get("ref_pos") is not None]
        if not ref_positions:
            self._debug_log(debug_lines, f"trim {label or '.'}: dropped because no reference positions remain in {keep_start}-{keep_end} for {self._segment_to_string(seg)}")
            return None
        trimmed = dict(seg)
        trimmed["query_start"] = keep_start
        trimmed["query_end"] = keep_end
        trimmed["aligned_bases"] = keep_end - keep_start + 1
        trimmed["start"] = min(ref_positions)
        trimmed["end"] = max(ref_positions)
        trimmed["coords"] = f"{trimmed['rname']}:{trimmed['start']}-{trimmed['end']}"
        trimmed["chrom_dist"] = float('inf') if trimmed["rname"] != entry.chrom else min(abs(trimmed["start"] - entry.pos), abs(trimmed["end"] - entry.pos))
        trimmed["perc"] = (trimmed["aligned_bases"] / trimmed["query_size"]) * 100.0 if trimmed["query_size"] else 0.0
        original_span = max(1, seg["query_end"] - seg["query_start"] + 1)
        trimmed["score"] = seg["score"] * (trimmed["aligned_bases"] / original_span)
        trimmed["summary"] = self._segment_to_string(trimmed)
        self._debug_log(debug_lines, f"trim {label or '.'}: {self._segment_to_string(seg)} -> {self._segment_to_string(trimmed)}")
        return trimmed

    def _pick_overlap_cut(self, left_seg, right_seg, overlap_start, overlap_end, debug_lines=None):
        overlap_positions = list(range(overlap_start, overlap_end + 1))
        left_profile = left_seg.get("profile", {})
        right_profile = right_seg.get("profile", {})
        left_costs = [left_profile.get(pos, {"mismatch": 1})["mismatch"] for pos in overlap_positions]
        right_costs = [right_profile.get(pos, {"mismatch": 1})["mismatch"] for pos in overlap_positions]
        n_overlap = len(overlap_positions)

        left_prefix = [0] * (n_overlap + 1)
        for idx in range(n_overlap):
            left_prefix[idx + 1] = left_prefix[idx] + left_costs[idx]
        right_suffix = [0] * (n_overlap + 1)
        for idx in range(n_overlap - 1, -1, -1):
            right_suffix[idx] = right_suffix[idx + 1] + right_costs[idx]

        candidate_cuts = []
        best_cost = None
        for keep_left_count in range(n_overlap + 1):
            total_cost = left_prefix[keep_left_count] + right_suffix[keep_left_count]
            if best_cost is None or total_cost < best_cost:
                best_cost = total_cost
                candidate_cuts = [keep_left_count]
            elif total_cost == best_cost:
                candidate_cuts.append(keep_left_count)

        prefer_left = (
            left_seg["aligned_bases"] > right_seg["aligned_bases"] or
            (left_seg["aligned_bases"] == right_seg["aligned_bases"] and left_seg["score"] >= right_seg["score"])
        )
        chosen_cut = max(candidate_cuts) if prefer_left else min(candidate_cuts)
        self._debug_log(
            debug_lines,
            "overlap cut: "
            f"left={self._segment_to_string(left_seg)} right={self._segment_to_string(right_seg)} "
            f"overlap=q{overlap_start}-{overlap_end} best_cost={best_cost} candidate_counts={candidate_cuts} "
            f"prefer_left={prefer_left} chosen_count={chosen_cut} chosen_q={overlap_start + chosen_cut - 1}"
        )
        return overlap_start + chosen_cut - 1

    def _resolve_pair_overlap(self, left_seg, right_seg, entry, debug_lines=None):
        overlap_start = max(left_seg["query_start"], right_seg["query_start"])
        overlap_end = min(left_seg["query_end"], right_seg["query_end"])
        if overlap_start > overlap_end:
            self._debug_log(debug_lines, f"no overlap: {self._segment_to_string(left_seg)} vs {self._segment_to_string(right_seg)}")
            return left_seg, right_seg

        cut_point = self._pick_overlap_cut(left_seg, right_seg, overlap_start, overlap_end, debug_lines=debug_lines)
        left_keep_end = min(left_seg["query_end"], cut_point)
        right_keep_start = max(right_seg["query_start"], cut_point + 1)

        self._debug_log(
            debug_lines,
            f"resolve overlap: left_keep={left_seg['query_start']}-{left_keep_end} right_keep={right_keep_start}-{right_seg['query_end']}"
        )
        trimmed_left = self._trim_segment_interval(left_seg, left_seg["query_start"], left_keep_end, entry, debug_lines=debug_lines, label="left")
        trimmed_right = self._trim_segment_interval(right_seg, right_keep_start, right_seg["query_end"], entry, debug_lines=debug_lines, label="right")
        return trimmed_left, trimmed_right

    def _resolve_chain_overlaps(self, chain, entry, debug_lines=None):
        if not chain:
            return []
        self._debug_log(debug_lines, "chain before overlap resolution: " + "|".join(self._segment_to_string(seg) for seg in sorted(chain, key=lambda s: (s["query_start"], s["query_end"]))))
        resolved = []
        for seg in sorted(chain, key=lambda s: (s["query_start"], s["query_end"])):
            current = seg
            while resolved and current is not None and resolved[-1]["query_end"] >= current["query_start"]:
                prior = resolved.pop()
                prior, current = self._resolve_pair_overlap(prior, current, entry, debug_lines=debug_lines)
                if prior is not None:
                    resolved.append(prior)
            if current is not None:
                resolved.append(current)
        self._debug_log(debug_lines, "chain after overlap resolution: " + "|".join(self._segment_to_string(seg) for seg in resolved))
        return resolved

    def _drop_nested_segments(self, segments, max_query_overlap):
        kept = []
        for seg in sorted(segments, key=lambda s: (s["query_start"], -(s["query_end"] - s["query_start"]), -s["score"])):
            nested = False
            for prev in kept:
                if seg["query_start"] >= prev["query_start"] and seg["query_end"] <= prev["query_end"]:
                    nested = True
                    break
                overlap = min(seg["query_end"], prev["query_end"]) - max(seg["query_start"], prev["query_start"]) + 1
                if overlap > max_query_overlap and seg["aligned_bases"] <= prev["aligned_bases"] and seg["score"] <= prev["score"]:
                    nested = True
                    break
            if not nested:
                kept.append(seg)
        return kept

    def _chain_query_segments(self, segments, seq_len, debug_lines=None):
        if not segments:
            return []
        max_query_overlap = max(25, int(seq_len * 0.05))
        max_query_gap = max(50, int(seq_len * 0.25))
        cleaned = self._drop_nested_segments(segments, max_query_overlap)
        self._debug_log(debug_lines, "candidate segments: " + "|".join(self._segment_to_string(seg) for seg in sorted(segments, key=lambda s: (s["query_start"], s["query_end"]))))
        self._debug_log(debug_lines, "cleaned segments: " + "|".join(self._segment_to_string(seg) for seg in sorted(cleaned, key=lambda s: (s["query_start"], s["query_end"]))))
        cleaned.sort(key=lambda s: (s["query_start"], s["query_end"], -s["score"]))
        best_score = [0.0] * len(cleaned)
        best_chain = [[] for _ in cleaned]
        for idx, seg in enumerate(cleaned):
            cov = float(seg["aligned_bases"])
            best_score[idx] = cov
            best_chain[idx] = [seg]
            self._debug_log(
                debug_lines,
                f"dp init idx={idx} seg={self._segment_to_string(seg)} cov={cov} best_score={best_score[idx]:.6f}"
            )
            for prev_idx in range(idx):
                prev = cleaned[prev_idx]
                if seg["query_start"] > prev["query_end"] + max_query_gap:
                    self._debug_log(
                        debug_lines,
                        f"dp skip idx={idx} prev_idx={prev_idx} reason=gap seg={self._segment_to_string(seg)} prev={self._segment_to_string(prev)} max_query_gap={max_query_gap}"
                    )
                    continue
                new_bases = seg["query_end"] - max(seg["query_start"], prev["query_end"] + 1) + 1
                if new_bases <= 0:
                    if (prev["query_end"] - seg["query_start"] + 1) > max_query_overlap:
                        self._debug_log(
                            debug_lines,
                            f"dp skip idx={idx} prev_idx={prev_idx} reason=overlap seg={self._segment_to_string(seg)} prev={self._segment_to_string(prev)} overlap={prev['query_end'] - seg['query_start'] + 1} max_query_overlap={max_query_overlap}"
                        )
                        continue
                    new_bases = 0
                candidate_score = best_score[prev_idx] + max(0, new_bases) + (seg["score"] / 1000000.0)
                self._debug_log(
                    debug_lines,
                    f"dp consider idx={idx} prev_idx={prev_idx} seg={self._segment_to_string(seg)} prev={self._segment_to_string(prev)} new_bases={new_bases} prev_score={best_score[prev_idx]:.6f} seg_bonus={(seg['score'] / 1000000.0):.6f} candidate_score={candidate_score:.6f} current_best={best_score[idx]:.6f}"
                )
                if candidate_score > best_score[idx]:
                    best_score[idx] = candidate_score
                    best_chain[idx] = best_chain[prev_idx] + [seg]
                    self._debug_log(
                        debug_lines,
                        f"dp update idx={idx} prev_idx={prev_idx} new_best={best_score[idx]:.6f} chain={'|'.join(self._segment_to_string(x) for x in best_chain[idx])}"
                    )
            self._debug_log(
                debug_lines,
                f"dp final idx={idx} best_score={best_score[idx]:.6f} chain={'|'.join(self._segment_to_string(x) for x in best_chain[idx])}"
            )
            chosen = max(best_chain, key=lambda chain: (self._query_coverage(chain), sum(seg["score"] for seg in chain), len(chain)))
        for idx, chain in enumerate(best_chain):
            self._debug_log(
                debug_lines,
                f"dp candidate_chain idx={idx} coverage={self._query_coverage(chain)} score_sum={sum(seg['score'] for seg in chain)} length={len(chain)} chain={'|'.join(self._segment_to_string(x) for x in chain)}"
            )
            self._debug_log(debug_lines, "chosen chain before overlap resolution: " + "|".join(self._segment_to_string(seg) for seg in chosen))
            return chosen

    def _query_coverage(self, segments):
        if not segments:
            return 0
        ordered = sorted(segments, key=lambda s: (s["query_start"], s["query_end"]))
        covered = 0
        cur_start = None
        cur_end = None
        for seg in ordered:
            if cur_start is None:
                cur_start = seg["query_start"]
                cur_end = seg["query_end"]
                continue
            if seg["query_start"] <= cur_end + 1:
                cur_end = max(cur_end, seg["query_end"])
            else:
                covered += cur_end - cur_start + 1
                cur_start = seg["query_start"]
                cur_end = seg["query_end"]
        covered += cur_end - cur_start + 1
        return covered

    def _summarize_segment_chain(self, entry, chain, seq_len):
        if not chain:
            return None
        sorted_chain = sorted(chain, key=lambda s: (s["query_start"], s["query_end"]))
        overall_start = min(seg["start"] for seg in sorted_chain)
        overall_end = max(seg["end"] for seg in sorted_chain)
        covered = self._query_coverage(sorted_chain)
        perc = (covered / seq_len) * 100.0 if seq_len else 0.0
        orientation_pattern = ''.join(seg["orientation"] for seg in sorted_chain)
        return {
            "coords": f"{entry.chrom}:{overall_start}-{overall_end}",
            "orientation": orientation_pattern,
            "perc": perc,
            "chrom_dist": min(seg["chrom_dist"] for seg in sorted_chain),
            "score": sum(seg["score"] for seg in sorted_chain),
            "summary": '|'.join(self._segment_to_string(seg) for seg in sorted_chain),
            "segments": '|'.join(self._segment_to_string(seg) for seg in sorted_chain),
            "segment_count": len(sorted_chain),
            "local_rescue": True,
        }

    def _local_complex_hit(self, entry, local_segments, seq_len, cov_threshold, debug_lines=None):
        min_segment_bases = max(20, min(100, int(seq_len * 0.1)))
        candidates = [seg for seg in local_segments if seg["aligned_bases"] >= min_segment_bases]
        self._debug_log(debug_lines, f"local complex input segments ({len(local_segments)} total, {len(candidates)} candidates): " + "|".join(self._segment_to_string(seg) for seg in sorted(candidates, key=lambda s: (s['query_start'], s['query_end']))))
        if len(candidates) < 2:
            return None
        chain = self._chain_query_segments(candidates, seq_len, debug_lines=debug_lines)
        chain = self._resolve_chain_overlaps(chain, entry, debug_lines=debug_lines)
        if len(chain) < 2:
            self._debug_log(debug_lines, "complex chain dropped because fewer than 2 segments remain after overlap resolution")
            return None
        covered = self._query_coverage(chain)
        pct_query = (covered / seq_len) if seq_len else 0.0
        self._debug_log(debug_lines, f"complex chain coverage={covered} pct_query={pct_query:.4f} threshold={max(0.5, cov_threshold * 0.75):.4f}")
        if pct_query < max(0.5, cov_threshold * 0.75):
            return None
        self._debug_log(debug_lines, "final complex chain: " + "|".join(self._segment_to_string(seg) for seg in chain))
        return self._summarize_segment_chain(entry, chain, seq_len)

    def _align_simple(self, seq, chrom=None, pos=None, qname=None):
        """Return list of AlignmentHit for given query name after batch alignment."""
        if seq is None:
            return []
        seq = str(seq)
        if qname is None:
            return []
        return self._batch_results.get(qname, [])

    def _safe_debug_name(self, name):
        return ''.join(ch if ch.isalnum() or ch in '._-' else '_' for ch in str(name))

    def _write_local_debug_alignment(self, entry, region, qname, ref_name, ref_seq, query_seq, sam_text):
        if not self._local_debug_dir:
            return
        base_name = self._safe_debug_name(qname)
        sam_path = os.path.join(self._local_debug_dir, f"{base_name}.sam")
        ref_path = os.path.join(self._local_debug_dir, f"{base_name}.ref.fa")
        query_path = os.path.join(self._local_debug_dir, f"{base_name}.query.fa")
        try:
            with open(sam_path, 'w', encoding='utf-8') as fout:
                fout.write(f"# variant={entry.chrom}:{entry.pos} id={entry.id if entry.id is not None else '.'} qname={qname} region={region}\n")
                fout.write(sam_text)
                if sam_text and not sam_text.endswith('\n'):
                    fout.write('\n')
            with open(ref_path, 'w', encoding='utf-8') as fout:
                fout.write(f">{ref_name}|region={region}|variant={entry.chrom}:{entry.pos}\n{ref_seq}\n")
            with open(query_path, 'w', encoding='utf-8') as fout:
                fout.write(f">{qname}|variant={entry.chrom}:{entry.pos}|id={entry.id if entry.id is not None else '.'}\n{query_seq}\n")
            logging.info("Saved local rescue debug outputs to %s", self._local_debug_dir)
        except Exception as exc:  # pragma: no cover - debug best effort
            logging.warning("Failed to save local rescue alignment outputs to %s: %s", self._local_debug_dir, exc)

    def _local_remap_rescue(self, entry, seq, cov_threshold):
        if not self.local_rescue:
            return None
        if entry.chrom not in self._contig_lengths:
            return None

        seq_len = len(seq)
        flank = max(1, self.local_window_mult * seq_len)
        contig_len = self._contig_lengths[entry.chrom]
        start_1based = max(1, entry.pos - flank)
        end_1based = min(contig_len, entry.pos + flank)
        if end_1based <= start_1based:
            return None

        region = f"{entry.chrom}:{start_1based}-{end_1based}"
        cmd = [self._samtools, "faidx", self.reference, region]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            log_subprocess_failure("samtools faidx(local)", proc)
            return None
        ref_lines = [line.strip() for line in proc.stdout.splitlines() if line and not line.startswith(">")]
        ref_seq = "".join(ref_lines)
        if not ref_seq:
            return None

        qname = self.make_qname(entry, seq)
        debug_lines = [
            f"variant={entry.chrom}:{entry.pos}",
            f"id={entry.id if entry.id is not None else '.'}",
            f"qname={qname}",
            f"region={region}",
            f"seq_len={seq_len}",
            f"cov_threshold={cov_threshold}",
        ]
        local_ref_name = entry.chrom
        local_hits = {qname: []}
        fd_ref, local_ref_fa = tempfile.mkstemp(prefix="remap_local_ref_", suffix=".fa")
        fd_q, local_q_fa = tempfile.mkstemp(prefix="remap_local_q_", suffix=".fa")
        try:
            with os.fdopen(fd_ref, 'w') as out_ref:
                out_ref.write(f">{local_ref_name}\n{ref_seq}\n")
            with os.fdopen(fd_q, 'w') as out_q:
                out_q.write(f">{qname}\n{seq}\n")

            cmd = [
                self.mm2_aligner.binary,
                "-x", "sr",
                "-a",
                "-Y",
                "--secondary=no",
                "-N", "200",
                "-k", "11",
                "-w", "5",
                "-t", str(self.threads),
                local_ref_fa,
                local_q_fa
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                log_subprocess_failure("minimap2(local)", proc)
                return None
            self._write_local_debug_alignment(entry, region, qname, local_ref_name, ref_seq, seq, proc.stdout)
            parse_sam_lines(proc.stdout.splitlines(), local_hits, source="minimap2_local")
        except Exception as exc:  # pragma: no cover - safety net
            logging.warning("Local remap rescue failed for %s:%d: %s", entry.chrom, entry.pos, exc)
            return None
        finally:
            for path in (local_ref_fa, local_q_fa):
                try:
                    os.unlink(path)
                except OSError:
                    pass
        best_hit = None
        local_segments = []
        self._debug_log(debug_lines, f"parsed_local_hits={len(local_hits.get(qname, []))}")
        for aln in local_hits.get(qname, []):
            end, _ = self.get_end_and_coverage(aln)
            if seq_len == 0:
                continue
            global_pos = start_1based + aln.pos - 1
            global_end = start_1based + end - 1
            hit_meta = self._build_hit_meta(aln, entry, seq_len, hit_pos=global_pos, hit_end=global_end, local_rescue=True)
            hit_meta["profile"] = self._build_alignment_profile(aln, seq, ref_seq, ref_offset=start_1based)
            local_segments.append(hit_meta)
            self._debug_log(
                debug_lines,
                "local hit: "
                f"local={aln.rname}:{aln.pos}-{end} global={hit_meta['coords']} q={aln.query_start}-{aln.query_end} "
                f"ori={aln.orientation} cigar={getattr(aln, 'cigar', None)} score={aln.align_score}"
            )
            if hit_meta["perc"] / 100.0 < cov_threshold:
                continue
            if hit_meta["chrom_dist"] > seq_len:
                continue
            if (best_hit is None or
                    hit_meta["chrom_dist"] < best_hit["chrom_dist"] or
                    (hit_meta["chrom_dist"] == best_hit["chrom_dist"] and hit_meta["score"] > best_hit["score"]) or
                    (hit_meta["chrom_dist"] == best_hit["chrom_dist"] and hit_meta["score"] == best_hit["score"] and hit_meta["perc"] > best_hit["perc"])):
                best_hit = hit_meta
        if best_hit is not None:
            self._debug_log(debug_lines, f"best_simple_local_hit={self._segment_to_string(best_hit)}")
        else:
            self._debug_log(debug_lines, "best_simple_local_hit=None")

        complex_hit = self._local_complex_hit(entry, local_segments, seq_len, cov_threshold, debug_lines=debug_lines)
        if complex_hit is not None:
            self._debug_log(debug_lines, f"complex_hit={complex_hit.get('segments', complex_hit.get('coords'))}")
            if best_hit is None:
                self._write_local_complex_debug(qname, debug_lines)
                return {"classification": "tandem_complex", "best_hit": complex_hit}
            if complex_hit["perc"] > best_hit["perc"] or (complex_hit["perc"] == best_hit["perc"] and complex_hit["score"] > best_hit["score"]):
                self._write_local_complex_debug(qname, debug_lines)
                return {"classification": "tandem_complex", "best_hit": complex_hit}
        else:
            self._debug_log(debug_lines, "complex_hit=None")

        if best_hit is not None:
            classification = "tandem" if best_hit["orientation"] == '+' else "tandem_inverted"
            self._write_local_complex_debug(qname, debug_lines)
            return {"classification": classification, "best_hit": best_hit}
        self._write_local_complex_debug(qname, debug_lines)
        return best_hit

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
            seq_len = len(seq)
            hit_meta = self._build_hit_meta(aln, entry, seq_len)
            pct_query = hit_meta["perc"] / 100.0
            passes_threshold = pct_query >= cov_threshold
            if not passes_threshold:
                partial_hits += 1
                if (best_partial_hit is None or
                        hit_meta["score"] > best_partial_hit["score"] or
                        (hit_meta["score"] == best_partial_hit["score"] and hit_meta["perc"] > best_partial_hit["perc"]) or
                        (hit_meta["score"] == best_partial_hit["score"] and hit_meta["perc"] == best_partial_hit["perc"] and hit_meta["chrom_dist"] < best_partial_hit["chrom_dist"])):
                    best_partial_hit = hit_meta
                continue
            hit_records.append(hit_meta)
            num_hits += 1
            if aln.rname == entry.chrom:
                dist = hit_meta["chrom_dist"]
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

        original_classification = classification
        local_hit = self._local_remap_rescue(entry, seq, cov_threshold)
        if local_hit is not None:
            if isinstance(local_hit, dict) and "best_hit" in local_hit:
                local_classification = local_hit["classification"]
                selected_hit = local_hit["best_hit"]
            else:
                local_classification = "tandem" if local_hit["orientation"] == '+' else "tandem_inverted"
                selected_hit = local_hit

            should_apply_local = (
                original_classification not in ("tandem", "tandem_inverted")
                or local_classification == "tandem_complex"
            )
            if should_apply_local:
                classification = local_classification
                retained = [h for h in hit_records if h.get("summary") != selected_hit.get("summary")]
                hit_records = [selected_hit] + retained

        if hit_records:
            if classification in ("tandem", "tandem_inverted"):
                hit_records.sort(key=lambda h: (h["chrom_dist"], -h["perc"], -h["score"]))
            elif classification == "tandem_complex":
                hit_records.sort(key=lambda h: (-h.get("segment_count", 1), -h["perc"], -h["score"], h["chrom_dist"]))
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
                if best_hit.get("segments"):
                    entry.info["remap_segments"] = best_hit["segments"]
                else:
                    try:
                        entry.info.pop("remap_segments", None)
                    except AttributeError:
                        if "remap_segments" in entry.info:
                            try:
                                del entry.info["remap_segments"]
                            except Exception:  # pragma: no cover - best effort cleanup
                                pass
                if best_hit.get("local_rescue"):
                    entry.info["remap_local_rescue"] = True
                else:
                    try:
                        entry.info.pop("remap_local_rescue", None)
                    except AttributeError:
                        if "remap_local_rescue" in entry.info:
                            try:
                                del entry.info["remap_local_rescue"]
                            except Exception:  # pragma: no cover - best effort cleanup
                                pass
            else:
                for key in ("remap_coords", "remap_ori", "remap_perc", "remap_segments", "remap_local_rescue"):
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
                raise RuntimeError(f"{self.aligner_choice} aligner failure") from e
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
                raise RuntimeError("minimap2 aligner failure") from e
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
                        help="Path prefix or directory stem; outputs are written under <prefix>alignments/")
    parser.add_argument("--disable-local-rescue", action="store_true",
                        help="Disable local rescue remapping for non-tandem calls")
    parser.add_argument("--local-window-mult", type=truvari.restricted_int, default=2,
                        help="Local rescue window multiplier on insertion length (%(default)s)")
    parser.add_argument("--local-window-min", type=truvari.restricted_int, default=2000,
                        help="Deprecated; retained for compatibility and currently unused")
    parser.add_argument("--local-window-max", type=truvari.restricted_int, default=50000,
                        help="Deprecated; retained for compatibility and currently unused")
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
             local_rescue=not args.disable_local_rescue,
             local_window_mult=args.local_window_mult,
             local_window_min=args.local_window_min,
             local_window_max=args.local_window_max,
             param_state=param_state)
        anno.annotate_vcf()
        logging.info("Finished remap")
    except Exception as e:
        logging.error(f"Error during remap: {e}")
        sys.exit(1)
