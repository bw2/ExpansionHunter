"""Integration tests for ExpansionHunter analysis modes on simulated STR data.

Two things are checked against the simulated read sets in
``../PolymorphicTandemRepeatFinder/simulated_data``:

1. ``seeking``, ``streaming`` and ``low-mem-streaming`` modes must produce
   identical genotype calls (GT / REPCN / REPCI) for every sample. These three
   modes share the genotyping core and are expected to agree exactly.

2. The accuracy of ``seeking`` vs ``optimized-streaming`` is compared against the
   known true allele sizes. ``optimized-streaming`` uses a faster, heuristic
   genotyping path (incl. the post-call STR over-call correction) so it is NOT
   expected to match the other modes; this test quantifies how its accuracy
   differs from the exact ``seeking`` path.

Each simulated ``sim_<N>x__10_150_450_50.bam`` holds reads from a single allele
with N repeats. A diploid sample with truth N1/N2 is built by merging two such
bams; read names are prefixed per source allele so coverage doubles correctly
even when N1 == N2 (ExpansionHunter deduplicates identical read names, so a
naive self-merge would not raise coverage).

Run with:  python3 -m unittest mode_consistency_and_accuracy_tests -v
"""

import glob
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest

import pysam


REPO_DIR = os.path.dirname(os.path.abspath(__file__))
SIM_DATA_DIR = os.path.join(REPO_DIR, "..", "PolymorphicTandemRepeatFinder", "simulated_data")
EH_BINARY = os.path.join(REPO_DIR, "ehunter", "build", "ExpansionHunter")
REFERENCE_CANDIDATES = [
    "/Users/weisburd/hg38.fa",
    "/Users/weisburd/p1/ref/GRCh38/hg38.fa",
    "/Users/weisburd/code/str-truth-set/ref/hg38.fa",
    "/Users/weisburd/code/tandem-repeat-catalogs/hg38.fa",
]
MERGED_BAM_CACHE_DIR = os.path.join(REPO_DIR, ".merged_bam_cache")

# Modes expected to be byte-for-byte equivalent at the genotype level.
IDENTICAL_MODES = ["seeking", "streaming", "low-mem-streaming"]
# Modes compared for accuracy against the known truth.
ACCURACY_MODES = ["seeking", "optimized-streaming"]
ALL_MODES = ["seeking", "streaming", "low-mem-streaming", "optimized-streaming"]

# How many representative (N1, N2) pairs to test per locus.
MAX_PAIRS_PER_LOCUS = 7

BAM_NAME_RE = re.compile(r"^sim_(\d+)x__10_150_450_50\.bam$")
LOCUS_DIR_RE = re.compile(r"^(chr[0-9XYM]+)-(\d+)-(\d+)-([ACGT]+)$")


def find_reference():
    """Returns the first reference FASTA (with .fai) that exists, else None."""
    return next((p for p in REFERENCE_CANDIDATES if os.path.exists(p) and os.path.exists(p + ".fai")), None)


def parse_locus_dir(name):
    """Parses a locus directory name into (locus_id, chrom, start_0based, end_1based, motif).

    Returns None if the name does not match the expected ``chrom-start-end-motif`` form.
    """
    m = LOCUS_DIR_RE.match(name)
    if not m:
        return None
    return (name, m.group(1), int(m.group(2)), int(m.group(3)), m.group(4))


def allele_sizes_in_dir(locus_dir):
    """Returns the sorted list of repeat counts N for which a sim bam exists in locus_dir."""
    sizes = []
    for path in glob.glob(os.path.join(locus_dir, "sim_*x__10_150_450_50.bam")):
        m = BAM_NAME_RE.match(os.path.basename(path))
        if m:
            sizes.append(int(m.group(1)))
    return sorted(set(sizes))


def select_pairs(sizes):
    """Selects up to MAX_PAIRS_PER_LOCUS representative (N1, N2) diploid genotypes.

    Spans homozygous and heterozygous cases across the reference-sized, medium and
    largest available alleles so both agreement and accuracy are exercised on a
    range of expansion sizes. Each pair is ordered N1 <= N2 and the list is
    deterministic for a given input.
    """
    if len(sizes) < 2:
        return []
    lo, hi, mid = sizes[0], sizes[-1], sizes[len(sizes) // 2]
    candidates = [
        (lo, lo),        # homozygous, ~reference sized
        (lo, sizes[1]),  # small heterozygous expansion
        (lo, mid),       # medium heterozygous expansion
        (lo, hi),        # large heterozygous expansion
        (mid, mid),      # homozygous medium expansion
        (mid, hi),       # large heterozygous expansion
        (hi, hi),        # homozygous large expansion
    ]
    pairs = []
    for pair in candidates:
        if pair not in pairs:
            pairs.append(pair)
    return pairs[:MAX_PAIRS_PER_LOCUS]


def build_merged_bam(bam_n1, bam_n2, out_bam):
    """Merges two single-allele bams into a diploid bam, prefixing read names per allele.

    Read names from the first/second source are prefixed with ``a1_``/``a2_`` so that
    reads remain distinct (and coverage doubles) even when the two source bams are the
    same file (homozygous genotype). The result is coordinate-sorted and indexed.
    Skips work if out_bam already exists (pre-generate once, reuse across runs).
    """
    if os.path.exists(out_bam) and os.path.exists(out_bam + ".bai"):
        return
    os.makedirs(os.path.dirname(out_bam), exist_ok=True)
    with pysam.AlignmentFile(bam_n1, "rb") as f:
        header = f.header.to_dict()
    unsorted_bam = out_bam + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as writer:
        for source_bam, tag in ((bam_n1, "a1"), (bam_n2, "a2")):
            with pysam.AlignmentFile(source_bam, "rb") as f:
                for read in f:
                    read.query_name = tag + "_" + read.query_name
                    writer.write(read)
    pysam.sort("-o", out_bam, unsorted_bam)
    pysam.index(out_bam)
    os.remove(unsorted_bam)


def write_catalog(locus_id, chrom, start_0based, end_1based, motif, out_path):
    """Writes a single-locus ExpansionHunter variant catalog JSON."""
    with open(out_path, "w") as f:
        json.dump([{
            "LocusId": locus_id,
            "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
            "LocusStructure": f"({motif})*",
            "VariantType": "Repeat",
            "VariantId": locus_id,
        }], f)


def parse_repcn(repcn):
    """Parses a REPCN field (e.g. '17/29', '29/29', './.') into a sorted tuple of ints.

    Returns None for a no-call.
    """
    parts = repcn.split("/")
    if not all(p.isdigit() for p in parts):
        return None
    return tuple(sorted(int(p) for p in parts))


def run_expansion_hunter(reads_bam, reference, catalog, output_prefix, mode):
    """Runs ExpansionHunter and returns the parsed call dict for the single locus.

    The dict has keys: ``GT``, ``REPCN`` (raw string), ``REPCI``, and ``alleles``
    (sorted int tuple or None for a no-call). Raises CalledProcessError on a non-zero exit.
    """
    subprocess.run(
        [EH_BINARY, "--reads", reads_bam, "--reference", reference,
         "--variant-catalog", catalog, "--output-prefix", output_prefix,
         "--analysis-mode", mode, "--sort-catalog-by", "position"],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(output_prefix + ".vcf") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            return {
                "GT": sample.get("GT"),
                "REPCN": sample.get("REPCN"),
                "REPCI": sample.get("REPCI"),
                "alleles": parse_repcn(sample.get("REPCN", "./.")),
            }
    raise AssertionError(f"No VCF record produced for mode={mode} prefix={output_prefix}")


def allele_match_count(call, truth):
    """Returns the number of alleles (0-2) in `call` that match `truth` as a multiset."""
    if call is None:
        return 0
    remaining = list(truth)
    matched = 0
    for allele in call:
        if allele in remaining:
            remaining.remove(allele)
            matched += 1
    return matched


def within_one(call, truth):
    """Returns True if both alleles match truth within +/-1 repeat (multiset, greedy)."""
    if call is None:
        return False
    remaining = list(truth)
    for allele in call:
        hit = next((t for t in remaining if abs(t - allele) <= 1), None)
        if hit is None:
            return False
        remaining.remove(hit)
    return not remaining


class ModeConsistencyAndAccuracyTest(unittest.TestCase):
    """Runs every analysis mode once per simulated sample and stores the parsed calls."""

    # Populated by setUpClass: {(locus_id, n1, n2): {mode: call_dict_or_None}}
    RESULTS = {}
    SKIP_REASON = None

    @classmethod
    def setUpClass(cls):
        cls.SKIP_REASON = None
        if not os.path.isdir(SIM_DATA_DIR):
            cls.SKIP_REASON = f"simulated data dir not found: {SIM_DATA_DIR}"
            return
        if not (os.path.exists(EH_BINARY) and os.access(EH_BINARY, os.X_OK)):
            cls.SKIP_REASON = f"ExpansionHunter binary not found/executable: {EH_BINARY}"
            return
        reference = find_reference()
        if reference is None:
            cls.SKIP_REASON = "no hg38 reference FASTA (+.fai) found among candidates"
            return

        cls.tmp_dir = tempfile.mkdtemp(prefix="eh_mode_test_")
        os.makedirs(MERGED_BAM_CACHE_DIR, exist_ok=True)

        loci = sorted(
            filter(None, (parse_locus_dir(d) for d in os.listdir(SIM_DATA_DIR)
                          if os.path.isdir(os.path.join(SIM_DATA_DIR, d)))))

        cls.RESULTS = {}
        cls.truth = {}
        n_samples = 0
        for locus_id, chrom, start, end, motif in loci:
            locus_dir = os.path.join(SIM_DATA_DIR, locus_id)
            pairs = select_pairs(allele_sizes_in_dir(locus_dir))
            if not pairs:
                continue
            catalog = os.path.join(cls.tmp_dir, f"{locus_id}.catalog.json")
            write_catalog(locus_id, chrom, start, end, motif, catalog)
            for n1, n2 in pairs:
                merged_bam = os.path.join(MERGED_BAM_CACHE_DIR, locus_id, f"a1_{n1}__a2_{n2}.bam")
                build_merged_bam(
                    os.path.join(locus_dir, f"sim_{n1}x__10_150_450_50.bam"),
                    os.path.join(locus_dir, f"sim_{n2}x__10_150_450_50.bam"),
                    merged_bam)
                key = (locus_id, n1, n2)
                cls.truth[key] = tuple(sorted((n1, n2)))
                cls.RESULTS[key] = {}
                for mode in ALL_MODES:
                    cls.RESULTS[key][mode] = run_expansion_hunter(
                        merged_bam, reference, catalog,
                        os.path.join(cls.tmp_dir, f"{locus_id}.{n1}_{n2}.{mode}"), mode)
                n_samples += 1
        print(f"\nGenotyped {n_samples} simulated diploid samples "
              f"across {len({k[0] for k in cls.RESULTS})} loci x {len(ALL_MODES)} modes.")

    @classmethod
    def tearDownClass(cls):
        if getattr(cls, "tmp_dir", None) and os.path.isdir(cls.tmp_dir):
            shutil.rmtree(cls.tmp_dir, ignore_errors=True)

    def setUp(self):
        if self.SKIP_REASON:
            self.skipTest(self.SKIP_REASON)
        if not self.RESULTS:
            self.skipTest("no simulated samples were genotyped")

    def test_seeking_streaming_lowmem_produce_identical_calls(self):
        """seeking, streaming and low-mem-streaming must agree on GT/REPCN/REPCI per locus."""
        mismatches = 0
        for key in sorted(self.RESULTS):
            locus_id, n1, n2 = key
            with self.subTest(locus=locus_id, truth=f"{n1}/{n2}"):
                calls = {m: self.RESULTS[key][m] for m in IDENTICAL_MODES}
                reference_call = calls["seeking"]
                signature = (reference_call["GT"], reference_call["REPCN"], reference_call["REPCI"])
                for mode in IDENTICAL_MODES[1:]:
                    other = (calls[mode]["GT"], calls[mode]["REPCN"], calls[mode]["REPCI"])
                    if other != signature:
                        mismatches += 1
                    self.assertEqual(
                        other, signature,
                        f"{mode} disagrees with seeking at {locus_id} (truth {n1}/{n2}): "
                        f"seeking={signature} {mode}={other}")
        self.assertEqual(mismatches, 0)

    def test_seeking_vs_optimized_streaming_accuracy(self):
        """Reports and compares accuracy of seeking vs optimized-streaming against truth."""
        stats = {mode: {"n": 0, "exact": 0, "alleles_matched": 0,
                        "within1": 0, "no_call": 0} for mode in ACCURACY_MODES}
        seeking_correct_optimized_wrong = []

        for key in sorted(self.RESULTS):
            locus_id, n1, n2 = key
            truth = self.truth[key]
            per_mode_exact = {}
            for mode in ACCURACY_MODES:
                call = self.RESULTS[key][mode]["alleles"]
                s = stats[mode]
                s["n"] += 1
                if call is None:
                    s["no_call"] += 1
                exact = call == truth
                per_mode_exact[mode] = exact
                if exact:
                    s["exact"] += 1
                s["alleles_matched"] += allele_match_count(call, truth)
                if within_one(call, truth):
                    s["within1"] += 1
            if per_mode_exact["seeking"] and not per_mode_exact["optimized-streaming"]:
                seeking_correct_optimized_wrong.append(
                    (locus_id, f"{n1}/{n2}", self.RESULTS[key]["optimized-streaming"]["REPCN"]))

        print("\n=== Accuracy vs known truth (seeking vs optimized-streaming) ===")
        print(f"{'mode':<20} {'N':>4} {'exact_GT':>12} {'allele_acc':>12} "
              f"{'within±1':>12} {'no_call':>9}")
        for mode in ACCURACY_MODES:
            s = stats[mode]
            print(f"{mode:<20} {s['n']:>4} "
                  f"{s['exact']:>6} ({100*s['exact']/s['n']:4.0f}%) "
                  f"{100*s['alleles_matched']/(2*s['n']):>10.0f}% "
                  f"{s['within1']:>6} ({100*s['within1']/s['n']:4.0f}%) "
                  f"{s['no_call']:>9}")

        if seeking_correct_optimized_wrong:
            print(f"\n{len(seeking_correct_optimized_wrong)} samples where seeking is exact but "
                  f"optimized-streaming is not (locus, truth, optimized REPCN):")
            for locus_id, truth_str, optimized_repcn in seeking_correct_optimized_wrong[:40]:
                print(f"  {locus_id:<40} truth={truth_str:<8} optimized={optimized_repcn}")
            if len(seeking_correct_optimized_wrong) > 40:
                print(f"  ... and {len(seeking_correct_optimized_wrong) - 40} more")

        # Sanity guard: the exact seeking path must get a non-trivial fraction right,
        # otherwise the harness (catalog/reference/coverage) is broken rather than the modes.
        self.assertGreater(stats["seeking"]["exact"], 0,
                           "seeking produced no exact genotype matches — check harness inputs")


if __name__ == "__main__":
    unittest.main(verbosity=2)
