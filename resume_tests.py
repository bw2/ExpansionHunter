"""End-to-end tests for ExpansionHunter's --resume flag.

Each test interrupts a run at a chosen point (via the internal --internal-abort-after-loci
hook, which kills the process without unwinding, exactly as an external kill would), resumes
it, and checks that the final output is identical to an uninterrupted run's.

Everything the tests need is in the repo: the fixture under
``ehunter/tests/data/parallel_processing_fixtures`` is a 5-contig synthetic reference with 4
loci, so no external reference genome or read data is required.

Run with:  python3 -m unittest resume_tests -v
"""

import glob
import gzip
import os
import re
import shutil
import subprocess
import tempfile
import unittest


REPO_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(REPO_DIR, "ehunter", "tests", "data", "parallel_processing_fixtures")


def is_sanitizer_build(binary_path):
    """True if the binary's build directory was configured with a sanitizer.

    This repo keeps ASan/TSan build directories alongside the normal ones. Those binaries run
    roughly 10x slower and are not what these tests mean to measure, so they are skipped unless
    named explicitly through $EH_BINARY.
    """
    cache = os.path.join(os.path.dirname(binary_path), "CMakeCache.txt")
    try:
        with open(cache) as f:
            return "-fsanitize" in f.read()
    except OSError:
        return False


def binary_candidates():
    """Returns the ExpansionHunter binaries to consider, best first.

    $EH_BINARY wins when it is set. Otherwise every ehunter/build*/ExpansionHunter is considered,
    since this repo keeps several build directories side by side, newest first so a fresh build
    beats a stale one.
    """
    if os.environ.get("EH_BINARY"):
        return [os.environ["EH_BINARY"]]
    paths = glob.glob(os.path.join(REPO_DIR, "ehunter", "build*", "ExpansionHunter"))
    return sorted(paths, key=os.path.getmtime, reverse=True)


def find_binary():
    """Returns the newest non-sanitizer binary that supports --resume, or None.

    The --resume check matters because a stale binary from an older build directory would
    otherwise be picked up and fail every test with an unhelpful "unrecognised option" error.
    """
    explicit = bool(os.environ.get("EH_BINARY"))
    for path in binary_candidates():
        if not (os.path.exists(path) and os.access(path, os.X_OK)):
            continue
        if not explicit and is_sanitizer_build(path):
            continue
        help_text = subprocess.run([path, "--help"], capture_output=True, text=True)
        if "--resume" in help_text.stdout + help_text.stderr:
            return path
    return None


def read_maybe_gzipped(path):
    """Returns the text content of path, decompressing it when it ends in .gz."""
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            return f.read()
    with open(path) as f:
        return f.read()


def output_body(path):
    """Returns an output file's content without the trailing RunInfo record.

    RunInfo holds timestamps, a runtime and a peak-memory figure, none of which can match
    across two separate runs, so comparisons are made against everything before it.
    """
    return read_maybe_gzipped(path).split('"RunInfo"')[0]


class ResumeTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.binary = find_binary()
        # Printed because several build directories coexist here and silently testing the wrong
        # binary would make every result meaningless.
        print(f"\nTesting binary: {cls.binary}")
        if cls.binary is None:
            raise unittest.SkipTest(
                "no ExpansionHunter binary supporting --resume found; looked at "
                + (", ".join(binary_candidates()) or "(nothing)"))

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix="eh_resume_test_")
        self.addCleanup(shutil.rmtree, self.work_dir, ignore_errors=True)

    def run_eh(self, prefix, extra_args, expect_success=True):
        """Runs ExpansionHunter on the fixture and returns the CompletedProcess."""
        command = [
            self.binary,
            "--reads", os.path.join(FIXTURE_DIR, "reads.bam"),
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--catalog", os.path.join(FIXTURE_DIR, "variant_catalog.json"),
            "--output-prefix", prefix,
        ] + extra_args
        if "--analysis-mode" not in extra_args:
            command += ["--analysis-mode", "optimized-streaming"]
        result = subprocess.run(command, capture_output=True, text=True)
        if expect_success and result.returncode != 0:
            self.fail(f"{' '.join(command)}\nfailed with exit code {result.returncode}:\n{result.stderr}")
        return result

    def baseline(self, extra_args=()):
        """Runs to completion without --resume and returns (json_body, vcf_text)."""
        prefix = os.path.join(self.work_dir, "baseline")
        self.run_eh(prefix, list(extra_args))
        suffix = ".gz" if "-z" in extra_args else ""
        return output_body(prefix + ".json" + suffix), read_maybe_gzipped(prefix + ".vcf" + suffix)

    def interrupt_and_resume(self, prefix, abort_after, extra_args=()):
        """Runs until abort_after loci are checkpointed, then resumes to completion."""
        self.run_eh(prefix, list(extra_args) + ["--resume", "--internal-abort-after-loci", str(abort_after)],
                    expect_success=False)
        return self.run_eh(prefix, list(extra_args) + ["--resume"])

    def assert_matches_baseline(self, prefix, extra_args=()):
        expected_json, expected_vcf = self.baseline(extra_args)
        suffix = ".gz" if "-z" in extra_args else ""
        self.assertEqual(expected_json, output_body(prefix + ".json" + suffix))
        self.assertEqual(expected_vcf, read_maybe_gzipped(prefix + ".vcf" + suffix))

    def assert_no_leftover_files(self, prefix):
        leftovers = [os.path.basename(p) for p in glob.glob(prefix + "*")
                     if ".unfinished" in p or ".contig" in p or ".resume_part" in p]
        self.assertEqual(leftovers, [], f"checkpoint/temp files left behind: {leftovers}")

    def test_resumed_output_matches_uninterrupted_output(self):
        # The fixture has 4 loci, so this covers interrupting before, during and after each of them,
        # at one thread and at several.
        for threads in (1, 2, 3):
            for abort_after in (1, 2, 3, 4):
                with self.subTest(threads=threads, abort_after=abort_after):
                    prefix = os.path.join(self.work_dir, f"t{threads}_a{abort_after}")
                    self.interrupt_and_resume(prefix, abort_after, ["--threads", str(threads)])
                    self.assert_matches_baseline(prefix)
                    self.assert_no_leftover_files(prefix)

    def test_resume_with_compressed_output(self):
        prefix = os.path.join(self.work_dir, "compressed")
        self.interrupt_and_resume(prefix, 2, ["--threads", "2", "-z"])
        self.assert_matches_baseline(prefix, ["--threads", "2", "-z"])
        self.assert_no_leftover_files(prefix)

    def test_resume_in_low_mem_streaming_mode(self):
        prefix = os.path.join(self.work_dir, "lowmem")
        self.interrupt_and_resume(prefix, 2, ["--threads", "2", "--analysis-mode", "low-mem-streaming"])
        self.assert_matches_baseline(prefix, ["--threads", "2", "--analysis-mode", "low-mem-streaming"])

    def test_resume_keeps_loci_that_produced_no_record(self):
        # --skip-hom-ref genotypes a locus and then emits nothing for it. Those loci have to be
        # remembered as finished anyway, or every resume would genotype them again.
        args = ["--threads", "2", "--skip-hom-ref"]
        prefix = os.path.join(self.work_dir, "skip_hom_ref")
        self.run_eh(prefix, args + ["--resume", "--internal-abort-after-loci", "2"], expect_success=False)
        with open(prefix + ".processed_loci.unfinished") as f:
            finished = [line.split("\t")[0] for line in f.read().splitlines()]
        # The abort hook fires once a whole batch has been checkpointed, so it can overshoot its N.
        self.assertGreaterEqual(len(finished), 2)
        self.assertLess(len(finished), 4)

        result = self.run_eh(prefix, args + ["--resume"])
        # Every finished locus here is one --skip-hom-ref emitted no record for, so seeing them counted as
        # already genotyped is the whole point: without the processed-loci list they would all be redone.
        self.assertIn(f"{len(finished)} of 4 loci were already genotyped", result.stdout + result.stderr)
        self.assert_matches_baseline(prefix, args)

    def test_resume_across_a_change_in_thread_count(self):
        prefix = os.path.join(self.work_dir, "mixed_threads")
        self.run_eh(prefix, ["--threads", "3", "--resume", "--internal-abort-after-loci", "2"],
                    expect_success=False)
        self.run_eh(prefix, ["--threads", "1", "--resume"])
        self.assert_matches_baseline(prefix)

    def test_resume_without_a_checkpoint_runs_normally(self):
        prefix = os.path.join(self.work_dir, "fresh")
        self.run_eh(prefix, ["--threads", "2", "--resume"])
        self.assert_matches_baseline(prefix)
        self.assert_no_leftover_files(prefix)

    def test_a_truncated_checkpoint_is_recovered(self):
        # Cut bytes off the end of the JSON checkpoint, which is what losing buffered writes to a
        # dying machine looks like. Whatever survives has to be used, and the rest genotyped again.
        for keep_fraction in (0.4, 0.7, 0.95):
            with self.subTest(keep_fraction=keep_fraction):
                prefix = os.path.join(self.work_dir, f"truncated_{keep_fraction}")
                self.run_eh(prefix, ["--threads", "2", "--resume", "--internal-abort-after-loci", "3"],
                            expect_success=False)
                checkpoint = prefix + ".json.unfinished"
                with open(checkpoint, "rb") as f:
                    content = f.read()
                with open(checkpoint, "wb") as f:
                    f.write(content[:int(len(content) * keep_fraction)])

                self.run_eh(prefix, ["--threads", "2", "--resume"])
                self.assert_matches_baseline(prefix)

    def test_a_multi_variant_locus_missing_one_vcf_line_is_regenotyped(self):
        # CHR4_MULTI in the fixture has two variants and so writes two VCF lines. Losing only the
        # second one must not leave the locus counted as finished: checking mere presence of the
        # locus in the VCF checkpoint would silently drop that variant from the final output.
        prefix = os.path.join(self.work_dir, "multi_variant")
        self.run_eh(prefix, ["--threads", "1", "--resume", "--internal-abort-after-loci", "4"],
                    expect_success=False)

        checkpoint = prefix + ".vcf.unfinished"
        with open(checkpoint) as f:
            lines = f.read().splitlines(keepends=True)
        self.assertTrue(lines[-1].startswith("chr4"), "expected the fixture's last VCF line to be CHR4_MULTI's")
        with open(checkpoint, "w") as f:
            f.write("".join(lines[:-1]))

        self.run_eh(prefix, ["--threads", "1", "--resume"])
        self.assert_matches_baseline(prefix)

    def test_a_truncated_checkpoint_still_preserves_earlier_progress(self):
        # Asserting only that the final output is correct would not catch a resume that silently threw
        # the whole checkpoint away and re-genotyped everything, which is what a truncated gzip tail used
        # to do. So check how many loci the resumed run actually skipped.
        for compress in (False, True):
            with self.subTest(compressed=compress):
                args = ["--threads", "1"] + (["-z"] if compress else [])
                prefix = os.path.join(self.work_dir, f"progress_{int(compress)}")
                self.run_eh(prefix, args + ["--resume", "--internal-abort-after-loci", "3"],
                            expect_success=False)

                checkpoint = prefix + (".json.gz.unfinished" if compress else ".json.unfinished")
                with open(checkpoint, "rb") as f:
                    content = f.read()
                with open(checkpoint, "wb") as f:
                    f.write(content[:-1])   # lose a single trailing byte

                result = self.run_eh(prefix, args + ["--resume"])
                output = result.stdout + result.stderr
                match = re.search(r"Resuming: (\d+) of \d+ loci were already genotyped", output)
                self.assertIsNotNone(match, f"resume recovered nothing from the checkpoint:\n{output}")
                self.assertGreater(int(match.group(1)), 0,
                                   "a one-byte truncation should not discard every finished locus")
                self.assert_matches_baseline(prefix, args)

    def test_a_checkpoint_from_a_different_run_is_rejected(self):
        prefix = os.path.join(self.work_dir, "mismatch")
        self.run_eh(prefix, ["--threads", "2", "--resume", "--internal-abort-after-loci", "2"],
                    expect_success=False)

        other_catalog = os.path.join(self.work_dir, "other_catalog.json")
        shutil.copy(os.path.join(FIXTURE_DIR, "variant_catalog.json"), other_catalog)
        result = subprocess.run([
            self.binary,
            "--reads", os.path.join(FIXTURE_DIR, "reads.bam"),
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--catalog", other_catalog,
            "--output-prefix", prefix,
            "--analysis-mode", "optimized-streaming",
            "--threads", "2", "--resume",
        ], capture_output=True, text=True)

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("written by a different run", result.stdout + result.stderr)
        # The checkpoint is left alone rather than silently discarded, so the run can still be resumed.
        self.assertTrue(os.path.exists(prefix + ".json.unfinished"))

    def test_an_unreadable_checkpoint_starts_over(self):
        prefix = os.path.join(self.work_dir, "corrupt")
        self.run_eh(prefix, ["--threads", "2", "--resume", "--internal-abort-after-loci", "2"],
                    expect_success=False)
        with open(prefix + ".json.unfinished", "w") as f:
            f.write("this is not an ExpansionHunter output file")

        result = self.run_eh(prefix, ["--threads", "2", "--resume"])
        self.assertIn("starting from the beginning", result.stdout + result.stderr)
        self.assert_matches_baseline(prefix)

    def test_resume_has_no_effect_in_seeking_mode(self):
        prefix = os.path.join(self.work_dir, "seeking")
        result = subprocess.run([
            self.binary,
            "--reads", os.path.join(FIXTURE_DIR, "reads.bam"),
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--catalog", os.path.join(FIXTURE_DIR, "variant_catalog.json"),
            "--output-prefix", prefix,
            "--analysis-mode", "seeking", "--resume",
        ], capture_output=True, text=True)

        self.assertEqual(result.returncode, 0)
        self.assertIn("--resume has no effect in seeking mode", result.stdout + result.stderr)
        self.assert_no_leftover_files(prefix)


if __name__ == "__main__":
    unittest.main()
