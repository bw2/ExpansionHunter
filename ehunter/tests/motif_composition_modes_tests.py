"""End-to-end tests for --output-motif-composition in each analysis mode.

Uses the fixture under ``ehunter/tests/data/parallel_processing_fixtures`` (see its README.md). Every read pair
there lies inside its locus's window, so seeking and optimized-streaming modes count motifs from the same reads
and must write the same MotifComposition records.

Run from the repo root with:  python3 ehunter/tests/motif_composition_modes_tests.py -v
"""

import json
import os
import shutil
import subprocess
import tempfile
import unittest

from resume_tests import find_binary


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(TESTS_DIR, "data", "parallel_processing_fixtures")


def motif_compositions(prefix):
    """Returns {VariantId: MotifComposition record or None} from a run's JSON output."""
    with open(prefix + ".json") as f:
        locus_results = json.load(f)["LocusResults"]
    return {
        variant_id: variant.get("MotifComposition")
        for locus in locus_results.values()
        for variant_id, variant in locus["Variants"].items()
    }


class MotifCompositionModesTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.binary = find_binary()
        print(f"\nTesting binary: {cls.binary}")
        if cls.binary is None:
            raise unittest.SkipTest("no ExpansionHunter binary found under ehunter/build*/")

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix="eh_motif_composition_modes_test_")
        self.addCleanup(shutil.rmtree, self.work_dir, ignore_errors=True)

    def run_eh(self, name, extra_args):
        """Runs ExpansionHunter on the fixture and returns (output prefix, CompletedProcess)."""
        prefix = os.path.join(self.work_dir, name)
        command = [
            self.binary,
            "--reads", os.path.join(FIXTURE_DIR, "reads.bam"),
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--variant-catalog", os.path.join(FIXTURE_DIR, "variant_catalog.json"),
            "--output-prefix", prefix,
            "--output-motif-composition", "all-loci",
        ] + extra_args
        return prefix, subprocess.run(command, capture_output=True, text=True)

    def test_seeking_mode_matches_optimized_streaming(self):
        seeking_prefix, result = self.run_eh("seeking", ["--analysis-mode", "seeking", "--threads", "2"])
        self.assertEqual(result.returncode, 0, result.stderr)
        optimized_prefix, result = self.run_eh("optimized", ["--analysis-mode", "optimized-streaming"])
        self.assertEqual(result.returncode, 0, result.stderr)

        seeking = motif_compositions(seeking_prefix)
        # CHR1_CAG, CHR2_CCG and the repeat of CHR4_MULTI have reads; CHR3_GAA has none and gets an empty record.
        self.assertTrue(seeking["CHR1_CAG"] and seeking["CHR2_CCG"] and seeking["CHR4_MULTI_chr4:1405-1475"])
        self.assertEqual(seeking["CHR3_GAA"], {})
        self.assertEqual(seeking, motif_compositions(optimized_prefix))

    def test_streaming_mode_is_rejected(self):
        prefix, result = self.run_eh("streaming", ["--analysis-mode", "streaming"])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("--output-motif-composition only works in seeking and optimized-streaming modes",
                      result.stdout + result.stderr)
        self.assertFalse(os.path.exists(prefix + ".json"))


if __name__ == "__main__":
    unittest.main()
