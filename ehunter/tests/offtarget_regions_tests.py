"""End-to-end tests for catalog loci that list OfftargetRegions.

In optimized-streaming mode such a locus is genotyped from reads collected the way seeking mode collects
them, so its record must match seeking mode's exactly, under every --genotyping-approach and --threads value.

The fixture under ``ehunter/tests/data/offtarget_regions_fixtures`` has a GGCCCC locus (C9LIKE) with an
expanded allele that is visible mostly through in-repeat read pairs placed in its off-target regions, and an
ordinary CAG locus (CHR2_CAG) without off-target regions. See its README.md.

Run from the repo root with:  python3 ehunter/tests/offtarget_regions_tests.py -v
"""

import json
import os
import shutil
import subprocess
import tempfile
import unittest

from resume_tests import find_binary


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(TESTS_DIR, "data", "offtarget_regions_fixtures")
OFFTARGET_LOCUS_ID = "C9LIKE"
ORDINARY_LOCUS_ID = "CHR2_CAG"
# Read pairs placed in C9LIKE's off-target regions (see the fixture's README.md).
OFFTARGET_READ_PAIR_COUNT = 67


def vcf_records_by_variant_id(path):
    """Returns {VARID: data line} for the records in a VCF.

    The ID column is "." for every record, so records are told apart by the VARID INFO field.
    """
    records = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            info = dict(field.split("=", 1) for field in line.split("\t")[7].split(";") if "=" in field)
            records[info["VARID"]] = line
    return records


class OfftargetRegionsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.binary = find_binary()
        print(f"\nTesting binary: {cls.binary}")
        if cls.binary is None:
            raise unittest.SkipTest("no ExpansionHunter binary found under ehunter/build*/")

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix="eh_offtarget_regions_test_")
        self.addCleanup(shutil.rmtree, self.work_dir, ignore_errors=True)

    def write_modified_catalog(self, name, remove_offtarget_regions, offtarget_locus_variant_type=None):
        """Returns the path of a modified copy of the fixture catalog.

        Args:
            name: file name for the copy, inside the test's work directory.
            remove_offtarget_regions: if True, every OfftargetRegions field is removed.
            offtarget_locus_variant_type: if given, the VariantType to set on C9LIKE.
        """
        with open(os.path.join(FIXTURE_DIR, "variant_catalog.json")) as f:
            catalog = json.load(f)
        for locus in catalog:
            if remove_offtarget_regions:
                locus.pop("OfftargetRegions", None)
            if offtarget_locus_variant_type and locus["LocusId"] == OFFTARGET_LOCUS_ID:
                locus["VariantType"] = offtarget_locus_variant_type
        path = os.path.join(self.work_dir, name)
        with open(path, "w") as f:
            json.dump(catalog, f, indent=4)
        return path

    def write_catalog_without_offtarget_regions(self):
        return self.write_modified_catalog("catalog_without_offtarget_regions.json", remove_offtarget_regions=True)

    def run_eh(self, name, extra_args, catalog=None, expect_success=True):
        """Runs ExpansionHunter on the fixture and returns the output prefix."""
        prefix = os.path.join(self.work_dir, name)
        command = [
            self.binary,
            "--reads", os.path.join(FIXTURE_DIR, "reads.bam"),
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--variant-catalog", catalog or os.path.join(FIXTURE_DIR, "variant_catalog.json"),
            "--output-prefix", prefix,
        ] + extra_args
        result = subprocess.run(command, capture_output=True, text=True)
        if expect_success and result.returncode != 0:
            self.fail(f"{' '.join(command)}\nfailed with exit code {result.returncode}:\n{result.stderr}")
        return prefix

    def locus_results(self, prefix):
        with open(prefix + ".json") as f:
            return json.load(f)["LocusResults"]

    def test_offtarget_read_pairs_change_the_seeking_mode_call(self):
        # Guards the fixture: the comparisons below only exercise off-target read pairs if dropping the
        # off-target regions changes the call.
        with_regions = self.locus_results(self.run_eh("seeking", ["--analysis-mode", "seeking"]))
        without_regions = self.locus_results(self.run_eh(
            "seeking_without_offtarget_regions", ["--analysis-mode", "seeking"],
            catalog=self.write_catalog_without_offtarget_regions()))
        self.assertNotEqual(
            with_regions[OFFTARGET_LOCUS_ID]["Variants"], without_regions[OFFTARGET_LOCUS_ID]["Variants"])

    def test_optimized_streaming_matches_seeking_for_locus_with_offtarget_regions(self):
        seeking_prefix = self.run_eh("seeking", ["--analysis-mode", "seeking"])
        expected_record = self.locus_results(seeking_prefix)[OFFTARGET_LOCUS_ID]
        expected_vcf_record = vcf_records_by_variant_id(seeking_prefix + ".vcf")[OFFTARGET_LOCUS_ID]

        for approach in ("auto", "only-full", "only-quick"):
            for thread_count in ("1", "2"):
                with self.subTest(genotyping_approach=approach, threads=thread_count):
                    prefix = self.run_eh(
                        f"optimized_{approach}_threads{thread_count}",
                        ["--analysis-mode", "optimized-streaming", "--genotyping-approach", approach,
                         "--threads", thread_count])
                    self.assertEqual(expected_record, self.locus_results(prefix)[OFFTARGET_LOCUS_ID])
                    self.assertEqual(
                        expected_vcf_record, vcf_records_by_variant_id(prefix + ".vcf")[OFFTARGET_LOCUS_ID])

    def test_optimized_streaming_matches_seeking_with_motif_composition(self):
        # Both modes compute the locus's MotifComposition with the same function from the same reads.
        args = ["--output-motif-composition", "all-loci"]
        expected_record = self.locus_results(
            self.run_eh("seeking", ["--analysis-mode", "seeking"] + args))[OFFTARGET_LOCUS_ID]
        self.assertIn("MotifComposition", list(expected_record["Variants"].values())[0])
        record = self.locus_results(
            self.run_eh("optimized", ["--analysis-mode", "optimized-streaming"] + args))[OFFTARGET_LOCUS_ID]
        self.assertEqual(expected_record, record)

    def test_offtarget_reads_are_counted_in_motif_composition(self):
        # Every read of the fixture's 67 off-target pairs is pure GGCCCC repeat sequence, stored on the strand of its
        # off-target region (half of them on the opposite strand to the locus). Each one must be turned to the
        # locus's strand and counted once.
        args = ["--analysis-mode", "seeking", "--output-motif-composition", "all-loci"]

        def reads_with_catalog_motif(locus_results):
            composition = locus_results[OFFTARGET_LOCUS_ID]["Variants"][OFFTARGET_LOCUS_ID]["MotifComposition"]
            return composition["Motifs"]["1:GGCCCC"]["reads"]

        with_regions = self.locus_results(self.run_eh("seeking", args))
        without_regions = self.locus_results(self.run_eh(
            "seeking_without_offtarget_regions", args, catalog=self.write_catalog_without_offtarget_regions()))
        self.assertEqual(
            reads_with_catalog_motif(with_regions) - reads_with_catalog_motif(without_regions),
            2 * OFFTARGET_READ_PAIR_COUNT)

    def test_resumed_run_matches_seeking_for_locus_with_offtarget_regions(self):
        expected_record = self.locus_results(self.run_eh("seeking", ["--analysis-mode", "seeking"]))[OFFTARGET_LOCUS_ID]
        resume_args = ["--analysis-mode", "optimized-streaming", "--resume"]
        # The catalog is sorted by position, so C9LIKE (chr1) is the first locus checkpointed.
        self.run_eh("resumed", resume_args + ["--internal-abort-after-loci", "1"], expect_success=False)
        prefix = self.run_eh("resumed", resume_args)
        self.assertEqual(expected_record, self.locus_results(prefix)[OFFTARGET_LOCUS_ID])

    def test_locus_without_offtarget_regions_is_unaffected(self):
        # CHR2_CAG shares chr2 with one of C9LIKE's off-target regions; its record must not depend on whether
        # C9LIKE lists them.
        args = ["--analysis-mode", "optimized-streaming"]
        with_regions = self.locus_results(self.run_eh("optimized", args))
        without_regions = self.locus_results(self.run_eh(
            "optimized_without_offtarget_regions", args, catalog=self.write_catalog_without_offtarget_regions()))
        self.assertEqual(with_regions[ORDINARY_LOCUS_ID], without_regions[ORDINARY_LOCUS_ID])

    def test_offtarget_regions_of_a_common_repeat_are_ignored(self):
        # Off-target read pairs count only toward a RareRepeat, so optimized-streaming mode genotypes a plain Repeat
        # from its streamed reads as if it listed no off-target regions.
        args = ["--analysis-mode", "optimized-streaming", "--genotyping-approach", "only-full"]
        with_regions = self.locus_results(self.run_eh("with_regions", args, catalog=self.write_modified_catalog(
            "common_repeat_with_regions.json", remove_offtarget_regions=False, offtarget_locus_variant_type="Repeat")))
        without_regions = self.locus_results(self.run_eh("without_regions", args, catalog=self.write_modified_catalog(
            "common_repeat_without_regions.json", remove_offtarget_regions=True, offtarget_locus_variant_type="Repeat")))
        self.assertEqual(with_regions[OFFTARGET_LOCUS_ID], without_regions[OFFTARGET_LOCUS_ID])


if __name__ == "__main__":
    unittest.main()
