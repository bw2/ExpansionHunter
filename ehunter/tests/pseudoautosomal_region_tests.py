"""End-to-end tests for the IsInsidePseudoautosomalRegion catalog field.

The fixture under ``ehunter/tests/data/pseudoautosomal_region_fixtures`` has 14 CAG loci, each heterozygous for a
12-unit and a 16-unit allele: PAR_CAG_01 ... PAR_CAG_12 on chrX, which its catalog marks
"IsInsidePseudoautosomalRegion": true, NONPAR_CAG on chrX, and AUTOSOMAL_CAG on chr1. See its README.md.

A marked locus is genotyped with two alleles in every sample, while other chrX loci keep one allele in a male sample.

Run from the repo root with:  python3 ehunter/tests/pseudoautosomal_region_tests.py -v
"""

import json
import os
import shutil
import subprocess
import tempfile
import unittest

from offtarget_regions_tests import vcf_records_by_variant_id
from resume_tests import find_binary


TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
FIXTURE_DIR = os.path.join(TESTS_DIR, "data", "pseudoautosomal_region_fixtures")
CATALOG = os.path.join(FIXTURE_DIR, "variant_catalog.json")
BAM = os.path.join(FIXTURE_DIR, "reads.bam")

PSEUDOAUTOSOMAL_LOCUS_IDS = ["PAR_CAG_%02d" % (i + 1) for i in range(12)]
UNMARKED_LOCUS_IDS = ["NONPAR_CAG", "AUTOSOMAL_CAG"]
CHRX_LOCUS_IDS = PSEUDOAUTOSOMAL_LOCUS_IDS + ["NONPAR_CAG"]

# (name, analysis-mode arguments): every way of genotyping a locus.
MODES = [
    ("seeking", ["--analysis-mode", "seeking"]),
    ("streaming", ["--analysis-mode", "streaming"]),
    ("optimized_auto", ["--analysis-mode", "optimized-streaming", "--genotyping-approach", "auto"]),
    ("optimized_only_full", ["--analysis-mode", "optimized-streaming", "--genotyping-approach", "only-full"]),
    ("optimized_only_quick", ["--analysis-mode", "optimized-streaming", "--genotyping-approach", "only-quick"]),
]


class PseudoautosomalRegionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.binary = find_binary()
        print(f"\nTesting binary: {cls.binary}")
        if cls.binary is None:
            raise unittest.SkipTest("no ExpansionHunter binary found under ehunter/build*/")

    def setUp(self):
        self.work_dir = tempfile.mkdtemp(prefix="eh_pseudoautosomal_region_test_")
        self.addCleanup(shutil.rmtree, self.work_dir, ignore_errors=True)

    def run_eh(self, name, extra_args, catalog=CATALOG, expect_success=True):
        """Runs ExpansionHunter on the fixture and returns (output prefix, completed process)."""
        prefix = os.path.join(self.work_dir, name)
        command = [
            self.binary,
            "--reads", BAM,
            "--reference", os.path.join(FIXTURE_DIR, "reference.fa"),
            "--variant-catalog", catalog,
            "--output-prefix", prefix,
        ] + extra_args
        result = subprocess.run(command, capture_output=True, text=True)
        if expect_success and result.returncode != 0:
            self.fail(f"{' '.join(command)}\nfailed with exit code {result.returncode}:\n{result.stdout}{result.stderr}")
        return prefix, result

    def write_catalog(self, name, change_record):
        """Writes a copy of the fixture catalog with change_record(record) applied to each record."""
        with open(CATALOG) as f:
            records = json.load(f)
        for record in records:
            change_record(record)
        path = os.path.join(self.work_dir, name)
        with open(path, "w") as f:
            json.dump(records, f)
        return path

    def catalog_without_the_field(self):
        return self.write_catalog(
            "catalog_without_the_field.json", lambda record: record.pop("IsInsidePseudoautosomalRegion", None))

    def locus_results(self, prefix):
        with open(prefix + ".json") as f:
            return json.load(f)["LocusResults"]

    @staticmethod
    def genotype(locus_record):
        return list(locus_record["Variants"].values())[0]["Genotype"]

    def test_fixture_is_heterozygous_at_every_locus_in_a_female_sample(self):
        # Guards the fixture: the male-sample tests below only show the allele count if each locus has two alleles.
        results = self.locus_results(self.run_eh("female", ["--sex", "female"], catalog=self.catalog_without_the_field())[0])
        for locus_id in CHRX_LOCUS_IDS + ["AUTOSOMAL_CAG"]:
            self.assertEqual("12/16", self.genotype(results[locus_id]), locus_id)

    def test_male_sample_without_the_field_has_one_allele_at_every_chrx_locus(self):
        catalog = self.catalog_without_the_field()
        for mode_name, mode_args in MODES:
            with self.subTest(mode=mode_name):
                results = self.locus_results(self.run_eh(mode_name, ["--sex", "male"] + mode_args, catalog=catalog)[0])
                for locus_id in CHRX_LOCUS_IDS:
                    self.assertEqual(1, results[locus_id]["AlleleCount"], locus_id)
                    self.assertNotIn("/", self.genotype(results[locus_id]), locus_id)
                self.assertEqual(2, results["AUTOSOMAL_CAG"]["AlleleCount"])

    def test_male_sample_genotypes_only_the_marked_loci_with_two_alleles(self):
        catalog_without_the_field = self.catalog_without_the_field()
        for mode_name, mode_args in MODES:
            with self.subTest(mode=mode_name):
                default_prefix, _ = self.run_eh(
                    mode_name + "_without_the_field", ["--sex", "male"] + mode_args, catalog=catalog_without_the_field)
                marked_prefix, _ = self.run_eh(mode_name, ["--sex", "male"] + mode_args)
                default_results = self.locus_results(default_prefix)
                marked_results = self.locus_results(marked_prefix)
                marked_vcf = vcf_records_by_variant_id(marked_prefix + ".vcf")

                for locus_id in PSEUDOAUTOSOMAL_LOCUS_IDS:
                    self.assertEqual(2, marked_results[locus_id]["AlleleCount"], locus_id)
                    self.assertEqual("12/16", self.genotype(marked_results[locus_id]), locus_id)
                    self.assertEqual("12/16", marked_vcf[locus_id].rstrip("\n").split("\t")[9].split(":")[2], locus_id)

                default_vcf = vcf_records_by_variant_id(default_prefix + ".vcf")
                for locus_id in UNMARKED_LOCUS_IDS:
                    self.assertEqual(default_results[locus_id], marked_results[locus_id], locus_id)
                    self.assertEqual(default_vcf[locus_id], marked_vcf[locus_id], locus_id)

    def test_field_changes_nothing_for_a_female_sample(self):
        default_prefix, _ = self.run_eh("female_without_the_field", ["--sex", "female"],
                                        catalog=self.catalog_without_the_field())
        marked_prefix, _ = self.run_eh("female", ["--sex", "female"])
        self.assertEqual(self.locus_results(default_prefix), self.locus_results(marked_prefix))
        self.assertEqual(vcf_records_by_variant_id(default_prefix + ".vcf"),
                         vcf_records_by_variant_id(marked_prefix + ".vcf"))

    def test_field_set_to_true_on_an_autosomal_locus_is_rejected(self):
        def mark_autosomal_locus(record):
            if record["LocusId"] == "AUTOSOMAL_CAG":
                record["IsInsidePseudoautosomalRegion"] = True
        prefix, result = self.run_eh("male_marked_autosome", ["--sex", "male"],
                                     catalog=self.write_catalog("marked_autosome.json", mark_autosomal_locus),
                                     expect_success=False)
        self.assertNotEqual(0, result.returncode)
        self.assertIn('Locus AUTOSOMAL_CAG has "IsInsidePseudoautosomalRegion": true, but is not on chrX',
                      result.stdout + result.stderr)
        self.assertFalse(os.path.exists(prefix + ".json"))

    def test_value_other_than_true_or_false_is_rejected(self):
        def set_string_value(record):
            if "IsInsidePseudoautosomalRegion" in record:
                record["IsInsidePseudoautosomalRegion"] = "true"
        _, result = self.run_eh("male_string_value", ["--sex", "male"],
                                catalog=self.write_catalog("string_value.json", set_string_value),
                                expect_success=False)
        self.assertNotEqual(0, result.returncode)
        self.assertIn('IsInsidePseudoautosomalRegion of locus PAR_CAG_01 must be true or false, not "true"',
                      result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
