import contextlib
import gzip
import io
import json
import os
import tempfile
import unittest
from unittest import mock

import expansion_hunter_json_streaming_utils as streaming_utils
import generate_multisample_qc_stats as m


def repeat_record(variant_id, region, motif, genotype=None, counts_of_spanning_reads=None):
    """Returns an ExpansionHunter repeat record."""
    record = {"VariantId": variant_id, "VariantType": "Repeat", "ReferenceRegion": region, "RepeatUnit": motif}
    if genotype is not None:
        record["Genotype"] = genotype
    if counts_of_spanning_reads is not None:
        record["CountsOfSpanningReads"] = counts_of_spanning_reads
    return record


def expansion_hunter_json(sample_id, sex, records, read_length=None):
    """Returns an optimized-streaming ExpansionHunter output JSON with one single-repeat locus per record, laid out the
    way ExpansionHunter writes it: SampleParameters, RunInfo, LocusResults, RunStats. Each locus gets the given
    ReadLength, if any."""
    locus_results = {}
    for record in records:
        locus_record = {"LocusId": record["VariantId"], "Variants": {record["VariantId"]: record}}
        if read_length is not None:
            locus_record["ReadLength"] = read_length
        locus_results[record["VariantId"]] = locus_record
    return {
        "SampleParameters": {"SampleId": sample_id, "Sex": sex},
        "RunInfo": {"AnalysisMode": "optimized-streaming"},
        "LocusResults": locus_results,
        "RunStats": {},
    }


# L1: D has the only allele of size 20, larger than the most common size (10).
# L2: C has the only allele of size 7, smaller than the most common size (8). D has no call there.
# L3: A and B share the largest allele (9).
SAMPLE_JSONS = [
    expansion_hunter_json("A", "Female", [
        repeat_record("L1", "chr1:100-130", "CAG", "10/10"),
        repeat_record("L2", "chr1:500-530", "GCC", "8/8"),
        repeat_record("L3", "chr2:100-120", "AAT", "9/9"),
    ]),
    expansion_hunter_json("B", "Male", [
        repeat_record("L1", "chr1:100-130", "CAG", "10/12"),
        repeat_record("L2", "chr1:500-530", "GCC", "8/8"),
        repeat_record("L3", "chr2:100-120", "AAT", "5/9"),
    ]),
    expansion_hunter_json("C", "Female", [
        repeat_record("L1", "chr1:100-130", "CAG", "10/10"),
        repeat_record("L2", "chr1:500-530", "GCC", "7/8"),
        repeat_record("L3", "chr2:100-120", "AAT", "5/5"),
    ]),
    expansion_hunter_json("D", "Male", [
        repeat_record("L1", "chr1:100-130", "CAG", "10/20"),
        repeat_record("L2", "chr1:500-530", "GCC"),
        repeat_record("L3", "chr2:100-120", "AAT", "5/5"),
    ]),
]


def read_text(path):
    with gzip.open(path, "rt") as f:
        return f.read()


def read_qc_stats_file(path):
    """Returns (SampleIds of the header, {VariantId: locus row}, {SampleId: sample row}) of a QC stats file."""
    sample_ids = None
    locus_rows = {}
    sample_rows = {}
    for line in read_text(path).splitlines():
        fields = line.split("\t")
        if fields[0] == m.SAMPLE_IDS_LINE_START:
            sample_ids = fields[1:]
        elif fields[0] == "L":
            row = dict(zip(m.LOCUS_COLUMNS, fields[1:]))
            locus_rows[row["VariantId"]] = row
        elif fields[0] == "S":
            row = dict(zip(m.SAMPLE_COLUMNS, fields[1:]))
            sample_rows[row["SampleId"]] = row
    return sample_ids, locus_rows, sample_rows


def write_json(directory, name, data):
    """Writes an ExpansionHunter JSON to <directory>/<name>.json.gz, and returns its path."""
    path = os.path.join(directory, name + ".json.gz")
    with gzip.open(path, "wt") as f:
        json.dump(data, f)
    return path


def write_sample_jsons(directory, sample_jsons):
    """Writes each ExpansionHunter JSON to <directory>/<its SampleId>.json.gz, and returns their paths."""
    return [write_json(directory, data["SampleParameters"]["SampleId"], data) for data in sample_jsons]


def run_generate_multisample_qc_stats(directory, output_name, *args):
    """Runs generate_multisample_qc_stats.py's main with -o <directory>/<output_name>, and returns the output path."""
    output_prefix = os.path.join(directory, output_name)
    with contextlib.redirect_stdout(io.StringIO()):
        m.main(["-o", output_prefix] + [str(arg) for arg in args])
    return output_prefix + m.QC_STATS_SUFFIX


class GenerateMultisampleQcStatsTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir_object = tempfile.TemporaryDirectory()
        self.temp_dir = self.temp_dir_object.name

    def tearDown(self):
        self.temp_dir_object.cleanup()

    def write_json(self, name, data):
        return write_json(self.temp_dir, name, data)

    def write_sample_jsons(self, sample_jsons):
        return write_sample_jsons(self.temp_dir, sample_jsons)

    def run_main(self, output_name, *args):
        return run_generate_multisample_qc_stats(self.temp_dir, output_name, *args)

    def test_file_layout(self):
        lines = read_text(self.run_main("out", *self.write_sample_jsons(SAMPLE_JSONS))).splitlines()
        self.assertEqual(lines[:3], [m.FILE_TYPE_LINE, "##SampleIds\tA\tB\tC\tD", m.LOCI_HEADER_LINE])
        self.assertEqual([line.split("\t")[0] for line in lines[3:]], ["L", "L", "L", "#Samples", "S", "S", "S", "S"])

    def test_locus_rows(self):
        _, locus_rows, _ = read_qc_stats_file(self.run_main("out", *self.write_sample_jsons(SAMPLE_JSONS)))
        self.assertEqual(list(locus_rows), ["L1", "L2", "L3"])
        l1 = locus_rows["L1"]
        self.assertEqual((l1["LocusId"], l1["ReferenceRegion"], l1["Motif"], l1["SampleCount"],
                          l1["GenotypedSampleCount"]), ("L1", "chr1:100-130", "CAG", "4", "4"))
        self.assertEqual(l1["AlleleSizeHistogram"], "10x:6,12x:1,20x:1")
        self.assertEqual((l1["LargestAlleleSize"], l1["NumSamplesWithLargestAllele"],
                          l1["OnlySampleWithLargestAllele"]), ("20", "1", "D"))
        # the smallest allele (10) is shared by all four samples
        self.assertEqual((l1["SmallestAlleleSize"], l1["NumSamplesWithSmallestAllele"],
                          l1["OnlySampleWithSmallestAllele"]), ("10", "4", ""))
        # B (10/12) and D (10/20) of the 4 two-allele calls are heterozygous
        self.assertEqual((l1["DiploidGenotypedSampleCount"], l1["HeterozygousSampleCount"]), ("4", "2"))
        l2 = locus_rows["L2"]
        self.assertEqual((l2["SampleCount"], l2["GenotypedSampleCount"], l2["OnlySampleWithSmallestAllele"]),
                         ("4", "3", "C"))
        self.assertEqual((locus_rows["L3"]["LargestAlleleSize"], locus_rows["L3"]["NumSamplesWithLargestAllele"]),
                         ("9", "2"))

    def test_sample_rows(self):
        _, _, sample_rows = read_qc_stats_file(self.run_main("out", *self.write_sample_jsons(SAMPLE_JSONS)))
        self.assertEqual(list(sample_rows), ["A", "B", "C", "D"])
        d = sample_rows["D"]
        # D has no call at L2
        self.assertEqual((d["Sex"], d["RepeatRecords"], d["GenotypedRepeats"], d["CallRate"]),
                         ("Male", "3", "2", "0.66666667"))
        self.assertTrue(d["Path"].endswith("D.json.gz"))
        # SAMPLE_JSONS has no Coverage, ReadLength or FragmentLength
        self.assertEqual((d["MedianLocusCoverage"], d["MostCommonReadLength"], d["MedianFragmentLength"]), ("", "", ""))
        self.assertEqual({sample_id: (row["AutosomalDiploidGenotypes"], row["HeterozygosityRate"])
                          for sample_id, row in sample_rows.items()},
                         {"A": ("3", "0"), "B": ("3", "0.66666667"), "C": ("3", "0.33333333"), "D": ("2", "0.5")})
        # The reference sizes of L1 and L2 are 30 / 3 = 10 repeat units, and of L3 20 / 3 = 6 (rounded down, as in
        # ExpansionHunter's VCF REF). A's alleles differ from them by 0, 0, -2, -2, 3 and 3, and D's by 0, 10, -1, -1.
        self.assertEqual((sample_rows["A"]["MeanDifferenceFromReference"], d["MeanDifferenceFromReference"]),
                         ("0.33333333", "2"))

    def test_coverage_read_length_and_heterozygosity_rate(self):
        sample = expansion_hunter_json("A", "Male", [
            repeat_record("L1", "chr1:100-130", "CAG", "10/12"),
            repeat_record("L2", "chr1:500-530", "CAG", "10/10"),
            repeat_record("L3", "chr2:100-130", "CAG", "9/10"),
            # chrX calls don't count toward HeterozygosityRate, even with two alleles (as in a pseudoautosomal region)
            repeat_record("X1", "chrX:100-130", "CAG", "7/8"),
        ], read_length=151)
        for locus_id, coverage in [("L1", 30.4), ("L2", 34.6), ("L3", 12.0), ("X1", 20.0)]:
            sample["LocusResults"][locus_id]["Coverage"] = coverage
        # ExpansionHunter writes a ReadLength of 0 at loci without reads, which isn't a read length
        for locus_id in ("L2", "L3"):
            sample["LocusResults"][locus_id]["ReadLength"] = 0
        _, _, sample_rows = read_qc_stats_file(self.run_main("out", self.write_json("A", sample)))
        row = sample_rows["A"]
        # coverages rounded to 30, 35, 12 and 20: the median is (20 + 30) / 2
        self.assertEqual((row["MostCommonReadLength"], row["MedianLocusCoverage"]), ("151", "25"))
        self.assertEqual((row["AutosomalDiploidGenotypes"], row["HeterozygosityRate"]), ("3", "0.66666667"))

    def test_sex_check_coverage_ratios(self):
        def sample_with_coverage(sample_id, sex, chr_x_coverage, chr_y_coverage):
            # autosomal repeats at 30x, and a chrX repeat in a pseudoautosomal region, also at 30x, which the sex
            # check leaves out (an XY sample's chrX median would otherwise be (15 + 30) / 2, or 0.75 of the autosomes')
            coverage_by_locus = {"A1": 30, "A2": 30, "A3": 30, "XPAR": 30, "X1": chr_x_coverage, "Y1": chr_y_coverage}
            sample = expansion_hunter_json(sample_id, sex, [
                repeat_record("A1", "chr1:1000-1030", "CAG", "10/10"),
                repeat_record("A2", "chr2:1000-1030", "CAG", "10/10"),
                repeat_record("A3", "chr3:1000-1030", "CAG", "10/10"),
                repeat_record("XPAR", "chrX:100000-100030", "CAG", "10/10"),
                repeat_record("X1", "chrX:5000000-5000030", "CAG", "10/10"),
                repeat_record("Y1", "chrY:5000000-5000030", "CAG"),
            ])
            for locus_id, coverage in coverage_by_locus.items():
                sample["LocusResults"][locus_id]["Coverage"] = coverage
            return sample

        _, _, sample_rows = read_qc_stats_file(self.run_main("out", *self.write_sample_jsons([
            sample_with_coverage("XY", "Male", 15, 15), sample_with_coverage("XX", "Female", 30, 0)])))
        self.assertEqual({sample_id: (row["ChrXCoverageRatio"], row["ChrYCoverageRatio"])
                          for sample_id, row in sample_rows.items()},
                         {"XY": ("0.5", "0.5"), "XX": ("1", "0")})

    def test_read_count_metrics_and_fragment_length(self):
        sample = expansion_hunter_json("A", "Male", [
            # 1 of the 20 spanning reads is more than 1 repeat unit from the call (the read of size 11 is stutter)
            repeat_record("HOM", "chr1:1000-1030", "CAG", "10/10",
                          counts_of_spanning_reads="(10, 18), (11, 1), (14, 1)"),
            # left out: 30 CAG repeats are 90 bp, more than half of the 150 bp reads
            repeat_record("HOM_LONG", "chr1:2000-2030", "CAG", "30/30", counts_of_spanning_reads="(25, 5)"),
            # the longer allele has 8 of the 20 reads at the two allele sizes
            repeat_record("HET", "chr2:1000-1030", "CAG", "10/14", counts_of_spanning_reads="(10, 12), (14, 8)"),
            # left out: alleles 1 unit apart, whose reads include each other's stutter
            repeat_record("HET_ADJACENT", "chr2:2000-2030", "CAG", "10/11",
                          counts_of_spanning_reads="(10, 10), (11, 10)"),
            # left out: fewer than 10 reads at the two allele sizes
            repeat_record("HET_FEW_READS", "chr3:1000-1030", "CAG", "10/12",
                          counts_of_spanning_reads="(10, 3), (12, 2)"),
            # left out: not autosomal
            repeat_record("X1", "chrX:5000000-5000030", "CAG", "10", counts_of_spanning_reads="(10, 5), (20, 5)"),
        ], read_length=150)
        # ExpansionHunter writes a FragmentLength of 0 at loci without reads
        for locus_id, fragment_length in [("HOM", 300), ("HOM_LONG", 310), ("HET", 320), ("HET_ADJACENT", 330),
                                          ("HET_FEW_READS", 340), ("X1", 0)]:
            sample["LocusResults"][locus_id]["FragmentLength"] = fragment_length
        _, _, sample_rows = read_qc_stats_file(self.run_main("out", self.write_json("A", sample)))
        row = sample_rows["A"]
        self.assertEqual((row["FractionOfReadsDisagreeingWithHomozygousCalls"],
                          row["MeanLongerAlleleReadFractionAtHeterozygousCalls"], row["MedianFragmentLength"]),
                         ("0.05", "0.4", "320"))

    def test_call_rates_leave_out_female_samples_at_chrY_repeats(self):
        sample_jsons = [
            expansion_hunter_json("F", "Female", [
                repeat_record("L1", "chr1:100-130", "CAG", "10/10"),
                repeat_record("L2", "chr1:500-530", "CAG"),
                repeat_record("Y1", "chrY:100-130", "CAG"),
            ]),
            expansion_hunter_json("M1", "Male", [
                repeat_record("L1", "chr1:100-130", "CAG", "10/10"),
                repeat_record("L2", "chr1:500-530", "CAG", "10/10"),
                repeat_record("Y1", "chrY:100-130", "CAG", "10"),
            ]),
            expansion_hunter_json("M2", "Male", [
                repeat_record("L1", "chr1:100-130", "CAG", "10/10"),
                repeat_record("L2", "chr1:500-530", "CAG", "10/10"),
                repeat_record("Y1", "chrY:100-130", "CAG"),
            ]),
        ]
        _, locus_rows, sample_rows = read_qc_stats_file(self.run_main("out", *self.write_sample_jsons(sample_jsons)))
        self.assertEqual({sample_id: row["CallRate"] for sample_id, row in sample_rows.items()},
                         {"F": "0.5", "M1": "1", "M2": "0.66666667"})
        self.assertEqual({variant_id: (row["SampleCount"], row["CallRateSampleCount"],
                                       row["CallRateGenotypedSampleCount"])
                          for variant_id, row in locus_rows.items()},
                         {"L1": ("3", "3", "3"), "L2": ("3", "3", "2"), "Y1": ("3", "2", "1")})

    def test_combining_qc_stats_files_gives_the_same_file_as_combining_the_jsons(self):
        inputs = self.write_sample_jsons(SAMPLE_JSONS)
        expected = read_text(self.run_main("all", *inputs))
        first_half = self.run_main("first_half", inputs[0], inputs[2])
        second_half = self.run_main("second_half", inputs[1], inputs[3])
        self.assertEqual(read_text(self.run_main("tree", first_half, second_half)), expected)
        # QC stats files and JSONs can be combined together
        self.assertEqual(read_text(self.run_main("mixed", inputs[3], first_half, inputs[1])), expected)

    def test_threads_and_merge_rounds_give_the_same_file(self):
        # With at most 2 inputs per merge, the 5 inputs make 3 shards, whose files are merged in a round first
        inputs = self.write_sample_jsons(SAMPLE_JSONS + [
            {**SAMPLE_JSONS[0], "SampleParameters": {"SampleId": "E", "Sex": "Female"}}])
        expected = read_text(self.run_main("direct", *inputs))
        with mock.patch.object(streaming_utils, "get_max_inputs_per_merge", lambda: 2):
            self.assertEqual(read_text(self.run_main("in_rounds", *inputs)), expected)
            self.assertEqual(read_text(self.run_main("in_rounds_threaded", "--threads", 2, *inputs)), expected)

    def test_duplicate_sample_ids_are_rejected(self):
        inputs = self.write_sample_jsons(SAMPLE_JSONS)
        with self.assertRaisesRegex(ValueError, "both have SampleId A"):
            self.run_main("out", inputs[0], self.write_json("A_copy", SAMPLE_JSONS[0]))
        qc_stats_path = self.run_main("first_two", inputs[0], inputs[1])
        with self.assertRaisesRegex(ValueError, "both have SampleId B"):
            self.run_main("out", qc_stats_path, inputs[1])

    def test_samples_to_exclude(self):
        inputs = self.write_sample_jsons(SAMPLE_JSONS)
        samples_to_exclude = os.path.join(self.temp_dir, "samples_to_exclude.tsv")
        with open(samples_to_exclude, "w") as f:
            f.write("SampleId\tReasonsForExclusion\nB\tCall rate 0.5 < 0.9\nZ\tnot an input\n")
        sample_ids, _, sample_rows = read_qc_stats_file(
            self.run_main("out", "--samples-to-exclude", samples_to_exclude, *inputs))
        self.assertEqual((sample_ids, list(sample_rows)), (["A", "C", "D"], ["A", "C", "D"]))
        # a sample in a QC stats file input can't be left out, since the file's totals include it
        qc_stats_path = self.run_main("first_two", inputs[0], inputs[1])
        with self.assertRaisesRegex(ValueError, "Only JSON inputs can be left out"):
            self.run_main("out", "--samples-to-exclude", samples_to_exclude, qc_stats_path, inputs[2])

    def test_incomplete_or_other_files_are_rejected(self):
        inputs = self.write_sample_jsons(SAMPLE_JSONS)
        qc_stats_path = self.run_main("first_two", inputs[0], inputs[1])
        lines = read_text(qc_stats_path).splitlines(keepends=True)
        for name, text, error_text in [
                ("truncated_loci", "".join(lines[:4]), "ends before its sample rows"),
                ("truncated_samples", "".join(lines[:-1]), "fewer sample rows than SampleIds"),
                ("other", "LocusId\tVariantId\n", "isn't a QC stats file")]:
            path = os.path.join(self.temp_dir, name + m.QC_STATS_SUFFIX)
            with gzip.open(path, "wt") as f:
                f.write(text)
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, error_text):
                self.run_main("out", path, inputs[2])

    def test_jsons_from_every_analysis_mode_give_the_same_output(self):
        optimized_streaming_directory = os.path.join(self.temp_dir, "optimized_streaming")
        mixed_modes_directory = os.path.join(self.temp_dir, "mixed_modes")
        os.mkdir(optimized_streaming_directory)
        os.mkdir(mixed_modes_directory)
        expected_text = read_text(run_generate_multisample_qc_stats(
            optimized_streaming_directory, "out", *write_sample_jsons(optimized_streaming_directory, SAMPLE_JSONS)))
        # Seeking and streaming mode write the loci sorted by LocusId: here, reversed, so that L3 (chr2) comes first
        mixed_mode_jsons = json.loads(json.dumps(SAMPLE_JSONS))
        for sample, analysis_mode in zip(mixed_mode_jsons, ["seeking", "streaming", "optimized-streaming", "seeking"]):
            sample["RunInfo"]["AnalysisMode"] = analysis_mode
            if analysis_mode != "optimized-streaming":
                sample["LocusResults"] = dict(reversed(list(sample["LocusResults"].items())))
        for threads in (1, 2):
            with self.subTest(threads=threads):
                output_path = run_generate_multisample_qc_stats(
                    mixed_modes_directory, f"out_{threads}", "--threads", str(threads),
                    *write_sample_jsons(mixed_modes_directory, mixed_mode_jsons))
                # the sample rows' Path column names each JSON's directory
                self.assertEqual(read_text(output_path).replace(mixed_modes_directory, optimized_streaming_directory),
                                 expected_text)
                self.assertEqual([name for name in os.listdir(mixed_modes_directory) if ".tmp_" in name], [])

    def test_json_without_some_repeat_records_is_rejected(self):
        for option, error_text in [("--skip-missing-genotypes", "--skip-missing-genotypes"),
                                   ("--skip-hom-ref", "--skip-hom-ref"),
                                   ("--genotyping-approach only-quick", "--genotyping-approach only-quick"),
                                   ("--genotyping-approach=only-quick", "--genotyping-approach only-quick")]:
            sample = expansion_hunter_json("A", "Female", [repeat_record("L1", "chr1:100-130", "CAG", "10/10")])
            sample["RunInfo"]["CommandLine"] = f"ExpansionHunter --analysis-mode optimized-streaming {option}"
            # Older ExpansionHunter versions wrote RunInfo after the loci, so its command line is only read at the end
            older_layout = {key: sample[key] for key in ("SampleParameters", "LocusResults", "RunInfo")}
            for layout, data in [("current", sample), ("older", older_layout)]:
                with self.subTest(option=option, layout=layout), \
                        self.assertRaisesRegex(ValueError, f"written with ExpansionHunter's {error_text}"):
                    self.run_main("out", self.write_json("A", data))
        # The default genotyping approach writes every repeat's record, and a word that only starts with an option's
        # name (a made-up one here) isn't that option
        sample["RunInfo"]["CommandLine"] = "ExpansionHunter --genotyping-approach auto --skip-hom-ref-made-up-option"
        self.run_main("out", self.write_json("A", sample))

    def test_input_paths_file_and_default_output_name(self):
        inputs = self.write_sample_jsons(SAMPLE_JSONS)
        input_paths_file = os.path.join(self.temp_dir, "inputs.txt")
        with open(input_paths_file, "w") as f:
            f.write("# JSONs\n" + "\n".join(inputs[1:]) + "\n\n")
        current_directory = os.getcwd()
        os.chdir(self.temp_dir)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                m.main(["--input-paths-file", input_paths_file, inputs[0]])
        finally:
            os.chdir(current_directory)
        self.assertEqual(read_text(os.path.join(self.temp_dir, "combined.4_samples" + m.QC_STATS_SUFFIX)),
                         read_text(self.run_main("from_command_line", *inputs)))

    def test_compute_median_of_histogram(self):
        self.assertIsNone(m.compute_median_of_histogram({}))
        self.assertEqual(m.compute_median_of_histogram({5: 1}), 5)
        self.assertEqual(m.compute_median_of_histogram({1: 1, 2: 1, 3: 1}), 2)
        self.assertEqual(m.compute_median_of_histogram({1: 1, 2: 1}), 1.5)

    def test_combine_extreme_alleles(self):
        self.assertEqual(m.combine_extreme_alleles((10, 1, "A"), (12, 1, "B"), keep_larger=True), (12, 1, "B"))
        self.assertEqual(m.combine_extreme_alleles((10, 1, "A"), (12, 1, "B"), keep_larger=False), (10, 1, "A"))
        # a tie means no single sample has it
        self.assertEqual(m.combine_extreme_alleles((10, 1, "A"), (10, 2, ""), keep_larger=True), (10, 3, ""))
        self.assertEqual(m.combine_extreme_alleles(None, (10, 1, "A"), keep_larger=True), (10, 1, "A"))
        self.assertIsNone(m.combine_extreme_alleles(None, None, keep_larger=True))


if __name__ == "__main__":
    unittest.main()
