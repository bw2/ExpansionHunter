import contextlib
import gzip
import io
import json
import os
import tempfile
import unittest
from unittest import mock

import generate_multisample_tsv as m


def repeat_record(variant_id, region, motif, genotype=None, motif_composition=None, consensus_sequences=None):
    record = {"VariantId": variant_id, "VariantType": "Repeat", "ReferenceRegion": region, "RepeatUnit": motif}
    if genotype is not None:
        record["Genotype"] = genotype
    if motif_composition is not None:
        record["MotifComposition"] = motif_composition
    if consensus_sequences is not None:
        record["ConsensusSequences"] = consensus_sequences
    return record


def expansion_hunter_json(sample_id, sex, loci):
    """Returns a JSON laid out like the optimized-streaming output of older ExpansionHunter versions, with RunInfo
    after the loci. loci: list of (locus_id, [variant records])."""
    return {
        "SampleParameters": {"SampleId": sample_id, "Sex": sex},
        "LocusResults": {
            locus_id: {"LocusId": locus_id, "Variants": {record["VariantId"]: record for record in records}}
            for locus_id, records in loci
        },
        "RunInfo": {"AnalysisMode": "optimized-streaming"},
    }


def counts(count, reads):
    return {"count": count, "reads": reads}


def in_older_seeking_mode_layout(data):
    """Returns the JSON the way older ExpansionHunter versions' seeking mode wrote it, with AnalysisMode set to seeking
    and the top-level keys sorted (nlohmann::json sorts them): LocusResults (sorted by LocusId), then RunInfo, then
    SampleParameters."""
    return {
        "LocusResults": dict(sorted(data["LocusResults"].items())),
        "RunInfo": {**data["RunInfo"], "AnalysisMode": "seeking"},
        "SampleParameters": data["SampleParameters"],
    }


def in_current_layout(data, analysis_mode="optimized-streaming"):
    """Returns the JSON the way ExpansionHunter now writes it in every mode: SampleParameters, RunInfo (with
    AnalysisMode set to analysis_mode), LocusResults (sorted by LocusId unless in optimized-streaming mode), then
    RunStats."""
    locus_results = data["LocusResults"]
    if analysis_mode != "optimized-streaming":
        locus_results = dict(sorted(locus_results.items()))
    return {
        "SampleParameters": data["SampleParameters"],
        "RunInfo": {**data["RunInfo"], "AnalysisMode": analysis_mode},
        "LocusResults": locus_results,
        "RunStats": {"Started": "2026-11-03 18:00:00", "PeakRssMemoryMb": 100.0},
    }


def with_sample_id(data, sample_id):
    copy = json.loads(json.dumps(data))
    copy["SampleParameters"]["SampleId"] = sample_id
    return copy


SAMPLE_A = expansion_hunter_json("A", "Male", [
    ("L1", [repeat_record("L1", "chr1:100-130", "CAG", "10/12",
                          {"Motifs": {"1:CAG": counts(30, 5), "2:CAA": counts(2, 2)},
                           "MotifPairs": {"[1][1]": counts(25, 5), "[1][2]": counts(2, 2)}},
                          ["CAG" * 10, "CAG" * 11 + "CAA"])]),
    ("M", [repeat_record("M_CAG", "chr2:200-230", "CAG", "10/10"),
           repeat_record("M_CCG", "chr2:233-260", "CCG", "9/9"),
           {"VariantId": "M_SNV", "VariantType": "SmallVariant", "ReferenceRegion": "chr2:231-232"}]),
    ("L2", [repeat_record("L2", "chrX:500-530", "GCC", "8",
                          {"Motifs": {"1:GCC": counts(8, 3)}, "MotifPairs": {"[1][1]": counts(7, 3)}},
                          ["GCC" * 8])]),
])
SAMPLE_B = expansion_hunter_json("B", "Female", [
    ("L1", [repeat_record("L1", "chr1:100-130", "CAG", "12/12",
                          {"Motifs": {"1:CAA": counts(40, 6), "2:CAG": counts(10, 4)},
                           "MotifPairs": {"[1][1]": counts(30, 6), "[2][1]": counts(5, 3)}},
                          ["CAG" * 11 + "CAA"])]),
    ("M", [repeat_record("M_CAG", "chr2:200-230", "CAG", "10/11"),
           repeat_record("M_CCG", "chr2:233-260", "CCG")]),
    ("L2", [repeat_record("L2", "chrX:500-530", "GCC", "8/9",
                          {"Motifs": {"1:GCC": counts(17, 6)}, "MotifPairs": {"[1][1]": counts(15, 6)}},
                          ["GCC" * 8, "GCC" * 9])]),
])
SAMPLE_C = expansion_hunter_json("C", "Male", [
    ("L1", [repeat_record("L1", "chr1:100-130", "CAG", "11/40",
                          {"Motifs": {"1:CAG": counts(50, 7)}, "MotifPairs": {"[1][1]": counts(45, 7)}},
                          ["CAG" * 11, "CAG" * 40])]),
    ("M", [repeat_record("M_CAG", "chr2:200-230", "CAG", "10/10"),
           repeat_record("M_CCG", "chr2:233-260", "CCG", "9/10")]),
    ("L2", [repeat_record("L2", "chrX:500-530", "GCC", "20",
                          {"Motifs": {}, "MotifPairs": {}}, ["GCC" * 20])]),
])
SAMPLE_WITHOUT_EXTRA_FIELDS = expansion_hunter_json("D", "Female", [
    ("L1", [repeat_record("L1", "chr1:100-130", "CAG", "10/10")]),
    ("M", [repeat_record("M_CAG", "chr2:200-230", "CAG", "10/10"),
           repeat_record("M_CCG", "chr2:233-260", "CCG", "9/9")]),
    ("L2", [repeat_record("L2", "chrX:500-530", "GCC", "8/8")]),
])
# Only the second repeat record has MotifComposition, as in ExpansionHunter's sparse motif composition modes.
SAMPLE_WITH_MOTIF_COMPOSITION_ONLY_IN_SECOND_RECORD = expansion_hunter_json("A", "Female", [
    ("first", [repeat_record("first", "chr1:100-130", "A", "30/30")]),
    ("second", [repeat_record("second", "chr1:200-230", "CA", "15/15",
                              {"Motifs": {"1:CA": counts(15, 2)}, "MotifPairs": {}})]),
])


class GenerateMultisampleTsvTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir_object = tempfile.TemporaryDirectory()
        self.temp_dir = self.temp_dir_object.name

    def tearDown(self):
        self.temp_dir_object.cleanup()

    def write_json(self, name, data, compress=True):
        path = os.path.join(self.temp_dir, name + (".json.gz" if compress else ".json"))
        with (gzip.open(path, "wt") if compress else open(path, "w")) as f:
            json.dump(data, f)
        return path

    def run_main_and_capture_output(self, output_name, *args):
        """Runs main and returns (output prefix, what it printed)."""
        output_prefix = os.path.join(self.temp_dir, output_name)
        with contextlib.redirect_stdout(io.StringIO()) as stdout:
            m.main(["-o", output_prefix] + [str(arg) for arg in args])
        return output_prefix, stdout.getvalue()

    def run_main(self, output_name, *args):
        return self.run_main_and_capture_output(output_name, *args)[0]

    def read_output(self, output_prefix):
        with gzip.open(output_prefix + ".tsv.gz", "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in f]
        return header, {row["VariantId"]: row for row in rows}

    def read_output_text(self, output_prefix):
        with gzip.open(output_prefix + ".tsv.gz", "rt") as f:
            return f.read()

    def test_histograms_and_counts(self):
        output_prefix = self.run_main(
            "out", self.write_json("A", SAMPLE_A), self.write_json("B", SAMPLE_B, compress=False))
        header, rows = self.read_output(output_prefix)
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(list(rows), ["L1", "M_CAG", "M_CCG", "L2"])

        self.assertEqual(rows["L1"]["SampleCount"], "2")
        self.assertEqual(rows["L1"]["CalledAlleleCount"], "4")
        self.assertEqual(rows["L1"]["AllAlleleHistogram"], "10x:1,12x:3")
        self.assertEqual(rows["L1"]["ShortAlleleHistogram"], "10x:1,12x:1")
        self.assertEqual(rows["L1"]["HemizygousAlleleHistogram"], "")
        self.assertEqual(rows["L1"]["BiallelicHistogram"], "10/12:1,12/12:1")
        self.assertEqual(rows["L1"]["LocusStructure"], "(CAG)*")

        # B has no genotype at M_CCG: it is counted as a sample but adds no called alleles
        self.assertEqual(rows["M_CCG"]["SampleCount"], "2")
        self.assertEqual(rows["M_CCG"]["CalledAlleleCount"], "2")
        self.assertEqual(rows["M_CCG"]["AllAlleleHistogram"], "9x:2")
        self.assertEqual(rows["M_CCG"]["BiallelicHistogram"], "9/9:1")
        self.assertEqual(rows["M_CAG"]["LocusStructure"], "")

        # A is male, so its one-allele chrX call is hemizygous; B is female and diploid there
        self.assertEqual(rows["L2"]["CalledAlleleCount"], "3")
        self.assertEqual(rows["L2"]["AllAlleleHistogram"], "8x:2,9x:1")
        self.assertEqual(rows["L2"]["ShortAlleleHistogram"], "8x:2")
        self.assertEqual(rows["L2"]["HemizygousAlleleHistogram"], "8x:1")
        # A's hemizygous 8 is written as 8/8
        self.assertEqual(rows["L2"]["BiallelicHistogram"], "8/8:1,8/9:1")
        with open(output_prefix + ".samples.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tSex\nA\tMale\nB\tFemale\n")

    def test_motif_composition_is_summed_by_sequence_and_renumbered(self):
        output_prefix = self.run_main("out", self.write_json("A", SAMPLE_A), self.write_json("B", SAMPLE_B))
        _, rows = self.read_output(output_prefix)
        # CAG: 30 + 10 = 40 occurrences, CAA: 2 + 40 = 42, so CAA becomes motif 1 although A called it motif 2.
        self.assertEqual(rows["L1"]["MotifComposition"], "1:CAA:42:8,2:CAG:40:9")
        # A: CAG-CAG 25, CAG-CAA 2; B: CAA-CAA 30, CAG-CAA 5 (its [2][1])
        self.assertEqual(rows["L1"]["MotifPairsComposition"], "[1][1]:30:6,[2][1]:7:5,[2][2]:25:5")
        self.assertEqual(rows["M_CAG"]["MotifComposition"], "")

    def test_allele_sequences_are_haplotype_counts(self):
        output_prefix = self.run_main("out", "--add-consensus-sequence-columns",
                                      self.write_json("A", SAMPLE_A), self.write_json("B", SAMPLE_B))
        header, rows = self.read_output(output_prefix)
        self.assertEqual(header[-1], m.ALLELE_SEQUENCES_COLUMN)
        # B's homozygous 12/12 call has one consensus sequence, which counts twice
        self.assertEqual(rows["L1"]["AlleleSequences"], f"{'CAG' * 11 + 'CAA'}:3,{'CAG' * 10}:1")
        self.assertEqual(rows["L2"]["AlleleSequences"], f"{'GCC' * 8}:2,{'GCC' * 9}:1")

    def test_consensus_sequence_columns_need_consensus_sequences_in_every_input(self):
        with self.assertRaisesRegex(ValueError, "ConsensusSequences"):
            self.run_main("out", "--add-consensus-sequence-columns", "--skip-motif-composition-columns",
                          self.write_json("A", SAMPLE_A), self.write_json("D", SAMPLE_WITHOUT_EXTRA_FIELDS))

    def test_inputs_that_differ_in_motif_composition_need_the_skip_option(self):
        inputs = [self.write_json("A", SAMPLE_A), self.write_json("D", SAMPLE_WITHOUT_EXTRA_FIELDS)]
        with self.assertRaises(ValueError):
            self.run_main("out", *inputs)
        self.assertFalse(os.path.exists(os.path.join(self.temp_dir, "out.tsv.gz")))
        header, rows = self.read_output(self.run_main("out", "--skip-motif-composition-columns", *inputs))
        self.assertEqual(header, m.BASE_COLUMNS)
        self.assertEqual(rows["L1"]["AllAlleleHistogram"], "10x:3,12x:1")

    def test_motif_composition_found_only_late_in_a_file_reruns_the_merge_with_its_columns(self):
        # The sample's RunInfo has no CommandLine, so only the full read finds its motif composition
        with mock.patch.object(m, "MOTIF_COMPOSITION_CHECK_RECORDS", 1):
            output_prefix, printed = self.run_main_and_capture_output(
                "out", self.write_json("A", SAMPLE_WITH_MOTIF_COMPOSITION_ONLY_IN_SECOND_RECORD))
        self.assertIn("merging again", printed)
        header, rows = self.read_output(output_prefix)
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(rows["second"]["MotifComposition"], "1:CA:15:2")
        self.assertEqual(rows["first"]["MotifComposition"], "")

    def test_command_line_of_a_sparse_motif_composition_mode_avoids_rerunning_the_merge(self):
        sample = with_sample_id(SAMPLE_WITH_MOTIF_COMPOSITION_ONLY_IN_SECOND_RECORD, "A")
        sample["RunInfo"]["CommandLine"] = "ExpansionHunter --output-motif-composition loci-with-known-motifs"
        with mock.patch.object(m, "MOTIF_COMPOSITION_CHECK_RECORDS", 1):
            output_prefix, printed = self.run_main_and_capture_output("out", self.write_json("A", sample))
        self.assertNotIn("merging again", printed)
        header, rows = self.read_output(output_prefix)
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(rows["second"]["MotifComposition"], "1:CA:15:2")

    def test_command_line_decides_whether_a_json_has_motif_composition(self):
        # In loci-with-non-ref-motifs mode, a sample without non-reference motifs has no MotifComposition record
        sample_without_records = with_sample_id(SAMPLE_WITHOUT_EXTRA_FIELDS, "D")
        sample_without_records["RunInfo"]["CommandLine"] = "ExpansionHunter --output-motif-composition loci-with-non-ref-motifs"
        header, rows = self.read_output(
            self.run_main("out", self.write_json("A", SAMPLE_A), self.write_json("D", sample_without_records)))
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(rows["L1"]["MotifComposition"], "1:CAG:30:5,2:CAA:2:2")

        sample_run_without_motif_composition = with_sample_id(SAMPLE_WITHOUT_EXTRA_FIELDS, "E")
        sample_run_without_motif_composition["RunInfo"]["CommandLine"] = "ExpansionHunter --analysis-mode seeking"
        with self.assertRaises(ValueError):
            self.run_main("out2", self.write_json("D", sample_without_records),
                          self.write_json("E", sample_run_without_motif_composition))

    def test_motif_composition_records_count_even_if_the_command_line_abbreviates_the_option(self):
        # A has motif composition only through its records, D only through its command line; if A's records didn't
        # count, the two would differ and main would raise.
        sample_with_abbreviated_option = with_sample_id(SAMPLE_A, "A")
        sample_with_abbreviated_option["RunInfo"]["CommandLine"] = "ExpansionHunter --output-motif-comp all-loci"
        sample_without_records = with_sample_id(SAMPLE_WITHOUT_EXTRA_FIELDS, "D")
        sample_without_records["RunInfo"]["CommandLine"] = "ExpansionHunter --output-motif-composition loci-with-non-ref-motifs"
        header, _ = self.read_output(self.run_main("out", self.write_json("A", sample_with_abbreviated_option),
                                                   self.write_json("D", sample_without_records)))
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)

    def test_a_command_line_of_hundreds_of_kilobytes_is_still_read(self):
        long_command_line = ("ExpansionHunter --output-motif-composition loci-with-non-ref-motifs --locus "
                             + ",".join(f"chr1-{index}-{index + 30}-CAG" for index in range(20_000)))
        self.assertGreater(len(long_command_line), 300_000)
        sample_with_records = with_sample_id(SAMPLE_A, "A")
        sample_with_records["RunInfo"]["CommandLine"] = long_command_line
        # In seeking-mode layout, SampleParameters comes after RunInfo
        seeking_sample_without_records = in_older_seeking_mode_layout(with_sample_id(SAMPLE_WITHOUT_EXTRA_FIELDS, "D"))
        seeking_sample_without_records["RunInfo"]["CommandLine"] = long_command_line
        header, rows = self.read_output(self.run_main(
            "out", self.write_json("A", sample_with_records), self.write_json("D", seeking_sample_without_records)))
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(rows["L1"]["MotifComposition"], "1:CAG:30:5,2:CAA:2:2")

    def test_run_info_keeping_reader_finds_a_key_split_between_reads(self):
        json_bytes = json.dumps({"LocusResults": {"L": {}}, "RunInfo": {"CommandLine": "x" * 1000}}).encode()
        reader = m.RunInfoKeepingReader(io.BytesIO(json_bytes))
        while reader.read(3):
            pass
        self.assertEqual(m.parse_run_info(reader.tail), {"CommandLine": "x" * 1000})
        with mock.patch.object(m, "MAX_RUN_INFO_BYTES", 100):
            reader = m.RunInfoKeepingReader(io.BytesIO(json_bytes))
            while reader.read(3):
                pass
        self.assertIsNone(m.parse_run_info(reader.tail))

    def test_consensus_sequence_columns_accept_an_input_without_any_called_repeat(self):
        uncalled = expansion_hunter_json("U", "Female", [("L1", [repeat_record("L1", "chr1:100-130", "CAG")])])
        _, rows = self.read_output(self.run_main(
            "out", "--add-consensus-sequence-columns", "--skip-motif-composition-columns",
            self.write_json("A", SAMPLE_A), self.write_json("U", uncalled)))
        self.assertEqual(rows["L1"]["SampleCount"], "2")
        self.assertEqual(rows["L1"]["AlleleSequences"], f"{'CAG' * 10}:1,{'CAG' * 11 + 'CAA'}:1")

    def test_no_motif_composition_columns_when_no_input_has_them(self):
        header, _ = self.read_output(self.run_main("out", self.write_json("D", SAMPLE_WITHOUT_EXTRA_FIELDS)))
        self.assertEqual(header, m.BASE_COLUMNS)

    def test_outlier_columns(self):
        inputs = [self.write_json(name, data) for name, data in [("A", SAMPLE_A), ("B", SAMPLE_B), ("C", SAMPLE_C)]]
        header, rows = self.read_output(self.run_main("out", *inputs))
        self.assertFalse(set(m.OUTLIER_COLUMNS) & set(header))
        header, rows = self.read_output(
            self.run_main("out", "--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "2", *inputs))
        self.assertTrue(set(m.OUTLIER_COLUMNS) <= set(header))
        # largest first; 12x has 2 samples, which is n, so the walk stops there
        self.assertEqual(rows["L1"]["OutlierSampleIds_AllAlleles"], "40x:C")
        self.assertEqual(rows["L2"]["OutlierSampleIds_HemizygousAlleles"], "20x:C,8x:A")

    def test_tsvs_from_earlier_runs_give_the_same_output_as_the_json_files(self):
        inputs = [self.write_json(name, data) for name, data in [("A", SAMPLE_A), ("B", SAMPLE_B), ("C", SAMPLE_C)]]
        options = ["--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "2", "--add-consensus-sequence-columns"]
        direct = self.run_main("direct", *options, *inputs)
        batch = self.run_main("batch", *options, inputs[0], inputs[1])
        combined = self.run_main("combined", *options, batch + ".tsv.gz", inputs[2])
        self.assertEqual(self.read_output_text(direct), self.read_output_text(combined))
        # 12x at L1 is in 2 samples of the batch, so the batch lists nobody there; the combined walk has to recognize
        # it as common from the histogram alone
        _, batch_rows = self.read_output(batch)
        self.assertEqual(batch_rows["L1"]["OutlierSampleIds_AllAlleles"], "")
        _, combined_rows = self.read_output(combined)
        self.assertEqual(combined_rows["L1"]["OutlierSampleIds_AllAlleles"], "40x:C")
        self.assertEqual(combined_rows["L2"]["OutlierSampleIds_HemizygousAlleles"], "20x:C,8x:A")
        with open(combined + ".samples.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tSex\nA\tMale\nB\tFemale\nC\tMale\n")

    def test_threads_give_the_same_output(self):
        inputs = [self.write_json(name, data) for name, data in [("A", SAMPLE_A), ("B", SAMPLE_B), ("C", SAMPLE_C)]]
        options = ["--add-outlier-sample-id-columns", "--add-consensus-sequence-columns"]
        single = self.run_main("single", *options, *inputs)
        threaded = self.run_main("threaded", "--threads", "2", *options, *inputs)
        self.assertEqual(self.read_output_text(single), self.read_output_text(threaded))
        self.assertEqual([name for name in os.listdir(self.temp_dir) if ".tmp_" in name], [])

    def test_older_seeking_mode_layout_gives_the_same_output(self):
        options = ["--add-outlier-sample-id-columns", "--add-consensus-sequence-columns"]
        optimized_streaming = self.run_main(
            "optimized_streaming", *options, self.write_json("A", SAMPLE_A), self.write_json("B", SAMPLE_B))
        # sorted by LocusId, L2 (chrX) comes before M (chr2)
        seeking = self.run_main(
            "seeking", *options, self.write_json("A_seeking", in_older_seeking_mode_layout(SAMPLE_A)),
            self.write_json("B_seeking", in_older_seeking_mode_layout(SAMPLE_B), compress=False))
        self.assertEqual(self.read_output_text(optimized_streaming), self.read_output_text(seeking))
        with open(seeking + ".samples.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tSex\nA\tMale\nB\tFemale\n")

    def test_current_layout_gives_the_same_output_as_older_layouts(self):
        options = ["--add-outlier-sample-id-columns", "--add-consensus-sequence-columns"]
        older = self.run_main("older", *options, self.write_json("A", SAMPLE_A),
                              self.write_json("B", in_older_seeking_mode_layout(SAMPLE_B)))
        # In seeking mode, sorted by LocusId, L2 (chrX) comes before M (chr2), so only RunInfo.AnalysisMode tells
        # that B's loci have to be sorted.
        current = self.run_main(
            "current", *options, self.write_json("A_current", in_current_layout(SAMPLE_A)),
            self.write_json("B_current", in_current_layout(SAMPLE_B, "seeking"), compress=False))
        self.assertEqual(self.read_output_text(older), self.read_output_text(current))
        with open(current + ".samples.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tSex\nA\tMale\nB\tFemale\n")

    def test_current_layout_command_line_is_read_before_the_merge_without_keeping_the_tail(self):
        sample = with_sample_id(SAMPLE_WITH_MOTIF_COMPOSITION_ONLY_IN_SECOND_RECORD, "A")
        sample["RunInfo"]["CommandLine"] = "ExpansionHunter --output-motif-composition loci-with-known-motifs"
        with mock.patch.object(m, "MOTIF_COMPOSITION_CHECK_RECORDS", 1), \
                mock.patch.object(m, "RunInfoKeepingReader") as run_info_keeping_reader:
            output_prefix, printed = self.run_main_and_capture_output(
                "out", self.write_json("A", in_current_layout(sample)))
        run_info_keeping_reader.assert_not_called()
        self.assertNotIn("merging again", printed)
        header, rows = self.read_output(output_prefix)
        self.assertEqual(header, m.BASE_COLUMNS + m.MOTIF_COMPOSITION_COLUMNS)
        self.assertEqual(rows["second"]["MotifComposition"], "1:CA:15:2")

    def test_many_inputs_are_merged_in_rounds_that_stay_under_the_open_file_limit(self):
        inputs = [self.write_json(f"S{index}", with_sample_id(data, f"S{index}"))
                  for index, data in enumerate([SAMPLE_A, SAMPLE_B, SAMPLE_C, SAMPLE_A, SAMPLE_B])]
        options = ["--add-outlier-sample-id-columns", "--add-consensus-sequence-columns"]
        direct = self.run_main("direct", *options, *inputs)
        with mock.patch.object(m, "get_max_inputs_per_merge", lambda: 2):
            in_rounds = self.run_main("in_rounds", *options, *inputs)
            in_rounds_threaded = self.run_main("in_rounds_threaded", "--threads", "2", *options, *inputs)
        self.assertEqual(self.read_output_text(direct), self.read_output_text(in_rounds))
        self.assertEqual(self.read_output_text(direct), self.read_output_text(in_rounds_threaded))

    def test_failing_worker_raises_its_error_and_leaves_no_temporary_files(self):
        bad_sample = expansion_hunter_json("bad", "Female", [
            ("on_chr2", [repeat_record("on_chr2", "chr2:100-130", "CA", "5/5")]),
            ("on_chr1", [repeat_record("on_chr1", "chr1:100-130", "CA", "6/6")]),
        ])
        with self.assertRaisesRegex(ValueError, "lists chr1 after a chromosome"):
            self.run_main("out", "--threads", "2", "--skip-motif-composition-columns", self.write_json("bad", bad_sample),
                          self.write_json("A", SAMPLE_A))
        self.assertEqual([name for name in os.listdir(self.temp_dir) if ".tmp_" in name], [])

    def test_tsv_made_with_a_smaller_outlier_limit_raises(self):
        inputs = [self.write_json(name, data) for name, data in [("A", SAMPLE_A), ("B", SAMPLE_B), ("C", SAMPLE_C)]]
        batch = self.run_main("batch", "--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "2",
                              inputs[0], inputs[1])
        with self.assertRaises(ValueError):
            self.run_main("out", "--add-outlier-sample-id-columns", batch + ".tsv.gz", inputs[2])

    def test_outlier_limit_below_2_is_rejected(self):
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.run_main("out", "--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "1",
                          self.write_json("A", SAMPLE_A))

    def test_tsv_made_with_a_larger_outlier_limit_gives_the_same_output(self):
        inputs = [self.write_json(name, data) for name, data in [("A", SAMPLE_A), ("B", SAMPLE_B), ("C", SAMPLE_C)]]
        direct = self.run_main("direct", "--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "2", *inputs)
        batch = self.run_main("batch", "--add-outlier-sample-id-columns", inputs[0], inputs[1])
        combined = self.run_main("combined", "--add-outlier-sample-id-columns", "--n-outlier-sample-ids", "2",
                                 batch + ".tsv.gz", inputs[2])
        self.assertEqual(self.read_output_text(direct), self.read_output_text(combined))
        _, rows = self.read_output(combined)
        self.assertEqual(rows["L1"]["OutlierSampleIds_AllAlleles"], "40x:C")

    def test_tsv_input_needs_its_samples_table(self):
        batch = self.run_main("batch", self.write_json("A", SAMPLE_A))
        os.remove(batch + ".samples.tsv")
        for extra_options in [[], ["--allow-overlapping-samples"]]:
            with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
                self.run_main("out", *extra_options, batch + ".tsv.gz", self.write_json("B", SAMPLE_B))

    def test_tsv_input_without_allele_sequences_raises_when_they_are_requested(self):
        batch = self.run_main("batch", self.write_json("A", SAMPLE_A))
        with self.assertRaises(ValueError):
            self.run_main("out", "--add-consensus-sequence-columns", batch + ".tsv.gz", self.write_json("B", SAMPLE_B))

    def test_sample_in_two_inputs_raises(self):
        inputs = [self.write_json("A", SAMPLE_A), self.write_json("A_again", SAMPLE_A)]
        with self.assertRaises(ValueError):
            self.run_main("out", *inputs)
        batch = self.run_main("batch", inputs[0])
        with self.assertRaises(ValueError):
            self.run_main("out", batch + ".tsv.gz", inputs[1])

    def test_tsv_input_without_outlier_columns_raises_when_they_are_requested(self):
        batch = self.run_main("batch", self.write_json("A", SAMPLE_A))
        with self.assertRaises(ValueError):
            self.run_main("out", "--add-outlier-sample-id-columns", batch + ".tsv.gz", self.write_json("B", SAMPLE_B))

    def test_variant_catalog_fills_in_locus_structure(self):
        catalog_path = os.path.join(self.temp_dir, "catalog.json")
        with open(catalog_path, "w") as f:
            json.dump([{"LocusId": "L1", "LocusStructure": "(CAG)*"},
                       {"LocusId": "M", "LocusStructure": "(CAG)*AA(CCG)*"}], f)
        _, rows = self.read_output(
            self.run_main("out", "--variant-catalog", catalog_path, self.write_json("A", SAMPLE_A)))
        self.assertEqual(rows["M_CAG"]["LocusStructure"], "(CAG)*AA(CCG)*")
        self.assertEqual(rows["M_CCG"]["LocusStructure"], "(CAG)*AA(CCG)*")
        self.assertEqual(rows["L1"]["LocusStructure"], "(CAG)*")

    def test_records_slightly_out_of_order_are_put_in_order(self):
        sample = expansion_hunter_json("A", "Female", [
            ("late", [repeat_record("late", "chr1:900-950", "CA", "5/5")]),
            ("early", [repeat_record("early", "chr1:100-130", "CA", "6/6")]),
        ])
        _, rows = self.read_output(self.run_main("out", self.write_json("A", sample)))
        self.assertEqual(list(rows), ["early", "late"])

    def test_records_far_out_of_order_raise(self):
        # "first" is written once "far" is more than the reorder window past it, so "earliest" comes too late
        sample = expansion_hunter_json("A", "Female", [
            ("first", [repeat_record("first", "chr1:1000-1030", "CA", "5/5")]),
            ("far", [repeat_record("far", "chr1:900000-900050", "CA", "5/5")]),
            ("earliest", [repeat_record("earliest", "chr1:100-130", "CA", "6/6")]),
        ])
        with self.assertRaises(ValueError):
            self.run_main("out", self.write_json("A", sample))
        sample = expansion_hunter_json("A", "Female", [
            ("on_chr2", [repeat_record("on_chr2", "chr2:100-130", "CA", "5/5")]),
            ("on_chr1", [repeat_record("on_chr1", "chr1:100-130", "CA", "6/6")]),
        ])
        with self.assertRaises(ValueError):
            self.run_main("out", self.write_json("B", sample))


class FormattingTests(unittest.TestCase):

    def test_parse_genotype(self):
        self.assertEqual(m.parse_genotype("15/17"), [15, 17])
        self.assertEqual(m.parse_genotype("15"), [15])
        self.assertIsNone(m.parse_genotype(None))
        self.assertIsNone(m.parse_genotype("./."))

    def test_histogram_round_trip(self):
        self.assertEqual(m.format_histogram({16: 3, 15: 120}), "15x:120,16x:3")
        self.assertEqual(m.parse_histogram("15x:120,16x:3"), {15: 120, 16: 3})
        self.assertEqual(m.parse_histogram(""), {})

    def test_biallelic_histogram_round_trip(self):
        self.assertEqual(m.format_biallelic_histogram({(15, 16): 3, (15, 15): 58, (8, 30): 1}), "8/30:1,15/15:58,15/16:3")
        self.assertEqual(m.parse_biallelic_histogram("8/30:1,15/15:58,15/16:3"), {(8, 30): 1, (15, 15): 58, (15, 16): 3})
        self.assertEqual(m.parse_biallelic_histogram(""), {})

    def test_motif_composition_round_trip(self):
        motif_counts = {"GCG": [561, 64], "GCA": [82, 50]}
        pair_counts = {("GCG", "GCG"): [398, 63], ("GCG", "GCA"): [82, 50]}
        motifs_value, pairs_value = m.format_motif_composition(motif_counts, pair_counts)
        self.assertEqual(motifs_value, "1:GCG:561:64,2:GCA:82:50")
        self.assertEqual(pairs_value, "[1][1]:398:63,[1][2]:82:50")
        self.assertEqual(m.parse_motif_composition_columns(motifs_value, pairs_value), (motif_counts, pair_counts))
        self.assertEqual(m.format_motif_composition({"AT": [5, 1], "AC": [5, 2]}, {})[0], "1:AC:5:2,2:AT:5:1")

    def test_allele_sequences_round_trip(self):
        value = m.format_allele_sequences({"TGCTGA": 56, "TGCTGC": 86})
        self.assertEqual(value, "TGCTGC:86,TGCTGA:56")
        self.assertEqual(m.parse_allele_sequences(value), {"TGCTGA": 56, "TGCTGC": 86})

    def test_outlier_walk_stops_at_common_allele_sizes(self):
        histogram = {10: 40, 30: 1, 31: 2}
        outlier_sample_ids = {10: ["s1", "s2"], 30: ["s3"], 31: ["s4"]}
        # 10x has 40 alleles but only 2 listed samples, so it is common in inputs that didn't list it
        self.assertEqual(m.get_outlier_sample_ids_string(outlier_sample_ids, histogram, 3), "31x:s4,30x:s3")

    def test_genomic_sort_key_orders_chromosomes_naturally(self):
        regions = ["chrX:1-2", "chr10:1-2", "chr2:5-6", "chr2:1-2", "chrM:1-2", "chrUn_x:1-2", "chrY:1-2"]
        ordered = sorted(regions, key=lambda region: m.get_genomic_sort_key("id", "id", region))
        self.assertEqual(ordered, ["chr2:1-2", "chr2:5-6", "chr10:1-2", "chrX:1-2", "chrY:1-2", "chrM:1-2",
                                   "chrUn_x:1-2"])


if __name__ == "__main__":
    unittest.main()
