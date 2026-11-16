import contextlib
import gzip
import http.client
import io
import json
import math
import os
import tempfile
import threading
import unittest
from unittest import mock

import numpy as np

import start_multisample_qc_dashboard as m
from generate_multisample_qc_stats_tests import (
    SAMPLE_JSONS, expansion_hunter_json, repeat_record, run_generate_multisample_qc_stats, write_sample_jsons)


def count_expansions_and_contractions_at_repeats_called_in(min_genotyped_samples, min_fraction_of_samples=0):
    """Returns a patch of the thresholds that decide which repeats LociWithLargestExpansion and
    LociWithLargestContraction count, since the test cohorts are too small for the dashboard's."""
    return mock.patch.multiple(
        m, MIN_GENOTYPED_SAMPLES_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS=min_genotyped_samples,
        MIN_FRACTION_OF_SAMPLES_GENOTYPED_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS=min_fraction_of_samples)


class StartMultisampleQcDashboardTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir_object = tempfile.TemporaryDirectory()
        self.temp_dir = self.temp_dir_object.name

    def tearDown(self):
        self.temp_dir_object.cleanup()

    def write_qc_stats_file(self, sample_jsons):
        """Writes the JSONs and combines them into a QC stats file, and returns its path."""
        return run_generate_multisample_qc_stats(self.temp_dir, "cohort",
                                                 *write_sample_jsons(self.temp_dir, sample_jsons))

    def test_locus_metrics_and_expansion_counts(self):
        with count_expansions_and_contractions_at_repeats_called_in(1):
            qc_stats = m.read_qc_stats(self.write_qc_stats_file(SAMPLE_JSONS), 1)
        self.assertEqual([(repeat.locus_id, repeat.reference_region, repeat.genotyped_sample_count, repeat.call_rate)
                          for repeat in qc_stats.repeats],
                         [("L1", "chr1:100-130", 4, 1), ("L2", "chr1:500-530", 3, 0.75), ("L3", "chr2:100-120", 4, 1)])
        self.assertEqual(qc_stats.repeat_counts_of_left_out_compound_loci, {})
        # L1: 2 of 4 calls heterozygous, and allele frequencies 6/8, 1/8 and 1/8 give an expected heterozygosity of
        # 1 - (0.5625 + 0.015625 + 0.015625); L2: 1 of 3, frequencies 5/6 and 1/6; L3: 1 of 4, frequencies 5/8 and 3/8
        self.assertEqual([round(repeat.observed_vs_expected_heterozygosity, 3) for repeat in qc_stats.repeats],
                         [1.231, 1.2, 0.533])
        # D has L1's largest expansion and C L2's largest contraction. At L3, A and B share the largest allele.
        self.assertEqual({row["SampleId"]: (row["LociWithLargestExpansion"], row["LociWithLargestContraction"],
                                            row["ExpansionToContractionRatio"]) for row in qc_stats.sample_rows},
                         {"A": (0, 0, 1), "B": (0, 0, 1), "C": (0, 1, 0.5), "D": (1, 0, 2)})

    def test_qc_stats_file_without_repeats(self):
        qc_stats = m.read_qc_stats(self.write_qc_stats_file([expansion_hunter_json("A", "Female", [])]), 1)
        self.assertEqual((qc_stats.repeats, [row["SampleId"] for row in qc_stats.sample_rows]), ([], ["A"]))

    def test_observed_vs_expected_heterozygosity_needs_enough_two_allele_calls(self):
        qc_stats = m.read_qc_stats(self.write_qc_stats_file(SAMPLE_JSONS), 4)
        # L2 is called in only 3 samples
        self.assertEqual([repeat.observed_vs_expected_heterozygosity is None for repeat in qc_stats.repeats],
                         [False, True, False])
        # a chrX repeat called in 3 male samples, with one allele, and in 2 female ones, with two: only the 2 count
        sample_jsons = [expansion_hunter_json(f"M{index}", "Male", [
            repeat_record("X1", "chrX:5000000-5000030", "CAG", genotype)]) for index, genotype in enumerate(
            ["10", "11", "12"])] + [expansion_hunter_json(f"F{index}", "Female", [
                repeat_record("X1", "chrX:5000000-5000030", "CAG", genotype)]) for index, genotype in enumerate(
                ["10/11", "12/12"])]
        data_path = self.write_qc_stats_file(sample_jsons)
        self.assertIsNone(m.read_qc_stats(data_path, 3).repeats[0].observed_vs_expected_heterozygosity)
        self.assertIsNotNone(m.read_qc_stats(data_path, 2).repeats[0].observed_vs_expected_heterozygosity)

    def test_expansions_and_contractions_need_enough_genotyped_samples(self):
        data_path = self.write_qc_stats_file(SAMPLE_JSONS)
        # L1 is called in all 4 samples and L2 in 3 (D has no call)
        for min_genotyped_samples, min_fraction_of_samples, expected_counts in [
            (1, 0, (1, 1)),
            (4, 0, (1, 0)),
            (5, 0, (0, 0)),
            # 0.8 of 4 samples rounds up to 4
            (1, 0.8, (1, 0)),
            # the dashboard's own thresholds, 10 samples and half of them
            (10, 0.5, (0, 0)),
        ]:
            with self.subTest(min_genotyped_samples=min_genotyped_samples,
                              min_fraction_of_samples=min_fraction_of_samples), \
                    count_expansions_and_contractions_at_repeats_called_in(min_genotyped_samples,
                                                                          min_fraction_of_samples):
                rows = {row["SampleId"]: row for row in m.read_qc_stats(data_path, 1).sample_rows}
                self.assertEqual((rows["D"]["LociWithLargestExpansion"], rows["C"]["LociWithLargestContraction"]),
                                 expected_counts)

    def test_sex_check(self):
        self.assertEqual([m.infer_sex_karyotype(x, y) for x, y in [(0.5, 0.5), (1, 0), (1, 0.5), (0.5, 0), (0, 0)]],
                         ["XY", "XX", "XXY", "X", "neither X nor Y"])
        self.assertIsNone(m.infer_sex_karyotype(None, 0.5))
        self.assertEqual([m.get_sex_check_problem(sex, karyotype) for sex, karyotype in [
            ("Male", "XY"), ("Female", "XY"), ("male", "XX"), ("Male", "XXY"), ("", "XY"), ("Female", None)]],
            [None, "mismatch", "mismatch", "unusual karyotype", None, None])

    def test_suggested_cutoffs(self):
        values = [10, 11, 12, 13, 14, 15, 16, 17, 18, 100]
        both_sides = m.MetricDefinition("Metric", "", "", "linear", "both", None, [])
        # Q1 = 12.25, Q3 = 16.75, IQR = 4.5
        self.assertEqual(m.compute_suggested_cutoffs(np.array(values, dtype=float), both_sides, 3), (-1.25, 30.25))
        self.assertEqual(m.compute_suggested_cutoffs(np.array(values, dtype=float), both_sides._replace(
            suggested_side="high"), 3), (None, 30.25))
        self.assertEqual(m.compute_suggested_cutoffs(np.array(values, dtype=float), both_sides._replace(
            suggested_side=None), 3), (None, None))
        # a metric whose middle half of the values are all the same gets no Tukey fences
        self.assertEqual(m.compute_suggested_cutoffs(np.array([1, 1, 1, 1, 0.5]), both_sides, 3), (None, None))
        # a fence beyond the end of the log odds scale, where a fraction of 0 is placed, isn't suggested, since in the
        # metric's own units it would exclude the 0 (heterozygosity rates of the test cohort)
        logit_scale = both_sides._replace(scale="logit")
        self.assertEqual(m.compute_suggested_cutoffs(m.transform_to_scale([0, 0.66666667, 0.33333333, 0.5], "logit"),
                                                     logit_scale, 3)[0], None)
        # the repeats' call rate gets the fixed cutoff, whatever its values
        locus_call_rate = next(metric for metric in m.LOCUS_METRICS if metric.name == "CallRate")
        self.assertEqual(m.compute_suggested_cutoffs(np.array([1, 1, 1, 1, 0.5]), locus_call_rate, 3), (0.9, None))
        # cutoffs computed on a log scale are transformed back to the metric's units: the logs of 1, 2, 4 and 8 have
        # Q1 = 0.75 log(2) and Q3 = 2.25 log(2), so with k = 1 the fences are -0.75 log(2) and 3.75 log(2)
        log_scale = both_sides._replace(scale="log")
        lower, upper = m.compute_suggested_cutoffs(m.transform_to_scale([1, 2, 4, 8], "log"), log_scale, 1)
        self.assertAlmostEqual(math.log(lower), -0.75 * math.log(2))
        self.assertAlmostEqual(math.log(upper), 3.75 * math.log(2))

    def test_suggested_cutoffs_are_only_on_sides_that_have_a_cutoff(self):
        sides_by_setting = {"both": {"low", "high"}, "low": {"low"}, "high": {"high"}, None: set()}
        for metric in m.SAMPLE_METRICS + m.LOCUS_METRICS:
            with self.subTest(metric=metric.name):
                self.assertLessEqual(sides_by_setting[metric.suggested_side], sides_by_setting[metric.cutoff_sides])

    def test_every_sample_metric_is_in_one_section_and_its_help_can_be_formatted(self):
        # MostCommonReadLength is shown as the sequencing section's table of read lengths
        metrics_in_sections = [metric for _, _, metrics, _ in m.SAMPLE_SECTIONS for metric in metrics]
        self.assertEqual(sorted(metrics_in_sections + ["MostCommonReadLength"]),
                         sorted(metric.name for metric in m.SAMPLE_METRICS))
        qc_stats = m.read_qc_stats(self.write_qc_stats_file(SAMPLE_JSONS), 1)
        help_values = m.get_help_text_values(qc_stats, 10)
        for metric in m.SAMPLE_METRICS + m.LOCUS_METRICS:
            for item in metric.help:
                item.format(**help_values)

    def test_axis_ticks_label_the_metric_units(self):
        # on a log odds scale, the values kept away from 0 and 1 are labeled 0 and 1
        low, high = m.transform_to_scale([0.5, 1.0], "logit")
        tick_positions, tick_labels = m.make_axis_ticks("logit", low, high)
        self.assertEqual((tick_labels[0], tick_labels[-1]), ("0.5", "1"))
        self.assertEqual(tick_positions, sorted(tick_positions))
        self.assertLessEqual(len(tick_labels), m.MAX_AXIS_TICKS + 1)
        _, tick_labels = m.make_axis_ticks("log1p", *m.transform_to_scale([0, 1000], "log1p"))
        self.assertEqual(tick_labels, ["0", "2", "10", "100", "1,000"])
        self.assertIsNone(m.make_axis_ticks("linear", 0, 1))
        # an axis of counts that are all 0 (padded to -1..1 on the log(1 + value) scale) isn't labeled with negative or
        # fractional counts
        self.assertEqual(m.make_axis_ticks("log1p", -1, 1)[1], ["0", "1"])
        # no candidate lies between 20 and 50, so a coverage axis from 21.5 to 42.9 gets evenly spaced ticks
        low, high = m.transform_to_scale([21.5, 42.9], "log1p")
        tick_positions, tick_labels = m.make_axis_ticks("log1p", low, high)
        self.assertEqual(len(tick_labels), m.MIN_AXIS_TICKS)
        self.assertTrue(all(low < position < high for position in tick_positions))
        # evenly spaced ticks between call rates close to 1 keep the digits that differ from 1
        low, high = m.transform_to_scale([0.99991, 0.99998], "logit")
        _, tick_labels = m.make_axis_ticks("logit", low, high)
        self.assertEqual(len(set(tick_labels)), len(tick_labels))
        self.assertNotIn("1", tick_labels)
        self.assertEqual([m.format_number_for_text(value) for value in (0.99877, 0.5, 1234.5, 0.0012345)],
                         ["0.99877", "0.5", "1,234", "0.001234"])

    def test_violin_outline(self):
        outline_values, density = m.compute_violin_outline(np.array([1, 2, 2, 3, 3, 3, 4, 4, 5], dtype=float))
        self.assertEqual((outline_values[0], outline_values[-1], len(density)), (1, 5, m.VIOLIN_OUTLINE_POINTS))
        self.assertAlmostEqual(density.max(), 1)
        self.assertAlmostEqual(outline_values[density.argmax()], 3, delta=0.1)
        self.assertIsNone(m.compute_violin_outline(np.array([2, 2, 2], dtype=float)))

    def test_page_data(self):
        sample_jsons = SAMPLE_JSONS + [expansion_hunter_json("E", "Male", [
            repeat_record("L1", "chr1:100-130", "CAG", "10/10"), repeat_record("L2", "chr1:500-530", "GCC", "8/8"),
            repeat_record("L3", "chr2:100-120", "AAT", "5/5")])]
        data_path = self.write_qc_stats_file(sample_jsons)
        qc_stats = m.read_qc_stats(data_path, 1)
        page_data = m.build_page_data(qc_stats, data_path, m.get_output_paths("prefix"), 3, 1)
        json.dumps(page_data, allow_nan=False)
        self.assertEqual(page_data["samples"]["ids"], ["A", "B", "C", "D", "E"])
        self.assertEqual([metric["name"] for metric in page_data["sampleMetrics"]],
                         [metric.name for metric in m.SAMPLE_METRICS])
        call_rate = next(metric for metric in page_data["sampleMetrics"] if metric["name"] == "CallRate")
        self.assertEqual(call_rate["values"], [1, 1, 1, 0.666667, 1])
        self.assertEqual(call_rate["scale"], "logit")
        # the call rates' middle half are all 1, so no Tukey fence is suggested
        self.assertEqual(call_rate["suggestedCutoffs"], [None, None])
        # SAMPLE_JSONS has no coverage, so no sample has a value
        coverage = next(metric for metric in page_data["sampleMetrics"] if metric["name"] == "MedianLocusCoverage")
        self.assertEqual((coverage["values"], coverage["axis"]), ([None] * 5, None))
        locus_call_rate = page_data["locusMetrics"][0]
        # the repeats' call rate gets the fixed suggestion, although its interquartile range is 0 here too
        self.assertEqual(locus_call_rate["suggestedCutoffs"], [0.9, None])
        self.assertEqual(len(locus_call_rate["box"]), 5)
        # the VariantId is "" where it is the LocusId
        self.assertEqual(page_data["loci"]["variantIds"], ["", "", ""])

    def test_clean_state(self):
        state = {"sampleCutoffs": {"CallRate": [0.95, None]}, "manualSampleStatus": {"A": "excluded", "Z": "kept"},
                 "locusReasonEdits": {"L1\tL1": "noisy", "L9\tL9": "gone"}, "selectedSamples": ["B", "Z"],
                 "excludeSexMismatches": True, "unknownKey": 1}
        cleaned_state, dropped_id_counts = m.clean_state(state, {"A", "B"}, {m.get_repeat_key("L1", "L1")})
        self.assertEqual(cleaned_state, {"sampleCutoffs": {"CallRate": [0.95, None]},
                                         "excludeSexMismatches": True, "manualSampleStatus": {"A": "excluded"},
                                         "locusReasonEdits": {"L1\tL1": "noisy"}, "selectedSamples": ["B"]})
        self.assertEqual(dropped_id_counts, {"manualSampleStatus": 1, "locusReasonEdits": 1, "selectedSamples": 1})
        for bad_state, error_text in [([], "JSON object"), ({"selectedLoci": "L1"}, "must be a list"),
                                      ({"manualLocusStatus": {"L1": "maybe"}}, "excluded"),
                                      ({"locusCutoffs": {"CallRate": [0.9]}}, "lower threshold"),
                                      ({"sampleReasonEdits": {"A": 3}}, "text")]:
            with self.subTest(bad_state=bad_state), self.assertRaisesRegex(ValueError, error_text):
                m.clean_state(bad_state, {"A"}, {"L1"})

    def test_loci_with_more_than_one_repeat_are_left_out(self):
        # as in catalogs that define HTT both on its own and with its adjacent CCG repeat. Only B has a CCG allele of
        # 30, which would be the largest expansion at HTT_CCG.
        sample_jsons = []
        for sample_id, ccg_genotype in [("A", "7/10"), ("B", "7/30")]:
            sample = expansion_hunter_json(sample_id, "Female", [])
            sample["LocusResults"] = {
                "HTT": {"LocusId": "HTT", "Variants": {
                    "HTT": repeat_record("HTT", "chr4:3074876-3074933", "CAG", "17/19")}},
                "HTT_WITH_ADJACENT_REPEATS": {"LocusId": "HTT_WITH_ADJACENT_REPEATS", "Variants": {
                    "HTT": repeat_record("HTT", "chr4:3074876-3074933", "CAG", "17/19"),
                    "HTT_CCG": repeat_record("HTT_CCG", "chr4:3074939-3074966", "CCG", ccg_genotype)}},
            }
            sample_jsons.append(sample)
        data_path = self.write_qc_stats_file(sample_jsons)
        with count_expansions_and_contractions_at_repeats_called_in(1):
            qc_stats = m.read_qc_stats(data_path, 1)
        self.assertEqual([(repeat.locus_id, repeat.variant_id) for repeat in qc_stats.repeats], [("HTT", "HTT")])
        self.assertEqual(qc_stats.repeat_counts_of_left_out_compound_loci, {"HTT_WITH_ADJACENT_REPEATS": 2})
        self.assertEqual([row["LociWithLargestExpansion"] for row in qc_stats.sample_rows], [0, 0])
        self.assertEqual(m.build_page_data(qc_stats, data_path, m.get_output_paths("prefix"), 3, 1)["loci"][
            "leftOutCompoundLocusIds"], ["HTT_WITH_ADJACENT_REPEATS"])
        with mock.patch("plotly.offline.get_plotlyjs", lambda: ""), contextlib.redirect_stdout(io.StringIO()) as output:
            dashboard = m.create_dashboard(data_path, os.path.join(self.temp_dir, "dashboard"), 3, 1)
        self.assertIn("Left out the 2 repeats of the 1 locus with more than one repeat: HTT_WITH_ADJACENT_REPEATS",
                      output.getvalue())
        self.assertEqual(dashboard.repeat_keys, {"HTT\tHTT"})
        with open(dashboard.export("loci", [["HTT\tHTT", "manually excluded"]])) as f:
            self.assertEqual(f.read(), "LocusId\tReasonsForExclusion\nHTT\tmanually excluded\n")

    def test_damaged_state_file_is_not_overwritten(self):
        data_path = self.write_qc_stats_file(SAMPLE_JSONS)
        output_prefix = os.path.join(self.temp_dir, "dashboard")
        with open(output_prefix + ".qc_dashboard_state.json", "w") as f:
            f.write("{not json")
        with mock.patch("plotly.offline.get_plotlyjs", lambda: ""), \
                self.assertRaisesRegex(ValueError, "Move it away to start over"):
            m.create_dashboard(data_path, output_prefix, 3, 1)

    def start_test_server(self, data_path, output_prefix):
        """Starts the dashboard's server on a free port in a thread, and returns (dashboard, port)."""
        with mock.patch("plotly.offline.get_plotlyjs", lambda: "window.Plotly = {};"):
            dashboard = m.create_dashboard(data_path, output_prefix, 3, 1)
        server = m.start_server(dashboard, 0)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        self.addCleanup(server.server_close)
        self.addCleanup(server.shutdown)
        return dashboard, server.server_address[1]

    def request(self, port, method, path, body=None, host=None, headers=None):
        """Sends a request to the test server and returns (status, body bytes, response headers)."""
        connection = http.client.HTTPConnection("127.0.0.1", port)
        all_headers = {"Host": host or f"127.0.0.1:{port}"}
        all_headers.update(headers or {})
        # the server prints each export
        with contextlib.redirect_stdout(io.StringIO()):
            connection.request(method, path, body=None if body is None else json.dumps(body), headers=all_headers)
            response = connection.getresponse()
            result = response.status, response.read(), dict(response.getheaders())
        connection.close()
        return result

    def test_server(self):
        data_path = self.write_qc_stats_file(SAMPLE_JSONS)
        output_prefix = os.path.join(self.temp_dir, "dashboard")
        dashboard, port = self.start_test_server(data_path, output_prefix)
        token_query = f"?token={dashboard.token}"

        # requests without the token, with a wrong one, or from another Host are refused
        for path, host in [("/", None), ("/data?token=wrong", None), ("/data?token=%C3%A9", None),
                           ("/" + token_query, "evil.example.com"),
                           ("/" + token_query, f"evil.example.com:{port}")]:
            with self.subTest(path=path, host=host):
                self.assertEqual(self.request(port, "GET", path, host=host)[0], 403)
        status, page, _ = self.request(port, "GET", "/" + token_query, host=f"localhost:{port}")
        self.assertEqual(status, 200)
        self.assertIn(f'<meta name="dashboard-token" content="{dashboard.token}">', page.decode())
        self.assertEqual(self.request(port, "GET", "/plotly.js" + token_query)[1], b"window.Plotly = {};")
        status, body, headers = self.request(port, "GET", "/data" + token_query, headers={"Accept-Encoding": "gzip"})
        self.assertEqual((status, headers["Content-Encoding"]), (200, "gzip"))
        page_data = json.loads(gzip.decompress(body))
        self.assertEqual(page_data["samples"]["ids"], ["A", "B", "C", "D"])
        self.assertEqual(json.loads(self.request(port, "GET", "/data" + token_query)[1]), page_data)

        # the state is saved, read back, and loaded again by a new dashboard on the same prefix
        self.assertEqual(json.loads(self.request(port, "GET", "/state" + token_query)[1]), {})
        state = {"sampleCutoffs": {"CallRate": [0.9, None]}, "manualLocusStatus": {"L2\tL2": "excluded"},
                 "sampleReasonEdits": {"D": "low call rate"}, "selectedSamples": ["A"], "excludeSexMismatches": True}
        self.assertEqual(self.request(port, "POST", "/state" + token_query, body=state)[0], 200)
        self.assertEqual(json.loads(self.request(port, "GET", "/state" + token_query)[1]), state)
        with open(output_prefix + ".qc_dashboard_state.json") as f:
            saved_state = json.load(f)
        self.assertEqual((saved_state["version"], saved_state["numSamples"], saved_state["numRepeats"],
                          saved_state["dataPath"]), (m.STATE_FILE_VERSION, 4, 3, os.path.abspath(data_path)))
        restarted_dashboard, restarted_port = self.start_test_server(data_path, output_prefix)
        self.assertEqual(json.loads(self.request(
            restarted_port, "GET", f"/state?token={restarted_dashboard.token}")[1]), state)
        self.assertEqual(self.request(port, "POST", "/state" + token_query, body=[])[0], 400)

        # exports
        status, body, _ = self.request(port, "POST", "/export" + token_query, body={
            "kind": "samples", "rows": [["D", "Call rate 0.6667 < 0.9"], ["B", "manually\texcluded\n"]]})
        self.assertEqual((status, json.loads(body)["path"]), (200, output_prefix + ".samples_to_exclude.tsv"))
        with open(output_prefix + ".samples_to_exclude.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tReasonsForExclusion\nD\tCall rate 0.6667 < 0.9\n"
                                       "B\tmanually excluded\n")
        # the export gets the usual permissions, so that others who run the next steps can read it
        umask = os.umask(0)
        os.umask(umask)
        self.assertEqual(os.stat(output_prefix + ".samples_to_exclude.tsv").st_mode & 0o777, 0o666 & ~umask)
        self.assertEqual(self.request(port, "POST", "/export" + token_query, body={
            "kind": "loci", "rows": [[m.get_repeat_key("L2", "L2"), "Call rate 0.75 < 0.9"]]})[0], 200)
        with open(output_prefix + ".loci_to_exclude.tsv") as f:
            self.assertEqual(f.read(), "LocusId\tReasonsForExclusion\nL2\tCall rate 0.75 < 0.9\n")
        for body, error_text in [({"kind": "samples", "rows": [["Z", "x"]]}, "Unknown SampleId"),
                                 ({"kind": "loci", "rows": [["L2", "x"]]}, "Unknown repeat"),
                                 ({"kind": "other", "rows": []}, "Unknown export kind")]:
            with self.subTest(body=body):
                status, response_body, _ = self.request(port, "POST", "/export" + token_query, body=body)
                self.assertEqual(status, 400)
                self.assertIn(error_text, json.loads(response_body)["error"])
        # an export of an empty list writes just the header
        self.request(port, "POST", "/export" + token_query, body={"kind": "samples", "rows": []})
        with open(output_prefix + ".samples_to_exclude.tsv") as f:
            self.assertEqual(f.read(), "SampleId\tReasonsForExclusion\n")


if __name__ == "__main__":
    unittest.main()
