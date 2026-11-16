"""Starts a local web server with an interactive dashboard for choosing which samples and loci to exclude from later
analyses, from a QC stats file written by generate_multisample_qc_stats.py.

The dashboard plots each sample metric and each repeat metric (see generate_multisample_qc_stats.py for how the QC stats
file's columns are computed, and the dashboard's Methods section for the metrics below), with:
    - a half violin and a box of the metric's distribution, and a point per sample or repeat
    - suggested thresholds (faint dotted lines): Tukey's fences, Q1 - k * IQR and Q3 + k * IQR with
      k = --iqr-multiplier, computed on a scale where the metric is closer to symmetric, or a fixed 0.9 for the
      repeats' call rate
    - a lower threshold, an upper threshold or both (dashed red lines) that you drag, on the sides where the metric's
      values point to a problem (MetricDefinition.cutoff_sides). They start at the plot's edges, so at first nothing is
      excluded, and points beyond a suggested threshold are only highlighted. The median fragment length and the
      observed vs. expected heterozygosity are plotted without thresholds.
A sample or repeat is excluded if any plot's threshold excludes it, unless its status was set by hand (right-click a
point or a selection, then "Toggle whether this sample is on the Excluded list", or these samples, this locus, these
loci), which overrides the thresholds. The reasons for each exclusion can be edited before exporting.

Metrics the dashboard computes from the QC stats file:
    LociWithLargestExpansion       repeats where this sample is the only one with the largest allele, and that allele is
                                   larger than the repeat's most common allele. Only repeats called in at least
                                   MIN_GENOTYPED_SAMPLES_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS (10) samples, and in at
                                   least MIN_FRACTION_OF_SAMPLES_GENOTYPED_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS (half) of
                                   all the samples, are counted, since at a repeat called in few samples, one of them is
                                   often the only one with the largest allele by chance.
    LociWithLargestContraction     the same for the smallest allele, which must be smaller than the most common one
    ExpansionToContractionRatio    (LociWithLargestExpansion + 1) / (LociWithLargestContraction + 1)
    InferredSexKaryotype           2 * ChrXCoverageRatio copies of chrX and 2 * ChrYCoverageRatio copies of chrY,
                                   rounded, such as XX or XY
    repeat CallRate                CallRateGenotypedSampleCount / CallRateSampleCount
    repeat ObservedVsExpectedHeterozygosity
                                   the fraction of two-allele calls whose alleles differ, divided by the fraction
                                   expected from the allele size frequencies (1 - the sum of their squares), for the
                                   repeats with at least --min-two-allele-calls-for-observed-vs-expected-heterozygosity
                                   two-allele calls whose expected fraction is at least
                                   MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED (0.1)

The repeat metrics include every sample in the QC stats file, since recomputing them without some samples would take
each sample's call at every repeat. To recompute them without the excluded samples, rerun
generate_multisample_qc_stats.py with --samples-to-exclude <the exported samples TSV>.

Loci with more than one repeat (compound locus definitions, such as HTT's CAG and CCG repeats) are left out of the
dashboard, including its counts of the loci where a sample has the largest expansion or contraction.
generate_multisample_tsv.py --loci-to-exclude leaves out whole loci, so excluding one repeat of such a locus would also
leave out its other repeats.

The server listens on 127.0.0.1 only, and answers only requests that carry the random token in the URL it prints, so
other web pages open in the browser can't read or change the dashboard. Its files, named after --output-prefix (by
default the QC stats file's path without .qc_stats.tsv.gz, or for a gs:// path its file name in the current
directory), are:

    <prefix>.samples_to_exclude.tsv       written by the "Export excluded samples" button: SampleId,
                                          ReasonsForExclusion. generate_multisample_tsv.py and
                                          generate_multisample_qc_stats.py read it with --samples-to-exclude.
    <prefix>.loci_to_exclude.tsv          written by the "Export excluded loci" button: LocusId,
                                          ReasonsForExclusion. generate_multisample_tsv.py reads it with
                                          --loci-to-exclude.
    <prefix>.qc_dashboard_state.json      the thresholds, the statuses set by hand, the edited reasons, the sex check
                                          switch and the selected points, saved after every change and loaded when the
                                          dashboard is started again with the same prefix

Example:
    python3 start_multisample_qc_dashboard.py gregor_genomes.qc_stats.tsv.gz
"""

import argparse
import collections
import errno
import gzip
import hmac
import html
import http.server
import importlib.util
import json
import math
import os
import secrets
import threading
import urllib.parse

import numpy as np

from expansion_hunter_json_streaming_utils import (
    check_output_prefix, format_count_with_noun, is_gcs_path, open_text, set_up_command_line_run)
from generate_multisample_qc_stats import (
    MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH, MAX_STUTTER_REPEAT_UNITS, MIN_SPANNING_READS_FOR_ALLELE_BALANCE,
    QC_STATS_SUFFIX, SAMPLE_COLUMNS, get_most_common_value, iterate_locus_rows, iterate_sample_lines,
    read_qc_stats_header)

# How each plotted metric is described and plotted: its title, the title of its value axis, the scale it is plotted
# and its suggested cutoffs computed on (see transform_to_scale), which side of its distribution gets a suggested
# cutoff ("low", "high", "both" or None for none), a fixed suggested cutoff that replaces Tukey's fences, the bullets
# of its help text, which are formatted with the values of get_help_text_values, and which sides of its plot have a
# cutoff that can be dragged ("low", "high", "both" or None for none): only the sides where its values point to a
# problem. A side with a suggested cutoff must have a cutoff.
MetricDefinition = collections.namedtuple(
    "MetricDefinition",
    ["name", "title", "axis_title", "scale", "suggested_side", "fixed_suggested_cutoff", "help", "cutoff_sides"],
    defaults=["both"])
# Repeats called in fewer than this fraction of the samples are suggested for exclusion. Most repeats have a call rate
# of 1, so their interquartile range is 0 and Tukey's fences would suggest every repeat with a missing call.
SUGGESTED_MIN_REPEAT_CALL_RATE = 0.9
SAMPLE_METRICS = [
    MetricDefinition("MedianLocusCoverage", "Median locus coverage", "Read depth", "log1p", "both", None, [
        "The median of the Coverage that ExpansionHunter reports at the sample's repeats.",
        "Low: allele sizes are less reliable with fewer reads.",
        "High: usually fine, but often a different sequencing batch."]),
    MetricDefinition("MedianFragmentLength", "Median fragment length", "bp", "linear", None, None, [
        "The median of the FragmentLength that ExpansionHunter reports at the sample's repeats: the length of the DNA "
        "fragments that its read pairs come from.",
        "Differences usually reflect library preparation batches, so the plot has no threshold."],
        cutoff_sides=None),
    # shown as a table of the number of samples with each read length, whose cutoff lines snap to the boundaries
    # between its rows
    MetricDefinition("MostCommonReadLength", "Read length", "Read length (bp)", "linear", None, None, [
        "The read length that ExpansionHunter reports at most of each sample's repeats.",
        "Samples with shorter reads have lower call rates, and their alleles at long repeats are sized short.",
        "Drag a dashed red line to a boundary between rows to exclude the samples with longer reads (the upper line) "
        "or shorter reads (the lower line)."]),
    MetricDefinition("CallRate", "Call rate", "Fraction of repeats genotyped", "logit", "low", None, [
        "The fraction of the sample's repeats that ExpansionHunter genotyped.",
        "A female sample's chrY repeats are left out, since she has no chrY.",
        "Low: low coverage, shorter reads or a failed library."],
        cutoff_sides="low"),
    MetricDefinition(
        "FractionOfReadsDisagreeingWithHomozygousCalls", "Reads disagreeing with homozygous calls",
        "Fraction of spanning reads", "logit", "high", None, [
            "At the sample's homozygous autosomal calls, the fraction of spanning reads whose size is more than "
            "{max_stutter_repeat_units} repeat unit from the call (reads 1 unit away are counted as stutter).",
            "DNA from another person adds such reads, so this estimates contamination, like gnomAD's CHARR estimate "
            "at SNVs. It isn't a contamination percentage: compare samples to each other.",
            "Only calls whose alleles are at most {max_allele_length_percent} of the read length are used, since reads "
            "span such alleles about equally often."]),
    MetricDefinition(
        "MeanLongerAlleleReadFractionAtHeterozygousCalls", "Allele balance at heterozygous calls",
        "Longer allele's share of reads", "logit", "both", None, [
            "At heterozygous autosomal calls whose alleles are more than {max_stutter_repeat_units} unit apart and "
            "have at least {min_spanning_reads_for_allele_balance} spanning reads between them, the longer allele's "
            "share of those reads, averaged over the calls.",
            "Usually a little below 0.5, since reads span shorter alleles more often.",
            "Far from the other samples: a mix of two people's DNA, or one allele often missed."]),
    MetricDefinition("MeanDifferenceFromReference", "Mean allele size minus reference size", "Repeat units", "linear",
                     "both", None, [
                         "The mean of (allele size - reference size) over the sample's called alleles.",
                         "The reference size is the reference region length divided by the motif length, rounded "
                         "down, as in ExpansionHunter's VCF REF field.",
                         "A sample whose alleles are too short or too long overall, such as one with shorter reads, "
                         "stands out."]),
    MetricDefinition("HeterozygosityRate", "Heterozygosity rate", "Fraction of two-allele autosomal calls", "logit",
                     "both", None, [
                         "The fraction of the sample's two-allele autosomal calls whose alleles differ.",
                         "Low: often the child of related parents. High: possibly contamination or noisy calls."]),
    MetricDefinition("LociWithLargestExpansion", "Loci where the sample has the largest expansion", "Loci", "log1p",
                     "high", None, [
                         "Repeats where this sample is the only one with the largest allele, and that allele is larger "
                         "than the repeat's most common allele.",
                         "Only repeats called in at least {min_genotyped_samples_to_count_expansions_and_contractions} "
                         "samples count, and ties give no credit.",
                         "High: noisy allele sizes, or many real rare expansions."]),
    MetricDefinition("LociWithLargestContraction", "Loci where the sample has the largest contraction", "Loci", "log1p",
                     "high", None, [
                         "Repeats where this sample is the only one with the smallest allele, and that allele is "
                         "smaller than the repeat's most common allele.",
                         "Only repeats called in at least {min_genotyped_samples_to_count_expansions_and_contractions} "
                         "samples count, and ties give no credit.",
                         "High: noisy allele sizes, or many real rare contractions."]),
    MetricDefinition("ExpansionToContractionRatio", "Expansion to contraction ratio", "Ratio", "log", "both", None, [
        "(loci with the largest expansion + 1) / (loci with the largest contraction + 1)",
        "A sample with noisy allele sizes has many of both and a typical ratio; a ratio far from the other samples' "
        "means that its extreme alleles go mostly one way."]),
]
LOCUS_METRICS = [
    MetricDefinition("CallRate", "Call rate", "Fraction of samples genotyped", "logit", "low",
                     SUGGESTED_MIN_REPEAT_CALL_RATE, [
                         "The fraction of samples with a genotype at the repeat, leaving out female samples at chrY "
                         "repeats.",
                         "Most repeats have a call rate of 1, so the suggested threshold is a fixed "
                         "{min_repeat_call_rate} rather than Tukey's fences."],
                     cutoff_sides="low"),
    MetricDefinition("ObservedVsExpectedHeterozygosity", "Observed vs. Expected (HWE) Heterozygosity",
                     "Observed / expected heterozygosity", "log", None, None, [
                         "The fraction of two-allele calls whose alleles differ, divided by the fraction expected "
                         "under Hardy-Weinberg equilibrium (HWE) from the repeat's allele size frequencies (1 - the "
                         "sum of their squares).",
                         "Only for repeats with at least "
                         "{min_two_allele_calls_for_observed_vs_expected_heterozygosity} two-allele calls whose "
                         "expected heterozygosity is at least {min_expected_heterozygosity}, since the ratio is noise "
                         "at repeats with few calls or little variation.",
                         "High (too many heterozygous calls): reads from a similar repeat elsewhere in the genome "
                         "likely map here.",
                         "Low (too few): one allele is often missed.",
                         "The plot has no threshold, since the ratio also falls at the most polymorphic repeats, and "
                         "it can't exceed 1 / the expected heterozygosity."],
                     cutoff_sides=None),
]
# The dashboard's sections of sample metrics: (anchor, title, metrics plotted as violins, extra card). The extra card is
# the table of read lengths (MostCommonReadLength) or the sex check plot.
SAMPLE_SECTIONS = [
    ("sequencing", "Sequencing and coverage", ["MedianLocusCoverage", "MedianFragmentLength"], "read lengths"),
    ("genotyping", "Genotyping", ["CallRate"], None),
    ("contamination", "Contamination and allele balance",
     ["FractionOfReadsDisagreeingWithHomozygousCalls", "MeanLongerAlleleReadFractionAtHeterozygousCalls"], None),
    ("sex-check", "Sex check", [], "sex check"),
    ("allele-sizes", "Allele sizes and heterozygosity",
     ["MeanDifferenceFromReference", "HeterozygosityRate", "LociWithLargestExpansion", "LociWithLargestContraction",
      "ExpansionToContractionRatio"], None),
]
# LociWithLargestExpansion and LociWithLargestContraction only count repeats called in at least this many samples and
# in at least this fraction of all the samples
MIN_GENOTYPED_SAMPLES_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS = 10
MIN_FRACTION_OF_SAMPLES_GENOTYPED_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS = 0.5
# At a repeat with little variation, the expected heterozygosity is near 0 and its ratio to the observed one is noise
MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED = 0.1
# Fractions are kept this far from 0 and 1 before their log odds are taken, so that 0 and 1 can be plotted. Call rates
# are fractions with large denominators (a genome's is out of over 200,000 repeats), so this is small.
LOGIT_DISTANCE_FROM_0_AND_1 = 1e-6
# Before the log, ratios are raised to at least this, so that an ObservedVsExpectedHeterozygosity of 0 (a repeat with no
# heterozygous calls) is plotted as 0.01 rather than minus infinity
MIN_RATIO_ON_LOG_SCALE = 0.01
# The violin is a kernel density estimate at this many values, from a histogram of this many bins
VIOLIN_OUTLINE_POINTS = 200
VIOLIN_HISTOGRAM_BINS = 512
# Values that can label an axis on each scale, in order of preference (the roundest first), and how many of them at
# most are shown
LOGIT_AXIS_TICK_CANDIDATES = [0.5, 0.9, 0.99, 0.999, 0.9999, 0.99999, 0.1, 0.01, 0.001, 0.0001, 0.8, 0.95, 0.98, 0.995,
                              0.998, 0.9995, 0.2, 0.05, 0.02, 0.005, 0.002, 0.0005, 0.0002, 0.3, 0.4, 0.6, 0.7, 0.85,
                              0.92, 0.94, 0.96, 0.97, 0.15, 0.25, 0.45, 0.55]
LOG1P_AXIS_TICK_CANDIDATES = [0, 1, 10, 100, 1000, 10000, 100000, 2, 20, 200, 2000, 20000, 5, 50, 500, 5000, 50000]
LOG_AXIS_TICK_CANDIDATES = [0.1, 1, 10, 100, 1000, 0.2, 2, 20, 200, 0.02, 0.05, 0.5, 5, 50, 500]
MAX_AXIS_TICKS = 8
# When fewer of the candidates than this fall on an axis (a narrow distribution, such as coverage between 22x and 43x),
# its ticks are this many evenly spaced positions instead, labeled with their values
MIN_AXIS_TICKS = 3
# Values sent to the page are rounded to this many significant digits, which keeps the data of 200,000 repeats small
SIGNIFICANT_DIGITS_SENT_TO_PAGE = 6
PAGE_TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "multisample_qc_dashboard.html")
# Version 2 identifies a repeat by its LocusId and VariantId (see get_repeat_key), since loci of some catalogs share
# VariantIds
STATE_FILE_VERSION = 2
# The keys of the state the page saves, and the type of each value
STATE_VALUE_TYPES = {
    "sampleCutoffs": dict,
    "locusCutoffs": dict,
    "excludeSexMismatches": bool,
    "manualSampleStatus": dict,
    "manualLocusStatus": dict,
    "sampleReasonEdits": dict,
    "locusReasonEdits": dict,
    "selectedSamples": list,
    "selectedLoci": list,
}
# The state keys that hold sample IDs and repeat keys, which are dropped if the QC stats file no longer has them
STATE_KEYS_OF_SAMPLE_IDS = ("manualSampleStatus", "sampleReasonEdits", "selectedSamples")
STATE_KEYS_OF_REPEAT_KEYS = ("manualLocusStatus", "locusReasonEdits", "selectedLoci")
EXPORT_COLUMNS_BY_KIND = {
    "samples": ["SampleId", "ReasonsForExclusion"],
    "loci": ["LocusId", "ReasonsForExclusion"],
}
# Largest request body the server reads: an export of every repeat of a large catalog with long reasons fits
MAX_REQUEST_BYTES = 512 * 1024 * 1024

QcStats = collections.namedtuple("QcStats", [
    "sample_rows", "repeats", "repeat_counts_of_left_out_compound_loci",
    "min_genotyped_samples_to_count_expansions_and_contractions"])
# A repeat's IDs and metrics, and the samples credited with its largest expansion and contraction (see
# get_samples_with_largest_expansion_and_contraction)
RepeatMetrics = collections.namedtuple("RepeatMetrics", [
    "locus_id", "variant_id", "reference_region", "motif", "genotyped_sample_count", "call_rate",
    "observed_vs_expected_heterozygosity", "sample_with_largest_expansion", "sample_with_largest_contraction"])


def parse_optional_number(value, number_type):
    return number_type(value) if value else None


def parse_sample_line(line):
    """Returns {column: value} for a sample row of a QC stats file, with numbers parsed."""
    row = dict(zip(SAMPLE_COLUMNS, line.rstrip("\n").split("\t")[1:]))
    for column in SAMPLE_COLUMNS[3:]:
        row[column] = parse_optional_number(
            row[column], int if column in ("RepeatRecords", "GenotypedRepeats", "MostCommonReadLength",
                                           "AutosomalDiploidGenotypes") else float)
    return row


def compute_expected_heterozygosity(allele_size_histogram):
    """Returns 1 - the sum of the squared allele size frequencies, or None if there are no alleles."""
    total = sum(allele_size_histogram.values())
    if total == 0:
        return None
    return 1 - sum((count / total) ** 2 for count in allele_size_histogram.values())


def compute_observed_vs_expected_heterozygosity(locus_totals, min_two_allele_calls):
    """Returns a repeat's ObservedVsExpectedHeterozygosity, or None if it has fewer than min_two_allele_calls (at
    least 1) two-allele calls, or its expected heterozygosity is below
    MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED. The two-allele calls are counted, rather than all calls, since
    the observed heterozygosity is a fraction of them: at a chrX repeat in a mostly male cohort, most calls have one
    allele."""
    if locus_totals.diploid_genotyped_sample_count < min_two_allele_calls:
        return None
    expected_heterozygosity = compute_expected_heterozygosity(locus_totals.allele_size_histogram)
    if (expected_heterozygosity is None
            or expected_heterozygosity < MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED):
        return None
    return locus_totals.heterozygous_sample_count / locus_totals.diploid_genotyped_sample_count / expected_heterozygosity


def get_samples_with_largest_expansion_and_contraction(locus_totals, min_genotyped_samples):
    """Returns (the sample credited with a repeat's largest expansion, the one credited with its largest contraction),
    each None if no sample is: the only sample with the largest (smallest) allele, if that allele is larger (smaller)
    than the most common one and the repeat was called in at least min_genotyped_samples samples."""
    if locus_totals.largest_allele is None or locus_totals.genotyped_sample_count < min_genotyped_samples:
        return None, None
    most_common_allele_size = get_most_common_value(locus_totals.allele_size_histogram)
    largest_allele_size, num_samples_with_largest_allele, sample_with_largest_allele = locus_totals.largest_allele
    smallest_allele_size, num_samples_with_smallest_allele, sample_with_smallest_allele = locus_totals.smallest_allele
    return (
        sample_with_largest_allele if num_samples_with_largest_allele == 1 and largest_allele_size >
        most_common_allele_size else None,
        sample_with_smallest_allele if num_samples_with_smallest_allele == 1 and smallest_allele_size <
        most_common_allele_size else None)


def infer_sex_karyotype(chr_x_coverage_ratio, chr_y_coverage_ratio):
    """Returns the sex chromosomes that the coverage ratios suggest, with 2 * ratio copies (rounded) of chrX and of
    chrY, such as "XX", "XY" or "XXY", or None if a ratio is missing."""
    if chr_x_coverage_ratio is None or chr_y_coverage_ratio is None:
        return None
    # math.floor(value + 0.5) rounds halves up, where round() rounds them to the nearest even number
    karyotype = "X" * math.floor(2 * chr_x_coverage_ratio + 0.5) + "Y" * math.floor(2 * chr_y_coverage_ratio + 0.5)
    return karyotype or "neither X nor Y"


def get_sex_check_problem(sex, inferred_sex_karyotype):
    """Returns "mismatch" if the inferred karyotype is XX or XY but doesn't match the sample's Sex (XX for female, XY
    for male), "unusual karyotype" if it is neither XX nor XY, or None if the sample passes the sex check (as one
    without an inferred karyotype does)."""
    if inferred_sex_karyotype is None:
        return None
    if inferred_sex_karyotype not in ("XX", "XY"):
        return "unusual karyotype"
    expected_karyotype = {"female": "XX", "male": "XY"}.get((sex or "").lower())
    if expected_karyotype is not None and inferred_sex_karyotype != expected_karyotype:
        return "mismatch"
    return None


def read_qc_stats(path, min_two_allele_calls_for_observed_vs_expected_heterozygosity):
    """Reads a QC stats file and computes the metrics the dashboard plots, leaving out the loci with more than one
    repeat (compound locus definitions), since generate_multisample_tsv.py --loci-to-exclude leaves out whole loci.

    Args:
        path (str): the QC stats file.
        min_two_allele_calls_for_observed_vs_expected_heterozygosity (int): see
            compute_observed_vs_expected_heterozygosity.

    Returns:
        QcStats: the sample rows ({column: value}, with LociWithLargestExpansion, LociWithLargestContraction,
            ExpansionToContractionRatio, InferredSexKaryotype and SexCheckProblem added), the RepeatMetrics of the
            repeats of loci with one repeat, in file order, and {LocusId: number of repeats} for the loci left out,
            sorted by LocusId.
    """
    repeats = []
    with open_text(path) as f:
        sample_ids = read_qc_stats_header(f, path)
        min_genotyped_samples_to_count_expansions_and_contractions = max(
            MIN_GENOTYPED_SAMPLES_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS,
            math.ceil(MIN_FRACTION_OF_SAMPLES_GENOTYPED_TO_COUNT_EXPANSIONS_AND_CONTRACTIONS * len(sample_ids)))
        for key, locus_totals in iterate_locus_rows(f, path):
            repeats.append(RepeatMetrics(
                key[3], key[4], locus_totals.reference_region, locus_totals.motif,
                locus_totals.genotyped_sample_count,
                locus_totals.call_rate_genotyped_sample_count / locus_totals.call_rate_sample_count
                if locus_totals.call_rate_sample_count else None,
                compute_observed_vs_expected_heterozygosity(
                    locus_totals, min_two_allele_calls_for_observed_vs_expected_heterozygosity),
                *get_samples_with_largest_expansion_and_contraction(
                    locus_totals, min_genotyped_samples_to_count_expansions_and_contractions)))
        sample_rows = [parse_sample_line(line) for _, line in iterate_sample_lines(f, path, sample_ids)]
    # the repeats of a locus aren't always next to each other in the file, which is in genomic order, since another
    # locus can lie between them
    repeat_count_by_locus_id = collections.Counter(repeat.locus_id for repeat in repeats)
    repeats = [repeat for repeat in repeats if repeat_count_by_locus_id[repeat.locus_id] == 1]
    expansion_count_by_sample_id = collections.Counter(
        repeat.sample_with_largest_expansion for repeat in repeats if repeat.sample_with_largest_expansion)
    contraction_count_by_sample_id = collections.Counter(
        repeat.sample_with_largest_contraction for repeat in repeats if repeat.sample_with_largest_contraction)
    for row in sample_rows:
        expansion_count = expansion_count_by_sample_id[row["SampleId"]]
        contraction_count = contraction_count_by_sample_id[row["SampleId"]]
        row["LociWithLargestExpansion"] = expansion_count
        row["LociWithLargestContraction"] = contraction_count
        row["ExpansionToContractionRatio"] = (expansion_count + 1) / (contraction_count + 1)
        row["InferredSexKaryotype"] = infer_sex_karyotype(row["ChrXCoverageRatio"], row["ChrYCoverageRatio"])
        row["SexCheckProblem"] = get_sex_check_problem(row["Sex"], row["InferredSexKaryotype"])
    return QcStats(
        sample_rows=sample_rows, repeats=repeats,
        repeat_counts_of_left_out_compound_loci={
            locus_id: count for locus_id, count in sorted(repeat_count_by_locus_id.items()) if count > 1},
        min_genotyped_samples_to_count_expansions_and_contractions=(
            min_genotyped_samples_to_count_expansions_and_contractions))


def transform_to_scale(values, scale):
    """Puts values on the scale they are plotted and their suggested cutoffs computed on.

    Args:
        values: array of values.
        scale (str): "logit" for fractions, whose log odds log(x / (1 - x)) spread out the values near 0 and 1 (after
            keeping them LOGIT_DISTANCE_FROM_0_AND_1 away from 0 and 1), "log1p" for counts and coverage, as
            log(1 + value), "log" for ratios (raised to at least MIN_RATIO_ON_LOG_SCALE first), or "linear" for the
            values themselves.

    Returns:
        numpy.ndarray: the transformed values.
    """
    values = np.asarray(values, dtype=float)
    if scale == "logit":
        values = np.clip(values, LOGIT_DISTANCE_FROM_0_AND_1, 1 - LOGIT_DISTANCE_FROM_0_AND_1)
        return np.log(values / (1 - values))
    if scale == "log1p":
        return np.log1p(values)
    if scale == "log":
        return np.log(np.maximum(values, MIN_RATIO_ON_LOG_SCALE))
    if scale == "linear":
        return values
    raise ValueError(f"Unknown scale: {scale}")


def transform_back_from_scale(value, scale):
    """Returns a value on a scale of transform_to_scale in the metric's own units."""
    if scale == "logit":
        return float(1 / (1 + np.exp(-value)))
    if scale == "log1p":
        return float(np.expm1(value))
    if scale == "log":
        return float(np.exp(value))
    if scale == "linear":
        return float(value)
    raise ValueError(f"Unknown scale: {scale}")


def get_ends_of_scale(scale):
    """Returns the lowest and highest values transform_to_scale gives on a scale, where it keeps values from going
    further (minus or plus infinity where it doesn't)."""
    if scale == "logit":
        lowest, highest = transform_to_scale([0, 1], "logit")
        return float(lowest), float(highest)
    if scale == "log":
        return float(transform_to_scale([0], "log")[0]), math.inf
    return -math.inf, math.inf


def compute_suggested_cutoffs(transformed_values, metric, iqr_multiplier):
    """Returns a metric's suggested (lower cutoff, upper cutoff) in its own units, each None if not suggested: Tukey's
    fences Q1 - k * IQR and Q3 + k * IQR of the transformed values (at least one), or the metric's fixed cutoff. A
    metric whose interquartile range is 0 (the middle half of its values are all the same) gets no Tukey fences, since
    they would suggest excluding every other value."""
    if metric.suggested_side is None:
        return None, None
    if metric.fixed_suggested_cutoff is not None:
        lower_cutoff = upper_cutoff = metric.fixed_suggested_cutoff
    else:
        first_quartile, third_quartile = np.percentile(transformed_values, [25, 75])
        interquartile_range = third_quartile - first_quartile
        if interquartile_range == 0:
            return None, None
        lower_fence = first_quartile - iqr_multiplier * interquartile_range
        upper_fence = third_quartile + iqr_multiplier * interquartile_range
        # transform_to_scale puts fractions of 0 and 1, and ratios at or below MIN_RATIO_ON_LOG_SCALE, at the ends of
        # their scale. A fence at or beyond an end excludes nothing on the scale, but transformed back to the metric's
        # units (and compared with the values themselves) it would exclude the values there, so it isn't suggested.
        lowest_on_scale, highest_on_scale = get_ends_of_scale(metric.scale)
        lower_cutoff = (None if lower_fence <= lowest_on_scale
                        else transform_back_from_scale(lower_fence, metric.scale))
        upper_cutoff = (None if upper_fence >= highest_on_scale
                        else transform_back_from_scale(upper_fence, metric.scale))
    return (lower_cutoff if metric.suggested_side in ("low", "both") else None,
            upper_cutoff if metric.suggested_side in ("high", "both") else None)


def format_number_for_text(value):
    """Formats a number for the dashboard's text the way the page's formatNumber does: with 4 significant digits, as a
    whole number with thousands separators from 1,000 up, and for fractions just below 1 with the digits that differ
    from 1 (such as 0.99877), since 4 significant digits would round neighbouring call rates to the same value."""
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    if 0.99 < value < 1:
        return f"{value:.{min(10, math.ceil(-math.log10(1 - value)) + 3)}f}".rstrip("0")
    return f"{value:.4g}"


def make_axis_ticks(scale, low, high):
    """Returns the ticks of an axis that spans low to high on a scale of transform_to_scale.

    Args:
        scale (str): the scale.
        low (float): the low end of the axis, on that scale.
        high (float): the high end.

    Returns:
        tuple: (tick positions on the scale, tick labels in the metric's units), at least 1 / MAX_AXIS_TICKS of the
            axis apart and picked in the order of the scale's tick candidates (or MIN_AXIS_TICKS evenly spaced ones
            if fewer candidates fall on the axis), or None for the linear scale, whose ticks Plotly picks.
    """
    if scale == "linear":
        return None
    if scale == "logit":
        # the values that were kept away from 0 and 1 are labeled 0 and 1, and come first so that they keep a tick
        labeled_values = [(LOGIT_DISTANCE_FROM_0_AND_1, "0"), (1 - LOGIT_DISTANCE_FROM_0_AND_1, "1")] + [
            (value, f"{value:g}") for value in LOGIT_AXIS_TICK_CANDIDATES]
    elif scale == "log1p":
        labeled_values = [(value, f"{value:,d}") for value in LOG1P_AXIS_TICK_CANDIDATES]
    else:
        labeled_values = [(MIN_RATIO_ON_LOG_SCALE, f"{MIN_RATIO_ON_LOG_SCALE:g}")] + [
            (value, f"{value:g}") for value in LOG_AXIS_TICK_CANDIDATES]
    positions = transform_to_scale([value for value, _ in labeled_values], scale)
    min_distance_between_ticks = (high - low) / MAX_AXIS_TICKS
    ticks = []
    for position, (_, label) in zip(positions, labeled_values):
        if low <= position <= high and all(abs(position - other) >= min_distance_between_ticks for other, _ in ticks):
            ticks.append((float(position), label))
    # An axis of counts that reaches below 0 (as when every count is 0) keeps its candidate ticks, since evenly spaced
    # ones would label negative and fractional counts
    if len(ticks) < MIN_AXIS_TICKS and not (scale == "log1p" and low < 0 and ticks):
        positions = [float(position) for position in np.linspace(low, high, MIN_AXIS_TICKS + 2)[1:-1]]
        return positions, [format_number_for_text(transform_back_from_scale(position, scale))
                           for position in positions]
    ticks.sort()
    return [position for position, _ in ticks], [label for _, label in ticks]


def compute_violin_outline(transformed_values):
    """Returns the outline of a violin of the values: (VIOLIN_OUTLINE_POINTS values from their minimum to their
    maximum, and a Gaussian kernel density estimate at each, scaled so that its maximum is 1), or None if the values
    are all the same. The estimate uses Silverman's bandwidth, and counts the values in VIOLIN_HISTOGRAM_BINS bins
    first, so that it stays fast for hundreds of thousands of repeats.
    """
    low, high = float(transformed_values.min()), float(transformed_values.max())
    if low == high:
        return None
    first_quartile, third_quartile = np.percentile(transformed_values, [25, 75])
    spread = min(transformed_values.std(), (third_quartile - first_quartile) / 1.349) or transformed_values.std()
    bandwidth = 0.9 * spread * len(transformed_values) ** -0.2
    counts, bin_edges = np.histogram(transformed_values, bins=VIOLIN_HISTOGRAM_BINS, range=(low, high))
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    outline_values = np.linspace(low, high, VIOLIN_OUTLINE_POINTS)
    density = np.exp(-0.5 * ((outline_values[:, None] - bin_centers[None, :]) / bandwidth) ** 2) @ counts
    return outline_values, density / density.max()


def round_for_page(value):
    """Returns a number rounded to SIGNIFICANT_DIGITS_SENT_TO_PAGE significant digits (None stays None)."""
    if value is None:
        return None
    return float(f"{value:.{SIGNIFICANT_DIGITS_SENT_TO_PAGE}g}")


def build_metric_payload(metric, values, iqr_multiplier, help_values):
    """Returns what the page needs to plot one metric.

    Args:
        metric (MetricDefinition): the metric.
        values (list): its value for each sample or repeat, None for those that have none.
        iqr_multiplier (float): k in Tukey's fences.
        help_values (dict): the values its help text is formatted with.

    Returns:
        dict: the metric's name, title, axis title, help bullets, scale, the sides of its plot that have a cutoff,
            values, suggested cutoffs (in its own units, each None if not suggested), median, and, if any sample or
            repeat has a value, its violin outline, the percentiles of its box (5th, 25th, 50th, 75th and 95th) and its
            axis range and ticks, all on its scale.
    """
    present_values = [value for value in values if value is not None]
    payload = {
        "name": metric.name, "title": metric.title, "axisTitle": metric.axis_title, "scale": metric.scale,
        "cutoffSides": metric.cutoff_sides, "help": [item.format(**help_values) for item in metric.help],
        "values": [round_for_page(value) for value in values],
        "suggestedCutoffs": [None, None], "median": None, "violin": None, "box": None, "axis": None,
    }
    if not present_values:
        return payload
    transformed_values = transform_to_scale(present_values, metric.scale)
    suggested_cutoffs = compute_suggested_cutoffs(transformed_values, metric, iqr_multiplier)
    payload["suggestedCutoffs"] = [round_for_page(cutoff) for cutoff in suggested_cutoffs]
    payload["median"] = round_for_page(float(np.median(present_values)))
    outline = compute_violin_outline(transformed_values)
    if outline is not None:
        payload["violin"] = {"y": [round_for_page(value) for value in outline[0]],
                             "density": [round(float(value), 4) for value in outline[1]]}
    payload["box"] = [round_for_page(value) for value in np.percentile(transformed_values, [5, 25, 50, 75, 95])]
    transformed_cutoffs = [float(transform_to_scale([cutoff], metric.scale)[0])
                           for cutoff in suggested_cutoffs if cutoff is not None]
    low = min([float(transformed_values.min())] + transformed_cutoffs)
    high = max([float(transformed_values.max())] + transformed_cutoffs)
    padding = (high - low) * 0.06 or 1
    payload["axis"] = {"range": [low - padding, high + padding], "dataRange": [float(transformed_values.min()),
                                                                             float(transformed_values.max())]}
    ticks = make_axis_ticks(metric.scale, low - padding, high + padding)
    if ticks is not None:
        payload["axis"]["tickvals"], payload["axis"]["ticktext"] = ticks
    return payload


def get_help_text_values(qc_stats, min_two_allele_calls_for_observed_vs_expected_heterozygosity):
    """Returns the values that the metrics' help text is formatted with."""
    return {
        "max_stutter_repeat_units": MAX_STUTTER_REPEAT_UNITS,
        "max_allele_length_percent": f"{MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH:.0%}",
        "min_spanning_reads_for_allele_balance": MIN_SPANNING_READS_FOR_ALLELE_BALANCE,
        "min_genotyped_samples_to_count_expansions_and_contractions": (
            f"{qc_stats.min_genotyped_samples_to_count_expansions_and_contractions:,d}"),
        "min_repeat_call_rate": f"{SUGGESTED_MIN_REPEAT_CALL_RATE:g}",
        "min_two_allele_calls_for_observed_vs_expected_heterozygosity": (
            f"{min_two_allele_calls_for_observed_vs_expected_heterozygosity:,d}"),
        "min_expected_heterozygosity": f"{MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED:g}",
    }


def build_page_data(qc_stats, data_path, output_paths, iqr_multiplier,
                    min_two_allele_calls_for_observed_vs_expected_heterozygosity):
    """Returns the data the page loads: the samples, the repeats, and each metric's plot data.

    Args:
        qc_stats (QcStats): from read_qc_stats.
        data_path (str): the QC stats file's path, as the page shows it.
        output_paths (dict): the dashboard's file paths, by kind ("samples", "loci", "state").
        iqr_multiplier (float): k in Tukey's fences.
        min_two_allele_calls_for_observed_vs_expected_heterozygosity (int): see
            compute_observed_vs_expected_heterozygosity.
    """
    help_values = get_help_text_values(qc_stats, min_two_allele_calls_for_observed_vs_expected_heterozygosity)
    sample_rows = qc_stats.sample_rows
    repeats = qc_stats.repeats
    locus_values_by_metric = {
        "CallRate": [repeat.call_rate for repeat in repeats],
        "ObservedVsExpectedHeterozygosity": [repeat.observed_vs_expected_heterozygosity for repeat in repeats]}
    return {
        "title": os.path.basename(data_path),
        "dataPath": data_path,
        "outputPaths": {kind: os.path.abspath(path) for kind, path in output_paths.items()},
        "settings": {
            "iqrMultiplier": iqr_multiplier,
            "logitDistanceFrom0And1": LOGIT_DISTANCE_FROM_0_AND_1,
            "minRatioOnLogScale": MIN_RATIO_ON_LOG_SCALE,
            "minTwoAlleleCallsForObservedVsExpectedHeterozygosity": (
                min_two_allele_calls_for_observed_vs_expected_heterozygosity),
            "minExpectedHeterozygosity": MIN_EXPECTED_HETEROZYGOSITY_TO_COMPARE_WITH_OBSERVED,
            "suggestedMinRepeatCallRate": SUGGESTED_MIN_REPEAT_CALL_RATE,
        },
        "samples": {
            "ids": [row["SampleId"] for row in sample_rows],
            "sex": [row["Sex"] for row in sample_rows],
            "karyotype": [row["InferredSexKaryotype"] for row in sample_rows],
            "sexCheckProblem": [row["SexCheckProblem"] for row in sample_rows],
            "chrXCoverageRatio": [round_for_page(row["ChrXCoverageRatio"]) for row in sample_rows],
            "chrYCoverageRatio": [round_for_page(row["ChrYCoverageRatio"]) for row in sample_rows],
        },
        "loci": {
            "locusIds": [repeat.locus_id for repeat in repeats],
            # "" where the VariantId is the LocusId, as at most loci, which keeps the data small
            "variantIds": ["" if repeat.variant_id == repeat.locus_id else repeat.variant_id for repeat in repeats],
            "regions": [repeat.reference_region for repeat in repeats],
            "motifs": [repeat.motif for repeat in repeats],
            "genotypedSampleCounts": [repeat.genotyped_sample_count for repeat in repeats],
            "leftOutCompoundLocusIds": list(qc_stats.repeat_counts_of_left_out_compound_loci),
        },
        "sampleSections": [{"id": anchor, "title": title, "metrics": metrics, "extraCard": extra_card}
                           for anchor, title, metrics, extra_card in SAMPLE_SECTIONS],
        "sampleMetrics": [build_metric_payload(metric, [row[metric.name] for row in sample_rows], iqr_multiplier,
                                               help_values) for metric in SAMPLE_METRICS],
        "locusMetrics": [build_metric_payload(metric, locus_values_by_metric[metric.name], iqr_multiplier,
                                              help_values) for metric in LOCUS_METRICS],
    }


def is_number_or_none(value):
    return value is None or (isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value))


def check_state_value(key, value):
    """Raises ValueError if a state value of a known key holds something the page doesn't save."""
    if key.endswith("Cutoffs") and not all(
            isinstance(cutoffs, list) and len(cutoffs) == 2 and all(map(is_number_or_none, cutoffs))
            for cutoffs in value.values()):
        raise ValueError(f"The state's {key} must map metrics to [lower threshold, upper threshold]")
    if key.startswith("manual") and not all(status in ("excluded", "kept") for status in value.values()):
        raise ValueError(f"The state's {key} must map IDs to \"excluded\" or \"kept\"")
    if key.endswith("ReasonEdits") and not all(isinstance(reason, str) for reason in value.values()):
        raise ValueError(f"The state's {key} must map IDs to text")
    if key.startswith("selected") and not all(isinstance(item_id, str) for item_id in value):
        raise ValueError(f"The state's {key} must be a list of IDs")


def get_repeat_key(locus_id, variant_id):
    """Returns the ID of a repeat in the saved state and the page's exports: its LocusId and VariantId joined by a tab
    (which neither can contain, since both come from TSV columns). The VariantId alone isn't enough, since loci of some
    catalogs share VariantIds."""
    return f"{locus_id}\t{variant_id}"


def clean_state(state, sample_ids, repeat_keys):
    """Returns the state the page saved with only the known keys of STATE_VALUE_TYPES, and without the sample IDs and
    repeats the QC stats file doesn't have.

    Args:
        state (dict): the state.
        sample_ids (set): the QC stats file's SampleIds.
        repeat_keys (set): the get_repeat_key of each of its repeats.

    Returns:
        tuple: (the cleaned state, {state key: number of IDs dropped})

    Raises:
        ValueError: if the state isn't a JSON object, or a value has the wrong type.
    """
    if not isinstance(state, dict):
        raise ValueError("The state must be a JSON object")
    cleaned_state = {}
    dropped_id_counts = {}
    for key, value_type in STATE_VALUE_TYPES.items():
        if key not in state:
            continue
        value = state[key]
        if not isinstance(value, value_type):
            raise ValueError(f"The state's {key} must be a {value_type.__name__}")
        if value_type is not bool:
            check_state_value(key, value)
        known_ids = (sample_ids if key in STATE_KEYS_OF_SAMPLE_IDS else
                     repeat_keys if key in STATE_KEYS_OF_REPEAT_KEYS else None)
        if known_ids is not None:
            kept = [item for item in value if item in known_ids]
            if len(kept) < len(value):
                dropped_id_counts[key] = len(value) - len(kept)
            value = {item: value[item] for item in kept} if isinstance(value, dict) else kept
        cleaned_state[key] = value
    return cleaned_state, dropped_id_counts


def read_state(state_path, sample_ids, repeat_keys):
    """Returns the state saved at state_path, cleaned by clean_state, or {} if there is none, and prints what was dropped.

    Raises:
        ValueError: if the state file can't be read, so that a damaged file isn't overwritten.
    """
    if not os.path.isfile(state_path):
        return {}
    try:
        with open(state_path) as f:
            saved_state = json.load(f)
        if saved_state.get("version") != STATE_FILE_VERSION:
            raise ValueError(f"version {saved_state.get('version')} rather than {STATE_FILE_VERSION}")
        state, dropped_id_counts = clean_state(saved_state, sample_ids, repeat_keys)
    except (ValueError, AttributeError) as error:
        raise ValueError(f"Couldn't read the saved dashboard state {state_path} ({error}). Move it away to start over.")
    for key, count in dropped_id_counts.items():
        print(f"WARNING: dropped {format_count_with_noun(count, 'ID')} from the saved {key}, since the QC stats file "
              f"doesn't have them")
    return state


def write_file_atomically(path, text):
    """Writes text to a temporary file next to path and renames it to path, so that path is never half-written. The
    temporary file is created with open() rather than tempfile.mkstemp, so that it gets the usual permissions (the umask
    applies) rather than being readable only by its owner."""
    temp_path = f"{path}.{secrets.token_hex(6)}.tmp"
    try:
        with open(temp_path, "x") as f:
            f.write(text)
        os.replace(temp_path, path)
    except BaseException:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise


def remove_tabs_and_newlines(text):
    return " ".join(str(text).replace("\t", " ").splitlines()).strip()


class Dashboard:
    """What the server serves: the page, Plotly's code, the page data, and the saved state, and where it writes the
    exports and the state."""

    def __init__(self, page, plotly_js, page_data, sample_ids, repeat_keys, output_paths, data_path, state):
        self.token = secrets.token_urlsafe(24)
        self.page = page
        self.plotly_js = plotly_js.encode()
        self.page_data_gzipped = gzip.compress(json.dumps(page_data, allow_nan=False, separators=(",", ":")).encode(),
                                               compresslevel=6)
        self.sample_ids = sample_ids
        self.repeat_keys = repeat_keys
        self.output_paths = output_paths
        self.data_path = data_path
        self.state = state
        self.lock = threading.Lock()

    def save_state(self, state):
        """Cleans and saves the state the page sent, with the QC stats file's path and its numbers of samples and
        repeats."""
        cleaned_state, _ = clean_state(state, self.sample_ids, self.repeat_keys)
        with self.lock:
            self.state = cleaned_state
            write_file_atomically(self.output_paths["state"], json.dumps(dict(
                {"version": STATE_FILE_VERSION, "dataPath": self.data_path,
                 "numSamples": len(self.sample_ids), "numRepeats": len(self.repeat_keys)},
                **cleaned_state), indent=1) + "\n")

    def export(self, kind, rows):
        """Writes the exported samples or repeats to their TSV, and returns its path.

        Args:
            kind (str): "samples" or "loci".
            rows (list): [SampleId or get_repeat_key of the repeat, reason] for each sample or repeat to exclude.

        Raises:
            ValueError: if kind is unknown or a row's ID isn't in the QC stats file.
        """
        if kind not in EXPORT_COLUMNS_BY_KIND:
            raise ValueError(f"Unknown export kind: {kind}")
        lines = ["\t".join(EXPORT_COLUMNS_BY_KIND[kind])]
        for row in rows:
            if not isinstance(row, list) or len(row) != 2:
                raise ValueError("Each exported row must be [ID, reason]")
            item_id, reason = row
            if kind == "samples":
                if item_id not in self.sample_ids:
                    raise ValueError(f"Unknown SampleId: {item_id}")
                values = [item_id]
            else:
                if item_id not in self.repeat_keys:
                    raise ValueError(f"Unknown repeat: {item_id}")
                # the LocusId alone, since generate_multisample_tsv.py --loci-to-exclude reads only that, and the
                # dashboard leaves out loci with more than one repeat
                values = [item_id.split("\t")[0]]
            lines.append("\t".join(values + [remove_tabs_and_newlines(reason)]))
        path = self.output_paths[kind]
        with self.lock:
            write_file_atomically(path, "\n".join(lines) + "\n")
        return os.path.abspath(path)


def make_request_handler_class(dashboard, allowed_hosts):
    """Returns the request handler class of a server for the dashboard.

    Args:
        dashboard (Dashboard): what the server serves.
        allowed_hosts (set): the Host header values the server answers (127.0.0.1 and localhost, with its port), so
            that a web page on another domain that resolves to 127.0.0.1 can't reach it.
    """

    class DashboardRequestHandler(http.server.BaseHTTPRequestHandler):
        server_version = "MultisampleQcDashboard"

        def log_message(self, message_format, *args):
            """Leaves the terminal free for the dashboard's own messages."""

        def send_body(self, status, body, content_type, gzipped=False):
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Cache-Control", "no-store")
            if gzipped:
                self.send_header("Content-Encoding", "gzip")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def send_json(self, status, value):
            self.send_body(status, json.dumps(value).encode(), "application/json")

        def check_request(self):
            """Returns the request's path if its Host is allowed and it carries the token, and otherwise answers it
            with 403 and returns None."""
            url = urllib.parse.urlsplit(self.path)
            token = urllib.parse.parse_qs(url.query).get("token", [""])[0]
            # compared as bytes, since compare_digest raises an error for text with non-ASCII characters
            if self.headers.get("Host") not in allowed_hosts or not hmac.compare_digest(token.encode(),
                                                                                          dashboard.token.encode()):
                self.send_body(403, b"Forbidden: open the URL that start_multisample_qc_dashboard.py printed",
                               "text/plain")
                return None
            return url.path

        def do_GET(self):
            path = self.check_request()
            if path is None:
                return
            if path == "/":
                self.send_body(200, dashboard.page.replace("__TOKEN__", dashboard.token).encode(),
                               "text/html; charset=utf-8")
            elif path == "/plotly.js":
                self.send_body(200, dashboard.plotly_js, "application/javascript")
            elif path == "/data":
                if "gzip" in self.headers.get("Accept-Encoding", ""):
                    self.send_body(200, dashboard.page_data_gzipped, "application/json", gzipped=True)
                else:
                    self.send_body(200, gzip.decompress(dashboard.page_data_gzipped), "application/json")
            elif path == "/state":
                self.send_json(200, dashboard.state)
            else:
                self.send_body(404, b"Not found", "text/plain")

        def do_POST(self):
            path = self.check_request()
            if path is None:
                return
            content_length = int(self.headers.get("Content-Length") or 0)
            if content_length > MAX_REQUEST_BYTES:
                self.send_json(413, {"error": "The request is too large"})
                return
            try:
                body = json.loads(self.rfile.read(content_length) or b"null")
                if path == "/state":
                    dashboard.save_state(body)
                    self.send_json(200, {"saved": True})
                elif path == "/export":
                    if not isinstance(body, dict):
                        raise ValueError("The export must be a JSON object")
                    output_path = dashboard.export(body.get("kind"), body.get("rows") or [])
                    print(f"Exported {format_count_with_noun(len(body.get('rows') or []), 'row')} to {output_path}")
                    self.send_json(200, {"path": output_path})
                else:
                    self.send_body(404, b"Not found", "text/plain")
            except ValueError as error:
                self.send_json(400, {"error": str(error)})
            except OSError as error:
                self.send_json(500, {"error": f"Couldn't write the file: {error}"})

    return DashboardRequestHandler


def start_server(dashboard, port):
    """Returns a server for the dashboard on 127.0.0.1, on the given port, or on a free port if that one is taken (or
    is 0)."""
    def create_server(port_to_try):
        server = http.server.ThreadingHTTPServer(("127.0.0.1", port_to_try), http.server.BaseHTTPRequestHandler)
        actual_port = server.server_address[1]
        server.RequestHandlerClass = make_request_handler_class(
            dashboard, {f"127.0.0.1:{actual_port}", f"localhost:{actual_port}"})
        return server

    try:
        return create_server(port)
    except OSError as error:
        if error.errno != errno.EADDRINUSE or port == 0:
            raise
        print(f"Port {port} is in use, so using a free port")
        return create_server(0)


def get_output_paths(output_prefix):
    return {"samples": f"{output_prefix}.samples_to_exclude.tsv", "loci": f"{output_prefix}.loci_to_exclude.tsv",
            "state": f"{output_prefix}.qc_dashboard_state.json"}


def create_dashboard(data_path, output_prefix, iqr_multiplier,
                     min_two_allele_calls_for_observed_vs_expected_heterozygosity):
    """Reads the QC stats file and the saved state, and returns the Dashboard."""
    from plotly.offline import get_plotlyjs

    # the path the page and the saved state show: absolute, unless it is a gs:// path
    data_path = data_path if is_gcs_path(data_path) else os.path.abspath(data_path)
    qc_stats = read_qc_stats(data_path, min_two_allele_calls_for_observed_vs_expected_heterozygosity)
    left_out_repeat_counts = qc_stats.repeat_counts_of_left_out_compound_loci
    if left_out_repeat_counts:
        print(f"Left out the {format_count_with_noun(sum(left_out_repeat_counts.values()), 'repeat')} of the "
              f"{len(left_out_repeat_counts):,d} {'locus' if len(left_out_repeat_counts) == 1 else 'loci'} with more "
              f"than one repeat: {', '.join(list(left_out_repeat_counts)[:10])}"
              f"{', ...' if len(left_out_repeat_counts) > 10 else ''}")
    output_paths = get_output_paths(output_prefix)
    sample_ids = {row["SampleId"] for row in qc_stats.sample_rows}
    repeat_keys = {get_repeat_key(repeat.locus_id, repeat.variant_id) for repeat in qc_stats.repeats}
    state = read_state(output_paths["state"], sample_ids, repeat_keys)
    page_data = build_page_data(qc_stats, data_path, output_paths, iqr_multiplier,
                                min_two_allele_calls_for_observed_vs_expected_heterozygosity)
    with open(PAGE_TEMPLATE_PATH) as f:
        page = f.read().replace("__TITLE__", html.escape(f"QC dashboard: {os.path.basename(data_path)}"))
    return Dashboard(page, get_plotlyjs(), page_data, sample_ids, repeat_keys, output_paths, data_path, state)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output-prefix",
                        help=f"Prefix of the exported TSVs and the saved state (see above). Default: the QC stats "
                             f"file's path without {QC_STATS_SUFFIX}, or for a gs:// path its file name in the current "
                             f"directory")
    parser.add_argument("-p", "--port", type=int, default=8765,
                        help="Port to listen on, on 127.0.0.1 only; a free port is used if it is taken, or with 0 "
                             "(default: %(default)s)")
    parser.add_argument("--iqr-multiplier", type=float, default=3.0,
                        help="k in the suggested thresholds, Q1 - k * IQR and Q3 + k * IQR (default: %(default)s, "
                             "Tukey's 'far out' fences)")
    parser.add_argument("--min-two-allele-calls-for-observed-vs-expected-heterozygosity", type=int, default=10,
                        help="Plot the ObservedVsExpectedHeterozygosity of the repeats with at least this many "
                             "two-allele calls (default: %(default)s)")
    parser.add_argument("data_path", help=f"A QC stats file written by generate_multisample_qc_stats.py "
                                          f"({QC_STATS_SUFFIX}), local or gs://")
    args = parser.parse_args(argv)
    if not args.data_path.endswith(QC_STATS_SUFFIX):
        parser.error(f"{args.data_path} must end with {QC_STATS_SUFFIX}")
    if args.iqr_multiplier <= 0:
        parser.error("--iqr-multiplier must be positive")
    if args.min_two_allele_calls_for_observed_vs_expected_heterozygosity < 1:
        parser.error("--min-two-allele-calls-for-observed-vs-expected-heterozygosity must be at least 1")
    if not 0 <= args.port <= 65535:
        parser.error("--port must be between 0 and 65535")
    # the page loads Plotly's code from the plotly package, so the package is checked for before the data is read
    if importlib.util.find_spec("plotly") is None:
        parser.error("The dashboard needs the plotly package: pip3 install plotly")
    if args.output_prefix is None:
        args.output_prefix = os.path.basename(args.data_path) if is_gcs_path(args.data_path) else args.data_path
        args.output_prefix = args.output_prefix[:-len(QC_STATS_SUFFIX)]
    check_output_prefix(parser, args.output_prefix)
    return args


def main(argv=None):
    # so that the URL reaches a log file as soon as it is printed
    set_up_command_line_run()
    args = parse_args(argv)
    print(f"Reading {args.data_path}")
    dashboard = create_dashboard(args.data_path, args.output_prefix, args.iqr_multiplier,
                                 args.min_two_allele_calls_for_observed_vs_expected_heterozygosity)
    server = start_server(dashboard, args.port)
    print(f"Dashboard of {format_count_with_noun(len(dashboard.sample_ids), 'sample')} and "
          f"{format_count_with_noun(len(dashboard.repeat_keys), 'repeat')}. Open:")
    print(f"    http://127.0.0.1:{server.server_address[1]}/?token={dashboard.token}")
    for description, kind in [("Exports go to", "samples"), ("and", "loci"), ("The state is saved to", "state")]:
        print(f"{description} {os.path.abspath(dashboard.output_paths[kind])}")
    print("Press Ctrl-C to stop.")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("Stopped")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
