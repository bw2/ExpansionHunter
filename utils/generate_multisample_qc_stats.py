"""Combines ExpansionHunter output JSON files from many samples into one QC stats file, which holds what
start_multisample_qc_dashboard.py needs to plot each sample's and each repeat's QC metrics.

Inputs are ExpansionHunter output JSON files (.json, .json.gz, .json.bgz) from any of its modes, and QC stats files
written by earlier runs of this script (.qc_stats.tsv.gz), as local or gs:// paths. Since outputs can be inputs, a
large cohort can be combined in a tree: batches of JSONs into QC stats files, then those into one. The output is
<prefix>.qc_stats.tsv.gz, a gzipped text file laid out as:

    ##ExpansionHunterMultisampleQcStats	version=1
    ##SampleIds	<SampleId>	<SampleId>	...    every sample in the file, sorted
    #Loci	<the locus columns below>
    L	...                                 one row per repeat, in genomic order
    #Samples	<the sample columns below>
    S	...                                 one row per sample, sorted by SampleId

A sample's row depends only on its own JSON, so it is complete only once the JSON has been read to the end, which is
why the sample rows come after the loci. Combining QC stats files merges their locus rows and copies their sample rows.

Locus columns (one row per repeat: a locus with several repeats, such as HTT, gets one row per repeat). Except for the
IDs, they are totals over the samples, which add up when files are combined:
    LocusId, VariantId, ReferenceRegion, Motif
    SampleCount                    samples with a record for the repeat
    GenotypedSampleCount           of those, the ones with a genotype
    CallRateSampleCount, CallRateGenotypedSampleCount
                                   the same two counts without the female samples at chrY repeats, since they have no
                                   chrY to call
    AlleleSizeHistogram            the number of called alleles of each size, such as "10x:180,11x:20"
    LargestAlleleSize, NumSamplesWithLargestAllele, OnlySampleWithLargestAllele
                                   the largest called allele size, the number of samples with an allele of that size,
                                   and the ID of that sample if there is only one
    SmallestAlleleSize, NumSamplesWithSmallestAllele, OnlySampleWithSmallestAllele
                                   the same for the smallest allele
    DiploidGenotypedSampleCount    samples with a two-allele call
    HeterozygousSampleCount        of those, the ones whose two alleles differ

Sample columns:
    SampleId, Sex, Path            from the JSON's SampleParameters, and the JSON's path
    RepeatRecords                  the JSON's repeat records
    GenotypedRepeats               of those, the ones with a genotype
    CallRate                       the fraction of the sample's repeat records that have a genotype, leaving out a female
                                   sample's chrY repeats
    MostCommonReadLength           the ReadLength that ExpansionHunter reports at most of the sample's repeats (the
                                   shorter one if tied)
    MedianLocusCoverage            the median of the locus Coverage, rounded to whole numbers
    MedianFragmentLength           the median of the locus FragmentLength (the length of the DNA fragments the sample's
                                   read pairs come from)
    ChrXCoverageRatio, ChrYCoverageRatio
                                   the median Coverage at the sample's chrX repeats, and at its chrY repeats, outside the
                                   pseudoautosomal regions (PSEUDOAUTOSOMAL_REGIONS_GRCH38), divided by the median at its
                                   autosomal repeats: about 1 and 0 in an XX sample, and 0.5 and 0.5 in an XY sample
    MeanDifferenceFromReference    the mean of (allele size - reference size) over the sample's called alleles, in
                                   repeat units. The reference size is the repeat's reference region length divided by
                                   its motif length, rounded down, as in ExpansionHunter's VCF REF.
    AutosomalDiploidGenotypes      two-allele calls on autosomes
    HeterozygosityRate             the fraction of those whose two alleles differ
    FractionOfReadsDisagreeingWithHomozygousCalls
                                   of the spanning reads (CountsOfSpanningReads) at the sample's homozygous autosomal
                                   calls, the fraction whose size is more than MAX_STUTTER_REPEAT_UNITS (1) repeat unit
                                   from the called allele (reads 1 unit away are counted as stutter). DNA from another
                                   person adds such reads, so this estimates contamination, like gnomAD's CHARR estimate
                                   from reference reads at homozygous alternate SNV calls.
    MeanLongerAlleleReadFractionAtHeterozygousCalls
                                   at the sample's heterozygous autosomal calls whose alleles are more than
                                   MAX_STUTTER_REPEAT_UNITS apart and have at least MIN_SPANNING_READS_FOR_ALLELE_BALANCE
                                   (10) spanning reads between them, the mean of the longer allele's share of those
                                   reads. It is about 0.5 when each call's two alleles are equally represented; a mix
                                   of two people's DNA, or alleles that are often missed, moves it away.
                                   Both read-count metrics use only the calls whose alleles are at most
                                   MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH (half) of the read length, since reads
                                   span such alleles about equally often whatever their size.

Each JSON is streamed, its repeats are put back into genomic order, and the inputs are merged one repeat at a time, the
same way as in generate_multisample_tsv.py (both use expansion_hunter_json_streaming_utils.py), so memory use doesn't
grow with the number of samples. ExpansionHunter's optimized-streaming mode writes the loci nearly in genomic order,
with the chromosomes in the order of the BAM or CRAM header; if the header doesn't list them in the standard order (as
one sorted by name doesn't: chr1, chr10, ..., chr19, chr2, ...), the JSON is sorted from its first chromosome out of
order on, one chromosome at a time, through temporary files. Its seeking and streaming modes sort the loci by LocusId,
so a JSON from those two modes is first sorted, one chromosome at a time, into a temporary single-sample QC stats file.
With --threads N, N worker processes each merge a shard of the inputs into a temporary QC stats file next to the
output, and those are then merged.

Example:
    python3 generate_multisample_qc_stats.py --threads 8 -o gregor_genomes --input-paths-file json_paths.txt
    python3 start_multisample_qc_dashboard.py gregor_genomes.qc_stats.tsv.gz
"""

import argparse
import collections
import gzip
import heapq
import operator
import os
import re
import tempfile

import tqdm

from expansion_hunter_json_streaming_utils import (
    CHR_X_NAMES, CHR_Y_NAMES, JSON_SUFFIXES, OUTPUT_COMPRESSION_LEVEL, TEMPORARY_FILE_COMPRESSION_LEVEL, add_count,
    check_output_prefix, create_temporary_directory_next_to, file_exists, format_count_with_noun,
    format_histogram, get_chromosome, get_file_size, get_genomic_sort_key, get_shard_tsv_paths,
    iterate_checking_that_keys_increase, iterate_in_genomic_order, iterate_repeat_records_of_json, leave_out_samples,
    map_in_threads, map_in_worker_processes, merge_files_in_rounds, merge_records_by_key, open_text,
    output_would_overwrite_an_input, parse_genotype, parse_histogram, read_expansion_hunter_json_header,
    read_ids_to_exclude, read_input_paths, remove_files, set_up_command_line_run, split_inputs_into_shards,
    write_lines_in_genomic_order)

QC_STATS_SUFFIX = ".qc_stats.tsv.gz"
FILE_TYPE_LINE = "##ExpansionHunterMultisampleQcStats\tversion=1"
SAMPLE_IDS_LINE_START = "##SampleIds"
LOCUS_COLUMNS = [
    "LocusId",
    "VariantId",
    "ReferenceRegion",
    "Motif",
    "SampleCount",
    "GenotypedSampleCount",
    "CallRateSampleCount",
    "CallRateGenotypedSampleCount",
    "AlleleSizeHistogram",
    "LargestAlleleSize",
    "NumSamplesWithLargestAllele",
    "OnlySampleWithLargestAllele",
    "SmallestAlleleSize",
    "NumSamplesWithSmallestAllele",
    "OnlySampleWithSmallestAllele",
    "DiploidGenotypedSampleCount",
    "HeterozygousSampleCount",
]
SAMPLE_COLUMNS = [
    "SampleId",
    "Sex",
    "Path",
    "RepeatRecords",
    "GenotypedRepeats",
    "CallRate",
    "MostCommonReadLength",
    "MedianLocusCoverage",
    "MedianFragmentLength",
    "ChrXCoverageRatio",
    "ChrYCoverageRatio",
    "MeanDifferenceFromReference",
    "AutosomalDiploidGenotypes",
    "HeterozygosityRate",
    "FractionOfReadsDisagreeingWithHomozygousCalls",
    "MeanLongerAlleleReadFractionAtHeterozygousCalls",
]
# The rows of each part start with a letter, so that a reader can tell them apart
LOCI_HEADER_LINE = "\t".join(["#Loci"] + LOCUS_COLUMNS)
SAMPLES_HEADER_LINE = "\t".join(["#Samples"] + SAMPLE_COLUMNS)
LOCUS_ROW_START = "L\t"
SAMPLE_ROW_START = "S\t"
# The sex check leaves out the repeats in these pseudoautosomal regions (GRCh38; 0-based start, end), where both sexes
# have two copies
PSEUDOAUTOSOMAL_REGIONS_GRCH38 = {
    "X": [(10000, 2781479), (155701382, 156030895)],
    "Y": [(10000, 2781479), (56887902, 57217415)],
}
# The read-count metrics use the counts in CountsOfSpanningReads, "(size, number of reads), ...". Reads within
# MAX_STUTTER_REPEAT_UNITS of an allele's size count as its stutter, and only calls whose alleles are at most
# MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH of the read length are used. A heterozygous call counts toward the allele
# balance if it has at least MIN_SPANNING_READS_FOR_ALLELE_BALANCE spanning reads at its two allele sizes.
SPANNING_READ_COUNT_REGEX = re.compile(r"\((\d+), (\d+)\)")
MAX_STUTTER_REPEAT_UNITS = 1
MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH = 0.5
MIN_SPANNING_READS_FOR_ALLELE_BALANCE = 10
# ExpansionHunter options that leave repeat records out of its JSON output: the ones without a call, the ones with
# reference calls only, or the ones its quick genotyper can't resolve (written as skipped loci without Variants)
OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS = ("--skip-missing-genotypes", "--skip-hom-ref",
                                         "--genotyping-approach only-quick")

# One input. kind is "json" or "qc_stats". sample_id and sex are a JSON's (None for a QC stats file), and sample_ids are
# all the samples in the input. command_line_was_checked says whether a JSON's RunInfo.CommandLine was checked for
# OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS with its header (in JSONs from older ExpansionHunter versions, RunInfo comes
# after the loci, so it is checked once they have been read). loci_sorted_by_locus_id says whether a JSON was written
# by ExpansionHunter's seeking or streaming mode, which sort the loci by LocusId rather than by position.
InputSpec = collections.namedtuple("InputSpec", ["path", "kind", "sample_id", "sex", "command_line_was_checked",
                                                 "sample_ids", "loci_sorted_by_locus_id"], defaults=[False])


def qc_stats_input_spec(path, sample_ids):
    return InputSpec(path, "qc_stats", None, None, True, tuple(sample_ids))


def check_command_line_keeps_every_repeat_record(path, command_line):
    """Raises ValueError if an ExpansionHunter command line has an option in OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS
    (an option's value may follow it after a space or an "=")."""
    # with a space around each word, so that "--skip-hom-ref" doesn't match a longer word that starts with it
    words = f" {' '.join(command_line.replace('=', ' ').split())} "
    options_used = [option for option in OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS if f" {option} " in words]
    if options_used:
        raise ValueError(f"{path} was written with ExpansionHunter's {' and '.join(options_used)}, which leaves repeat "
                         f"records out, and this script needs every repeat's record: a missing call counts against "
                         f"the call rates, and reference calls are part of the allele size distributions")


class SampleTotals:
    """One sample's totals over its repeat records. coverage_histogram maps the locus Coverage, rounded to a whole
    number, to its number of repeat records, and read_length_histogram and fragment_length_histogram map ReadLength and
    FragmentLength to theirs. The sex check's coverage histograms (autosomal, and chrX and chrY outside the
    pseudoautosomal regions) are like coverage_histogram. The call rate counts leave out a female sample's chrY repeats.
    difference_from_reference_sum adds up (allele size - reference size) over the called alleles.
    longer_allele_read_fraction_sum adds up the longer allele's share of the spanning reads over the heterozygous calls
    counted in heterozygous_calls_with_allele_balance."""

    __slots__ = ("sample_id", "sex", "is_female", "path", "repeat_records", "genotyped_repeats",
                 "repeats_counted_in_call_rate", "genotyped_repeats_counted_in_call_rate", "called_alleles",
                 "coverage_histogram", "read_length_histogram", "fragment_length_histogram",
                 "autosomal_coverage_histogram", "chr_x_coverage_histogram", "chr_y_coverage_histogram",
                 "difference_from_reference_sum", "autosomal_diploid_genotypes", "autosomal_heterozygous_genotypes",
                 "spanning_reads_at_homozygous_calls", "spanning_reads_disagreeing_with_homozygous_calls",
                 "longer_allele_read_fraction_sum", "heterozygous_calls_with_allele_balance")

    def __init__(self, input_spec):
        self.sample_id = input_spec.sample_id
        self.sex = input_spec.sex
        self.is_female = (input_spec.sex or "").lower() == "female"
        self.path = input_spec.path
        self.repeat_records = 0
        self.genotyped_repeats = 0
        self.repeats_counted_in_call_rate = 0
        self.genotyped_repeats_counted_in_call_rate = 0
        self.called_alleles = 0
        self.coverage_histogram = {}
        self.read_length_histogram = {}
        self.fragment_length_histogram = {}
        self.autosomal_coverage_histogram = {}
        self.chr_x_coverage_histogram = {}
        self.chr_y_coverage_histogram = {}
        self.difference_from_reference_sum = 0
        self.autosomal_diploid_genotypes = 0
        self.autosomal_heterozygous_genotypes = 0
        self.spanning_reads_at_homozygous_calls = 0
        self.spanning_reads_disagreeing_with_homozygous_calls = 0
        self.longer_allele_read_fraction_sum = 0.0
        self.heterozygous_calls_with_allele_balance = 0


class LocusTotals:
    """Totals for one repeat over some samples: one sample's record, one locus row, or a merge of those.

    allele_size_histogram maps allele size (repeat units) to its number of called alleles. largest_allele and
    smallest_allele are (allele size, number of samples with an allele of that size, and the ID of that sample if there
    is only one, otherwise ""), or None if no sample has a call. diploid_genotyped_sample_count counts the calls with
    two alleles, and heterozygous_sample_count those whose two alleles differ. call_rate_sample_count and
    call_rate_genotyped_sample_count are sample_count and genotyped_sample_count without the female samples at a chrY
    repeat.
    """

    __slots__ = ("motif", "reference_region", "sample_count", "genotyped_sample_count", "call_rate_sample_count",
                 "call_rate_genotyped_sample_count", "allele_size_histogram", "largest_allele", "smallest_allele",
                 "diploid_genotyped_sample_count", "heterozygous_sample_count")

    def __init__(self, motif, reference_region):
        self.motif = motif
        self.reference_region = reference_region
        self.sample_count = 0
        self.genotyped_sample_count = 0
        self.call_rate_sample_count = 0
        self.call_rate_genotyped_sample_count = 0
        self.allele_size_histogram = {}
        self.largest_allele = None
        self.smallest_allele = None
        self.diploid_genotyped_sample_count = 0
        self.heterozygous_sample_count = 0


def compute_median_of_histogram(histogram):
    """Returns the median of the values counted in a histogram ({value: count}), or None if it is empty. With an even
    total count, the median is the mean of the two middle values."""
    total = sum(histogram.values())
    if total == 0:
        return None
    sorted_items = sorted(histogram.items())

    def get_value_at_position(position):
        """Returns the value at a 0-based position in the sorted list of the counted values."""
        cumulative_count = 0
        for value, count in sorted_items:
            cumulative_count += count
            if cumulative_count > position:
                return value

    return (get_value_at_position((total - 1) // 2) + get_value_at_position(total // 2)) / 2


def get_most_common_value(histogram):
    """Returns the most common value of a histogram ({value: count}), such as an allele size or a read length (the
    smallest one if several are tied), or None if it is empty."""
    if not histogram:
        return None
    return min(histogram, key=lambda value: (-histogram[value], value))


def combine_extreme_alleles(first, second, keep_larger):
    """Combines two (allele size, number of samples with it, ID of the only sample with it or "") of the largest or the
    smallest allele of the same repeat in two groups of samples.

    Args:
        first (tuple): one group's largest or smallest allele, or None if it has no call.
        second (tuple): the other group's, or None.
        keep_larger (bool): whether these are largest alleles (otherwise smallest).

    Returns:
        tuple: the largest or smallest allele of both groups combined, or None if neither has a call.
    """
    if first is None or second is None:
        return first if second is None else second
    if first[0] == second[0]:
        return first[0], first[1] + second[1], ""
    larger, smaller = (first, second) if first[0] > second[0] else (second, first)
    return larger if keep_larger else smaller


def get_sex_check_coverage_histogram(sample_totals, chromosome, reference_start, reference_end):
    """Returns the sample's coverage histogram that a repeat's Coverage counts toward in the sex check: the autosomal,
    chrX or chrY one, or None for a repeat that overlaps a pseudoautosomal region."""
    if chromosome in CHR_X_NAMES:
        chromosome_name, histogram = "X", sample_totals.chr_x_coverage_histogram
    elif chromosome in CHR_Y_NAMES:
        chromosome_name, histogram = "Y", sample_totals.chr_y_coverage_histogram
    else:
        return sample_totals.autosomal_coverage_histogram
    if any(reference_start < region_end and reference_end > region_start
           for region_start, region_end in PSEUDOAUTOSOMAL_REGIONS_GRCH38[chromosome_name]):
        return None
    return histogram


def add_spanning_read_counts(sample_totals, alleles, motif, read_length, counts_of_spanning_reads):
    """Adds an autosomal two-allele call's spanning reads to the sample's FractionOfReadsDisagreeingWithHomozygousCalls
    (if homozygous) or MeanLongerAlleleReadFractionAtHeterozygousCalls (if heterozygous) totals.

    Args:
        sample_totals (SampleTotals): the sample's totals, to update.
        alleles (list): the call's two allele sizes.
        motif (str): the repeat's motif.
        read_length (int): the locus ReadLength, or 0 or None if ExpansionHunter didn't report one.
        counts_of_spanning_reads (str): the call's CountsOfSpanningReads, such as "(10, 12), (14, 8)", or None.
    """
    if not read_length or not counts_of_spanning_reads or (
            max(alleles) * len(motif) > MAX_ALLELE_LENGTH_AS_FRACTION_OF_READ_LENGTH * read_length):
        return
    read_count_by_size = {int(size): int(count)
                          for size, count in SPANNING_READ_COUNT_REGEX.findall(counts_of_spanning_reads)}
    shorter_allele, longer_allele = sorted(alleles)
    if shorter_allele == longer_allele:
        sample_totals.spanning_reads_at_homozygous_calls += sum(read_count_by_size.values())
        sample_totals.spanning_reads_disagreeing_with_homozygous_calls += sum(
            count for size, count in read_count_by_size.items()
            if abs(size - shorter_allele) > MAX_STUTTER_REPEAT_UNITS)
    elif longer_allele - shorter_allele > MAX_STUTTER_REPEAT_UNITS:
        shorter_allele_reads = read_count_by_size.get(shorter_allele, 0)
        longer_allele_reads = read_count_by_size.get(longer_allele, 0)
        if shorter_allele_reads + longer_allele_reads >= MIN_SPANNING_READS_FOR_ALLELE_BALANCE:
            sample_totals.longer_allele_read_fraction_sum += longer_allele_reads / (
                shorter_allele_reads + longer_allele_reads)
            sample_totals.heterozygous_calls_with_allele_balance += 1


def summarize_repeat_record(locus_record, variant_record, reference_start, reference_end, sample_totals):
    """Turns one sample's repeat record from an ExpansionHunter JSON into LocusTotals, and adds it to the sample's
    totals.

    Args:
        locus_record (dict): the record under LocusResults.<locus>, which has the locus Coverage, ReadLength and
            FragmentLength.
        variant_record (dict): the record under LocusResults.<locus>.Variants.
        reference_start (int): the start of the repeat's ReferenceRegion.
        reference_end (int): its end.
        sample_totals (SampleTotals): the sample's totals, to update.

    Returns:
        LocusTotals: the sample's contribution to the repeat's totals.
    """
    locus_totals = LocusTotals(variant_record["RepeatUnit"], variant_record["ReferenceRegion"])
    chromosome = get_chromosome(locus_totals.reference_region)
    is_autosomal = chromosome not in CHR_X_NAMES | CHR_Y_NAMES
    # a female sample has no chrY, so a missing call there doesn't count against either call rate
    counts_in_call_rate = int(not (sample_totals.is_female and chromosome in CHR_Y_NAMES))
    locus_totals.sample_count = 1
    locus_totals.call_rate_sample_count = counts_in_call_rate
    sample_totals.repeat_records += 1
    sample_totals.repeats_counted_in_call_rate += counts_in_call_rate
    if locus_record.get("Coverage") is not None:
        rounded_coverage = round(locus_record["Coverage"])
        add_count(sample_totals.coverage_histogram, rounded_coverage, 1)
        sex_check_coverage_histogram = get_sex_check_coverage_histogram(
            sample_totals, chromosome, reference_start, reference_end)
        if sex_check_coverage_histogram is not None:
            add_count(sex_check_coverage_histogram, rounded_coverage, 1)
    # ExpansionHunter writes a ReadLength and FragmentLength of 0 at loci without reads
    if locus_record.get("ReadLength"):
        add_count(sample_totals.read_length_histogram, locus_record["ReadLength"], 1)
    if locus_record.get("FragmentLength"):
        add_count(sample_totals.fragment_length_histogram, locus_record["FragmentLength"], 1)
    alleles = parse_genotype(variant_record.get("Genotype"))
    if alleles is None:
        return locus_totals

    locus_totals.genotyped_sample_count = 1
    locus_totals.call_rate_genotyped_sample_count = counts_in_call_rate
    sample_totals.genotyped_repeats += 1
    sample_totals.genotyped_repeats_counted_in_call_rate += counts_in_call_rate
    sample_totals.called_alleles += len(alleles)
    if len(alleles) == 2:
        is_heterozygous = alleles[0] != alleles[1]
        locus_totals.diploid_genotyped_sample_count = 1
        locus_totals.heterozygous_sample_count = int(is_heterozygous)
        if is_autosomal:
            sample_totals.autosomal_diploid_genotypes += 1
            sample_totals.autosomal_heterozygous_genotypes += int(is_heterozygous)
            add_spanning_read_counts(sample_totals, alleles, locus_totals.motif, locus_record.get("ReadLength"),
                                     variant_record.get("CountsOfSpanningReads"))
    for allele_size in alleles:
        add_count(locus_totals.allele_size_histogram, allele_size, 1)
    locus_totals.largest_allele = (max(alleles), 1, sample_totals.sample_id)
    locus_totals.smallest_allele = (min(alleles), 1, sample_totals.sample_id)
    # the reference size the way ExpansionHunter's VCF REF computes it
    reference_size = (reference_end - reference_start) // len(locus_totals.motif)
    sample_totals.difference_from_reference_sum += sum(alleles) - len(alleles) * reference_size
    return locus_totals


def merge_locus_totals(target, source, description):
    """Adds source's totals to target's.

    Raises:
        ValueError: if the two have different motifs or reference regions.
    """
    if (target.motif, target.reference_region) != (source.motif, source.reference_region):
        raise ValueError(f"{description} has motif {target.motif} at {target.reference_region} in some inputs and "
                         f"motif {source.motif} at {source.reference_region} in others")
    target.sample_count += source.sample_count
    target.genotyped_sample_count += source.genotyped_sample_count
    target.call_rate_sample_count += source.call_rate_sample_count
    target.call_rate_genotyped_sample_count += source.call_rate_genotyped_sample_count
    target.diploid_genotyped_sample_count += source.diploid_genotyped_sample_count
    target.heterozygous_sample_count += source.heterozygous_sample_count
    for value, count in source.allele_size_histogram.items():
        add_count(target.allele_size_histogram, value, count)
    target.largest_allele = combine_extreme_alleles(target.largest_allele, source.largest_allele, keep_larger=True)
    target.smallest_allele = combine_extreme_alleles(target.smallest_allele, source.smallest_allele, keep_larger=False)


def format_value(value):
    """Formats a value for the QC stats file: "" for None, and a float with 8 significant digits."""
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.8g}"
    return str(value)


def compute_ratio(numerator, denominator):
    """Returns numerator / denominator, or None if the denominator is 0."""
    return numerator / denominator if denominator else None


def compute_coverage_ratio(coverage_histogram, autosomal_median_coverage):
    """Returns the median of a coverage histogram divided by the median autosomal coverage, or None if the histogram is
    empty or the autosomal median is missing or 0 (as in an exome, where most catalog repeats have no reads)."""
    median_coverage = compute_median_of_histogram(coverage_histogram)
    if median_coverage is None or not autosomal_median_coverage:
        return None
    return median_coverage / autosomal_median_coverage


def compute_sample_row(sample_totals):
    """Returns the values of SAMPLE_COLUMNS for a sample whose JSON has been read to the end."""
    autosomal_median_coverage = compute_median_of_histogram(sample_totals.autosomal_coverage_histogram)
    return [
        sample_totals.sample_id,
        sample_totals.sex,
        sample_totals.path,
        sample_totals.repeat_records,
        sample_totals.genotyped_repeats,
        compute_ratio(sample_totals.genotyped_repeats_counted_in_call_rate, sample_totals.repeats_counted_in_call_rate),
        get_most_common_value(sample_totals.read_length_histogram),
        compute_median_of_histogram(sample_totals.coverage_histogram),
        compute_median_of_histogram(sample_totals.fragment_length_histogram),
        compute_coverage_ratio(sample_totals.chr_x_coverage_histogram, autosomal_median_coverage),
        compute_coverage_ratio(sample_totals.chr_y_coverage_histogram, autosomal_median_coverage),
        compute_ratio(sample_totals.difference_from_reference_sum, sample_totals.called_alleles),
        sample_totals.autosomal_diploid_genotypes,
        compute_ratio(sample_totals.autosomal_heterozygous_genotypes, sample_totals.autosomal_diploid_genotypes),
        compute_ratio(sample_totals.spanning_reads_disagreeing_with_homozygous_calls,
                      sample_totals.spanning_reads_at_homozygous_calls),
        compute_ratio(sample_totals.longer_allele_read_fraction_sum,
                      sample_totals.heterozygous_calls_with_allele_balance),
    ]


def format_sample_line(sample_totals):
    """Returns the sample row of a sample whose JSON has been read to the end, with its newline."""
    return SAMPLE_ROW_START + "\t".join(format_value(value) for value in compute_sample_row(sample_totals)) + "\n"


def format_locus_line(key, locus_totals):
    """Returns the locus row of a repeat, with its newline."""
    largest_allele = locus_totals.largest_allele or (None, 0, None)
    smallest_allele = locus_totals.smallest_allele or (None, 0, None)
    values = [
        key[3], key[4], locus_totals.reference_region, locus_totals.motif, locus_totals.sample_count,
        locus_totals.genotyped_sample_count, locus_totals.call_rate_sample_count,
        locus_totals.call_rate_genotyped_sample_count, format_histogram(locus_totals.allele_size_histogram),
        largest_allele[0], largest_allele[1], largest_allele[2] or None,
        smallest_allele[0], smallest_allele[1], smallest_allele[2] or None,
        locus_totals.diploid_genotyped_sample_count, locus_totals.heterozygous_sample_count,
    ]
    return LOCUS_ROW_START + "\t".join(format_value(value) for value in values) + "\n"


def parse_locus_line(line):
    """Returns (sort key, LocusTotals) of a locus row."""
    row = dict(zip(LOCUS_COLUMNS, line[len(LOCUS_ROW_START):].rstrip("\n").split("\t")))
    locus_totals = LocusTotals(row["Motif"], row["ReferenceRegion"])
    locus_totals.sample_count = int(row["SampleCount"])
    locus_totals.genotyped_sample_count = int(row["GenotypedSampleCount"])
    locus_totals.call_rate_sample_count = int(row["CallRateSampleCount"])
    locus_totals.call_rate_genotyped_sample_count = int(row["CallRateGenotypedSampleCount"])
    locus_totals.allele_size_histogram = parse_histogram(row["AlleleSizeHistogram"])
    if row["LargestAlleleSize"]:
        locus_totals.largest_allele = (int(row["LargestAlleleSize"]), int(row["NumSamplesWithLargestAllele"]),
                                       row["OnlySampleWithLargestAllele"])
        locus_totals.smallest_allele = (int(row["SmallestAlleleSize"]), int(row["NumSamplesWithSmallestAllele"]),
                                        row["OnlySampleWithSmallestAllele"])
    locus_totals.diploid_genotyped_sample_count = int(row["DiploidGenotypedSampleCount"])
    locus_totals.heterozygous_sample_count = int(row["HeterozygousSampleCount"])
    return get_genomic_sort_key(row["LocusId"], row["VariantId"], row["ReferenceRegion"]), locus_totals


def write_qc_stats_header(f, sample_ids):
    """Writes the lines that start a QC stats file, through the locus header."""
    f.write(FILE_TYPE_LINE + "\n")
    f.write("\t".join([SAMPLE_IDS_LINE_START] + sorted(sample_ids)) + "\n")
    f.write(LOCI_HEADER_LINE + "\n")


def read_qc_stats_header(f, path):
    """Reads the lines that start a QC stats file, through the locus header, from the open file.

    Returns:
        tuple: the file's SampleIds.

    Raises:
        ValueError: if the file wasn't written by this version of this script.
    """
    file_type_line = f.readline().rstrip("\n")
    if file_type_line != FILE_TYPE_LINE:
        raise ValueError(f"{path} isn't a QC stats file written by this version of generate_multisample_qc_stats.py: "
                         f"its first line is {file_type_line[:100]!r} rather than {FILE_TYPE_LINE!r}")
    sample_ids = f.readline().rstrip("\n").split("\t")
    if sample_ids[0] != SAMPLE_IDS_LINE_START:
        raise ValueError(f"{path}'s second line doesn't start with {SAMPLE_IDS_LINE_START}")
    if f.readline().rstrip("\n") != LOCI_HEADER_LINE:
        raise ValueError(f"{path}'s locus columns differ from the ones this version of the script writes")
    return tuple(sample_ids[1:])


def read_sample_ids_of_qc_stats_file(path):
    """Returns the SampleIds of a QC stats file, read from its header."""
    with open_text(path) as f:
        return read_qc_stats_header(f, path)


def iterate_locus_rows(f, path):
    """Yields (sort key, LocusTotals) for each locus row of a QC stats file open after its header, and stops after the
    sample header, so that the file is left at the first sample row.

    Raises:
        ValueError: if a line among the locus rows is neither a locus row nor the sample header, or the file ends before
            the sample header (as an incomplete file does).
    """
    for line in f:
        if line.startswith(LOCUS_ROW_START):
            yield parse_locus_line(line)
        elif line.rstrip("\n") == SAMPLES_HEADER_LINE:
            return
        else:
            raise ValueError(f"{path} has an unexpected line among its locus rows: {line[:100]!r}")
    raise ValueError(f"{path} ends before its sample rows, so it is incomplete")


def iterate_sample_lines(f, path, sample_ids):
    """Yields (SampleId, line) for each sample row of a QC stats file open at its first sample row.

    Raises:
        ValueError: if the sample rows aren't the file's sample IDs in sorted order.
    """
    expected_sample_ids = iter(sorted(sample_ids))
    for line in f:
        if not line.startswith(SAMPLE_ROW_START):
            raise ValueError(f"{path} has an unexpected line among its sample rows: {line[:100]!r}")
        sample_id = line[len(SAMPLE_ROW_START):].split("\t", 1)[0]
        if sample_id != next(expected_sample_ids, None):
            raise ValueError(f"{path}'s sample rows don't match the SampleIds in its header (at {sample_id})")
        yield sample_id, line
    if next(expected_sample_ids, None) is not None:
        raise ValueError(f"{path} has fewer sample rows than SampleIds in its header, so it is incomplete")


def iterate_locus_totals_of_json(input_spec, sample_totals):
    """Yields (sort key, LocusTotals) for each repeat record of an ExpansionHunter JSON, in file order, and adds the
    records to the sample's totals.

    Raises:
        ValueError: if the JSON's RunInfo, which in JSONs from older ExpansionHunter versions comes after the loci, has
            an option in OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS.
    """
    def check_run_info_after_the_loci(run_info):
        if run_info is not None:
            check_command_line_keeps_every_repeat_record(input_spec.path, run_info.get("CommandLine", ""))

    for locus_id, variant_id, locus_record, variant_record in iterate_repeat_records_of_json(
            input_spec.path, not input_spec.command_line_was_checked, check_run_info_after_the_loci):
        key = get_genomic_sort_key(locus_id, variant_id, variant_record["ReferenceRegion"])
        # key[1] and key[2] are the repeat's start and end
        yield key, summarize_repeat_record(locus_record, variant_record, key[1], key[2], sample_totals)


def get_sort_key_of_locus_line(line):
    """Returns get_genomic_sort_key of a locus row from its LocusId, VariantId and ReferenceRegion columns."""
    locus_id, variant_id, reference_region = line[len(LOCUS_ROW_START):].split("\t", 3)[:3]
    return get_genomic_sort_key(locus_id, variant_id, reference_region)


def write_sorted_qc_stats_file_from_json(input_spec, output_path, temp_dir):
    """Writes a single-sample QC stats file from a seeking or streaming mode JSON, whose loci are sorted by LocusId,
    with its locus rows in genomic order and memory use of about one chromosome's rows (see
    write_lines_in_genomic_order)."""
    sample_totals = SampleTotals(input_spec)
    with gzip.open(output_path, "wt", compresslevel=TEMPORARY_FILE_COMPRESSION_LEVEL) as f:
        write_qc_stats_header(f, [input_spec.sample_id])
        write_lines_in_genomic_order(
            ((key, format_locus_line(key, locus_totals))
             for key, locus_totals in iterate_locus_totals_of_json(input_spec, sample_totals)),
            f, temp_dir, get_sort_key_of_locus_line)
        f.write(SAMPLES_HEADER_LINE + "\n")
        f.write(format_sample_line(sample_totals))


def merge_into_qc_stats_file(input_specs, output_path, compression_level, temp_dir, show_progress_bar=False):
    """Merges inputs into a QC stats file: the locus rows one repeat at a time, then the sample rows.

    Args:
        input_specs (list): InputSpec of each input.
        output_path (str): the QC stats file to write.
        compression_level (int): its gzip compression level.
        temp_dir (str): directory for temporary files.
        show_progress_bar (bool): whether to show a progress bar of the repeats merged.

    Returns:
        int: the number of locus rows written.
    """
    sample_totals_list = []
    qc_stats_files = []
    iterators = []
    sorted_qc_stats_paths = []
    num_locus_rows = 0
    try:
        for input_spec in input_specs:
            if input_spec.kind == "json" and input_spec.loci_sorted_by_locus_id:
                # ExpansionHunter's seeking and streaming modes sort the loci by LocusId, so the JSON is first sorted
                # into a temporary single-sample QC stats file, which is then merged like any other
                file_descriptor, sorted_qc_stats_path = tempfile.mkstemp(suffix=QC_STATS_SUFFIX, dir=temp_dir)
                os.close(file_descriptor)
                sorted_qc_stats_paths.append(sorted_qc_stats_path)
                write_sorted_qc_stats_file_from_json(input_spec, sorted_qc_stats_path, temp_dir)
                input_spec = qc_stats_input_spec(sorted_qc_stats_path, input_spec.sample_ids)
            if input_spec.kind == "json":
                sample_totals = SampleTotals(input_spec)
                sample_totals_list.append(sample_totals)
                iterators.append(iterate_in_genomic_order(
                    iterate_locus_totals_of_json(input_spec, sample_totals), input_spec.path, temp_dir))
            else:
                f = open_text(input_spec.path)
                qc_stats_files.append((input_spec, f))
                read_qc_stats_header(f, input_spec.path)
                iterators.append(iterate_checking_that_keys_increase(
                    iterate_locus_rows(f, input_spec.path), input_spec.path))
        with gzip.open(output_path, "wt", compresslevel=compression_level) as output_file:
            write_qc_stats_header(output_file, [sample_id for input_spec in input_specs
                                               for sample_id in input_spec.sample_ids])
            merged = merge_records_by_key(iterators, merge_locus_totals)
            if show_progress_bar:
                merged = tqdm.tqdm(merged, unit=" repeats", unit_scale=True, desc="Merging")
            for key, locus_totals in merged:
                output_file.write(format_locus_line(key, locus_totals))
                num_locus_rows += 1
            output_file.write(SAMPLES_HEADER_LINE + "\n")
            # Every JSON has now been read to the end, so its sample row is complete. Each input's sample rows are
            # sorted, so merging them keeps the output's sorted.
            sample_line_iterators = [iterate_sample_lines(f, input_spec.path, input_spec.sample_ids)
                                     for input_spec, f in qc_stats_files]
            sample_line_iterators += [[(sample_totals.sample_id, format_sample_line(sample_totals))]
                                      for sample_totals in sample_totals_list]
            for _, line in heapq.merge(*sample_line_iterators, key=operator.itemgetter(0)):
                output_file.write(line)
    finally:
        for iterator in iterators:
            iterator.close()
        for _, f in qc_stats_files:
            f.close()
        remove_files(sorted_qc_stats_paths)
    return num_locus_rows


def merge_input_shard(task):
    """Merges one shard of the inputs into a temporary QC stats file; the worker-process entry point for --threads.

    Args:
        task (tuple): (list of InputSpec, output path, temp dir)
    """
    input_specs, output_path, temp_dir = task
    merge_into_qc_stats_file(input_specs, output_path, TEMPORARY_FILE_COMPRESSION_LEVEL, temp_dir)


def merge_qc_stats_files(paths, output_path, temp_dir, show_progress_bar):
    """Merges QC stats files into one, in rounds if there are more than one merge can read at once, and deletes them.
    With show_progress_bar, the final merge shows a progress bar of the repeats merged.

    Returns:
        int: the number of locus rows written.
    """
    def merge_files(paths_to_merge, merged_path, is_final_merge):
        return merge_into_qc_stats_file(
            [qc_stats_input_spec(path, read_sample_ids_of_qc_stats_file(path)) for path in paths_to_merge], merged_path,
            OUTPUT_COMPRESSION_LEVEL if is_final_merge else TEMPORARY_FILE_COMPRESSION_LEVEL, temp_dir,
            show_progress_bar=show_progress_bar and is_final_merge)

    return merge_files_in_rounds(paths, output_path, temp_dir, merge_files)


def collect_input_specs(input_paths):
    """Reads each JSON's SampleId and Sex and each QC stats file's SampleIds, and each input's file size, in parallel
    threads.

    Returns:
        tuple: (list of InputSpec, {path: file size})

    Raises:
        ValueError: if a JSON was written with an option that leaves repeat records out
            (OPTIONS_THAT_LEAVE_OUT_REPEAT_RECORDS), or a SampleId is in more than one input.
    """
    json_paths = [path for path in input_paths if path.endswith(JSON_SUFFIXES)]
    qc_stats_paths = [path for path in input_paths if not path.endswith(JSON_SUFFIXES)]
    header_by_json_path = dict(zip(json_paths, map_in_threads(read_expansion_hunter_json_header, json_paths)))
    sample_ids_by_qc_stats_path = dict(zip(qc_stats_paths, map_in_threads(read_sample_ids_of_qc_stats_file,
                                                                        qc_stats_paths)))
    file_size_by_path = dict(zip(input_paths, map_in_threads(get_file_size, input_paths)))

    input_specs = []
    path_by_sample_id = {}
    for path in input_paths:
        if path in header_by_json_path:
            header = header_by_json_path[path]
            if header.command_line is not None:
                check_command_line_keeps_every_repeat_record(path, header.command_line)
            input_spec = InputSpec(path, "json", header.sample_id, header.sex, header.command_line is not None,
                                   (header.sample_id,), header.loci_sorted_by_locus_id)
        else:
            input_spec = qc_stats_input_spec(path, sample_ids_by_qc_stats_path[path])
        for sample_id in input_spec.sample_ids:
            if sample_id in path_by_sample_id:
                raise ValueError(f"{path} and {path_by_sample_id[sample_id]} both have SampleId {sample_id}")
            path_by_sample_id[sample_id] = path
        input_specs.append(input_spec)
    return input_specs, file_size_by_path


def check_samples_to_exclude_are_not_in_qc_stats_files(input_specs, sample_ids_to_exclude, samples_to_exclude_path):
    """Raises ValueError if a sample to exclude is in a QC stats file input, whose totals already include it."""
    for input_spec in input_specs:
        if input_spec.kind != "qc_stats":
            continue
        listed_sample_ids = sorted(sample_ids_to_exclude.intersection(input_spec.sample_ids))
        if listed_sample_ids:
            raise ValueError(
                f"{samples_to_exclude_path} lists {format_count_with_noun(len(listed_sample_ids), 'sample')} in the QC "
                f"data file {input_spec.path} (e.g. {listed_sample_ids[0]}), whose totals already include them. Only "
                f"JSON inputs can be left out, so rerun the runs that made that file with --samples-to-exclude.")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output-prefix",
                        help=f"Writes <prefix>{QC_STATS_SUFFIX} (a local path). Default: combined.<N>_samples in the "
                             f"current directory, N being the number of samples")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of worker processes. Each merges a share of the inputs into a temporary QC stats "
                             "file next to the output, and those are then merged (default: %(default)s)")
    parser.add_argument("--input-paths-file",
                        help="A file listing more input paths, one per line, for cohorts too large to list on the "
                             "command line")
    parser.add_argument("--samples-to-exclude",
                        help="A file of samples to leave out: a TSV whose header has a SampleId column, such as the "
                             "samples_to_exclude.tsv that start_multisample_qc_dashboard.py exports, or one SampleId "
                             "per line. Only JSON inputs can be left out.")
    parser.add_argument("--show-progress-bar", action="store_true",
                        help="Show progress bars of the shards and the repeats merged")
    parser.add_argument("input_paths", nargs="*",
                        help=f"ExpansionHunter output JSON files (.json, .json.gz, .json.bgz) from any of its modes, "
                             f"and/or QC stats files written by earlier runs of this script ({QC_STATS_SUFFIX}), local "
                             f"or gs:// paths")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.samples_to_exclude and not file_exists(args.samples_to_exclude):
        parser.error(f"--samples-to-exclude {args.samples_to_exclude} not found")
    if args.output_prefix:
        check_output_prefix(parser, args.output_prefix)
    args.input_paths = read_input_paths(parser, args.input_paths, args.input_paths_file,
                                        JSON_SUFFIXES + (QC_STATS_SUFFIX,))
    return args


def main(argv=None):
    set_up_command_line_run()
    args = parse_args(argv)
    json_count = sum(path.endswith(JSON_SUFFIXES) for path in args.input_paths)
    print(f"Reading the headers of {format_count_with_noun(json_count, 'ExpansionHunter JSON file')} and "
          f"{format_count_with_noun(len(args.input_paths) - json_count, 'QC stats file')}")
    input_specs, file_size_by_path = collect_input_specs(args.input_paths)
    if args.samples_to_exclude:
        sample_ids_to_exclude = read_ids_to_exclude(args.samples_to_exclude, "SampleId")
        check_samples_to_exclude_are_not_in_qc_stats_files(input_specs, sample_ids_to_exclude, args.samples_to_exclude)
        input_specs = leave_out_samples(input_specs, sample_ids_to_exclude, args.samples_to_exclude)
    num_samples = sum(len(input_spec.sample_ids) for input_spec in input_specs)
    output_path = (args.output_prefix or f"combined.{num_samples}_samples") + QC_STATS_SUFFIX
    if output_would_overwrite_an_input(output_path, args.input_paths):
        raise ValueError(f"The output, {output_path}, would overwrite an input; pass a different -o")
    samples_without_sex = sorted(input_spec.sample_id for input_spec in input_specs
                                 if input_spec.kind == "json" and input_spec.sex.lower() not in ("male", "female"))
    if samples_without_sex:
        print(f"WARNING: no Sex in SampleParameters for {format_count_with_noun(len(samples_without_sex), 'sample')} "
              f"(e.g. {samples_without_sex[0]}). Their chrY repeats count toward their call rates, and the dashboard "
              f"can't compare their sex chromosome coverage to their sex.")

    input_specs_by_shard = split_inputs_into_shards(input_specs, args.threads, file_size_by_path)
    num_shards = len(input_specs_by_shard)
    num_workers = min(args.threads, num_shards)
    print(f"Merging {format_count_with_noun(len(input_specs), 'input')} with "
          f"{format_count_with_noun(num_samples, 'sample')}"
          + (f" in {format_count_with_noun(num_shards, 'shard')}" if num_shards > 1 else "")
          + (f" using {num_workers:,d} threads" if num_workers >= 2 else ""))
    with create_temporary_directory_next_to(output_path) as temp_dir:
        merged_path = os.path.join(temp_dir, "merged" + QC_STATS_SUFFIX)
        if num_shards == 1:
            num_locus_rows = merge_into_qc_stats_file(input_specs_by_shard[0], merged_path, OUTPUT_COMPRESSION_LEVEL,
                                                      temp_dir, args.show_progress_bar)
        else:
            shard_paths = get_shard_tsv_paths(temp_dir, num_shards)
            map_in_worker_processes(merge_input_shard, [(input_specs, shard_path, temp_dir) for input_specs, shard_path
                                                        in zip(input_specs_by_shard, shard_paths)],
                                    num_workers, args.show_progress_bar, "shards")
            num_locus_rows = merge_qc_stats_files(shard_paths, merged_path, temp_dir, args.show_progress_bar)
        os.replace(merged_path, output_path)
    print(f"Wrote {format_count_with_noun(num_locus_rows, 'repeat')} and "
          f"{format_count_with_noun(num_samples, 'sample')} to {output_path}")


if __name__ == "__main__":
    main()
