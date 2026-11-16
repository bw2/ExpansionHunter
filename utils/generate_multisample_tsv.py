"""Combine ExpansionHunter output JSON files from many samples into one TSV with one row per locus.

Each input is either an ExpansionHunter output JSON (.json, or .json.gz) or a TSV written by an earlier
run of this script (.tsv, or .tsv.gz), so a large cohort can be combined in batches and the batch tables
combined afterwards. Combining TSVs gives the same output as combining the JSON files they were made from.

All inputs must come from the same catalog: rows are matched by LocusId, VariantId and ReferenceRegion together, so a
repeat whose coordinates differ between two inputs ends up in two rows.

With --threads N, N worker processes each merge a share of the inputs into a temporary TSV, and those are merged into
the output.

To keep memory use low, the inputs are read side by side and merged one locus at a time. That requires that every 
input be pre-sorted into genomic order:

    - ExpansionHunter's optimized-streaming mode already writes the loci nearly in genomic order, with the chromosomes
      in the order of the BAM or CRAM header. If the header doesn't list them in the standard order (as one sorted by
      name doesn't: chr1, chr10, ..., chr19, chr2, ...), the JSON is sorted from its first chromosome out of order
      on, one chromosome at a time, through temporary files.
    - ExpansionHunter's seeking and streaming modes write the loci sorted by LocusId. This script pre-sorts their JSON output 
      into a temporary single-sample TSV, one chromosome at a time, before processing it.


The output has one row per repeat: for loci with a single repeat (every locus in the TRExplorer catalog) that is one
row per LocusId, and a locus with several repeats (for example HTT and HTT_CCG) gets one
row per repeat. Rows are sorted by position. Columns include:

    LocusId, VariantId, ReferenceRegion, Motif   from the ExpansionHunter output
    SampleCount                 samples whose ExpansionHunter output has a record for this locus, regardless of whether 
                                 a genotype was called or not. ExpansionHunter leaves a locus out of a sample's output 
                                 if its --skip-hom-ref or --skip-missing-genotypes options were set, and also skips loci 
                                 that are too wide (more than 2x the sample's read length) or whose motif is too long 
                                 (more than half the read length). A sample without a record for the repeat is not 
                                 counted.
    CalledAlleleCount            called alleles: 2 per diploid genotype, 1 per hemizygous one, 0 per no-call
    AllAlleleHistogram           every called allele, in repeat units: "15x:120,16x:3"
    ShortAlleleHistogram         the shorter allele of each called sample (the only allele of a hemizygous call)
    HemizygousAlleleHistogram    the allele of each hemizygous call: a repeat that ExpansionHunter genotyped with one
                                 allele (the locus's AlleleCount in its JSON output), which is a chrY repeat in any
                                 sample, or a chrX repeat outside the pseudoautosomal regions in a male sample. JSONs
                                 without AlleleCount fall back to the chromosome and SampleParameters.Sex.
    BiallelicHistogram           genotype counts, "short/long:count", sorted: "15/15:58,15/16:3". A hemizygous call's
                                 single allele is written twice ("8" counts as "8/8"). Missing genotypes are not counted.

With --add-outlier-sample-id-columns, the MaxOutlierSampleIds column is set to the value of the --n-outlier-sample-ids argument, 
while the OutlierSampleIds_AllAlleles, OutlierSampleIds_ShortAlleles and OutlierSampleIds_HemizygousAlleles columns are set to a list of the
samples with the largest alleles (eg. "20x:sampleA,18x:sampleB"). Allele sizes are taken from largest to smallest, stopping at the 
first allele size that is shared by --n-outlier-sample-ids or more samples, or after more than that many samples have been listed. 

When the inputs have ExpansionHunter's motif composition (enabled by ExpansionHunter's --output-motif-composition option), these columns are also added:

    MotifComposition    "1:GCG:23,2:GCA:20" (representing motif ID, motif, number of samples in which the motif was observed)
    MotifPairs          "[1][1]:23,[1][2]:18" (representing the IDs of a motif and the motif after it, number of samples in which the motif pair was observed)

Both of these columns are computed from the locus-level Motifs and MotifPairs records of ExpansionHunter's MotifComposition output field (ExpansionHunter's
per-allele Allele1/Allele2 output blocks are not used). ExpansionHunter numbers motifs separately in each sample, so motif 1 can be GCG in one sample and GCA in
another: samples are counted by motif sequence, and both columns renumber motifs by their number of samples, highest first, with ties broken by motif sequence.
If only some inputs have motif composition, a warning is printed and these columns count only the samples in those
inputs. Pass --skip-motif-composition-columns to leave them out.

With --output-motif-occurrence-distributions, a MotifOccurrenceDistributions column is also added. For each motif, it lists how many
samples had each number of occurrences of that motif (summed over the sample's reads), using the motif IDs of the MotifComposition column,
with motifs separated by "|": "1:GCG:10x:5,11x:12,235x:6|2:GCA:1x:18,2x:2" means GCG was seen 10 times in 5 samples, 11 times in 12 samples
and 235 times in 6 samples. A TSV input with motif composition columns must then have this column too.

With -L (--loci), only the repeats of the given LocusIds are kept, for example -L HTT -L ATXN1 or -L HTT,ATXN1.

With --add-consensus-sequence-column, an AlleleSequences column is also added: each distinct consensus allele sequence from
ExpansionHunter's ConsensusSequences field is listed, along with its total haplotype count, most common first, "TGCTGCTGA:86,TGCTGCTGC:56". A
homozygous diploid call's one consensus sequence counts twice, and a hemizygous call's once. Calls without
ConsensusSequences add nothing, so the counts can add up to less than CalledAlleleCount. A warning is printed for
inputs without any consensus sequences.

JSON inputs are checked to make sure that no sample id is present in more than one input file.

Samples and loci can be left out of the output with --samples-to-exclude and --loci-to-exclude, for example the
<prefix>.samples_to_exclude.tsv and <prefix>.loci_to_exclude.tsv that start_multisample_qc_dashboard.py exports. Each
takes a TSV whose header has a SampleId (or LocusId) column, or a list of IDs, one per line. Loci are matched by LocusId
only, so listing a locus with several repeats (such as HTT) leaves out all of its repeats. Samples can only be left out
of JSON inputs, since a TSV from an earlier run has its samples merged already.

Inputs can be local paths or gs:// paths in Google Cloud Storage, which are streamed rather
than downloaded, given on the command line or, for a cohort too large to list there, in a file passed with
--input-paths-file. Reading gs:// paths needs the google-cloud-storage package and application default credentials
(gcloud auth application-default login). The output is written locally.

The code that reads, orders and shards the inputs is in expansion_hunter_json_streaming_utils.py, which
generate_multisample_qc_stats.py uses too.

Example:
    python3 generate_multisample_tsv.py --threads 4 --show-progress-bar --add-outlier-sample-id-columns \\
        -o batch1 sample1.json.gz sample2.json.gz sample3.json.gz
    python3 generate_multisample_tsv.py --add-outlier-sample-id-columns -o cohort \\
        batch1.tsv.gz batch2.tsv.gz sample4.json.gz
"""

import argparse
import collections
import gzip
import itertools
import math
import os
import re
import tempfile

import tqdm

from expansion_hunter_json_streaming_utils import (
    CHR_X_NAMES, CHR_Y_NAMES, JSON_SUFFIXES, MOTIF_COMPOSITION_OPTION, OUTPUT_COMPRESSION_LEVEL,
    TEMPORARY_FILE_COMPRESSION_LEVEL, add_count, check_output_prefix, create_temporary_directory_next_to,
    file_exists, format_count_with_noun, format_histogram, get_chromosome, get_file_size, get_genomic_sort_key,
    get_shard_tsv_paths, is_gcs_path, iterate_checking_that_keys_increase, iterate_in_genomic_order,
    iterate_repeat_records, iterate_repeat_records_of_json, iterate_tsv_rows, leave_out_samples, map_in_threads,
    map_in_worker_processes, merge_files_in_rounds, merge_records_by_key, open_binary, open_text,
    output_would_overwrite_an_input, parse_genotype, parse_histogram, parse_records_from_run_info_on,
    read_expansion_hunter_json_header, read_ids_to_exclude, read_input_paths, read_tail_from_run_info, remove_files,
    run_info_has_motif_composition, set_up_command_line_run, split_inputs_into_shards, strip_suffix,
    write_lines_in_genomic_order)

TSV_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.bgz")

BASE_COLUMNS = [
    "LocusId",
    "VariantId",
    "ReferenceRegion",
    "Motif",
    "SampleCount",
    "CalledAlleleCount",
    "AllAlleleHistogram",
    "ShortAlleleHistogram",
    "HemizygousAlleleHistogram",
    "BiallelicHistogram",
]
MAX_OUTLIER_SAMPLE_IDS_COLUMN = "MaxOutlierSampleIds"
OUTLIER_COLUMNS = [
    "OutlierSampleIds_AllAlleles",
    "OutlierSampleIds_ShortAlleles",
    "OutlierSampleIds_HemizygousAlleles",
]
MOTIF_COMPOSITION_COLUMNS = ["MotifComposition", "MotifPairs"]
MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN = "MotifOccurrenceDistributions"
ALLELE_SEQUENCES_COLUMN = "AlleleSequences"

MOTIF_PAIR_KEY_REGEX = re.compile(r"\[(\d+)\]\[(\d+)\]")
EARLIER_MOTIF_COMPOSITION_FORMAT_EXPLANATION = (
    "TSVs made by earlier versions of this script, which wrote total occurrences and reads (such as "
    "\"1:GCG:561:64\") rather than numbers of samples, can't be combined with this version's. Remake it from the "
    "ExpansionHunter JSON files, or pass --skip-motif-composition-columns.")
# Repeat records of each JSON checked for MotifComposition before the merge, to decide whether to track it.
# ExpansionHunter's all-loci mode puts a MotifComposition record on every eligible repeat, so one shows up within the
# first few records. Its other modes put it on few loci, so for those the check relies on RunInfo.CommandLine.
MOTIF_COMPOSITION_CHECK_RECORDS = 10_000
# Per-input stats keys. The JSON ones are printed as "<count> <key>" after "<count> repeat records".
REPEAT_RECORDS = "repeat records"
WITH_GENOTYPE = "with genotype"
WITH_MOTIF_COMPOSITION = "with motif composition"
WITH_CONSENSUS_SEQUENCES = "with consensus sequences"
WITH_ALLELE_COUNT_NOT_MATCHING_SEX = "with an allele count that doesn't match the sample's sex"
WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE = "whose number of consensus sequences doesn't match the genotype"
TSV_ROWS = "rows"
# These describe the input rather than its records, so they aren't printed. The first is set from
# RunInfo.CommandLine, the second once a JSON has been read to its end (with -n it may not be), and the third counts
# the called repeats that could have ConsensusSequences: ExpansionHunter leaves them out at loci with several repeats.
COMMAND_LINE_HAS_MOTIF_COMPOSITION = "ExpansionHunter command line has --output-motif-composition"
READ_TO_THE_END = "read to the end"
WITH_GENOTYPE_AT_LOCUS_WITH_ONE_REPEAT = "with genotype at a locus with one repeat"

# What every merge needs to know about the run, passed to worker processes. track_motif_occurrence_distributions is
# --output-motif-occurrence-distributions, which only takes effect along with track_motif_composition.
# max_rows_to_write is -n (None for all), locus_ids_to_exclude the LocusIds whose repeats --loci-to-exclude leaves out
# of the output, and locus_ids_to_include the LocusIds given with -L (None to keep every locus).
MergeSettings = collections.namedtuple("MergeSettings", [
    "n_outlier_sample_ids",
    "track_outlier_sample_ids",
    "track_motif_composition",
    "track_motif_occurrence_distributions",
    "track_allele_sequences",
    "max_rows_to_write",
    "locus_ids_to_exclude",
    "locus_ids_to_include",
], defaults=[frozenset(), None])

# One input: kind is "json" or "tsv". sample_id, sex, loci_sorted_by_locus_id (whether the JSON was written by
# ExpansionHunter's seeking or streaming mode) and command_line_has_motif_composition (whether its RunInfo.CommandLine
# has --output-motif-composition, or None if RunInfo couldn't be read before the loci are streamed, as in older
# optimized-streaming output, whose RunInfo is its last record) are set for JSON inputs only. tsv_columns is the header
# of a TSV input, read once so that a gs:// TSV isn't reopened for it (None for JSON inputs and temporary TSVs).
InputSpec = collections.namedtuple(
    "InputSpec", ["path", "kind", "sample_id", "sex", "loci_sorted_by_locus_id", "command_line_has_motif_composition",
                  "tsv_columns"], defaults=[None])


def tsv_input_spec(path, tsv_columns=None):
    return InputSpec(path, "tsv", None, None, False, None, tsv_columns)


class RepeatSummary:
    """Totals for one repeat over some samples: one sample's record, one TSV row, or a merge of those.

    Histograms map allele size (repeat units) to allele count, and the biallelic histogram maps a genotype
    (short allele, long allele) to its number of samples. Outlier sample IDs map allele size to the IDs of up
    to n_outlier_sample_ids samples with that allele, and are None when outlier columns weren't requested. The
    optional counts stay None unless the repeat has data for them: motif sample counts map a motif to the number of
    samples in which it was observed, motif pair sample counts do the same for (motif, next motif), motif occurrence
    distributions map a motif to {number of occurrences in a sample: number of samples}, and allele sequence counts
    map a consensus sequence to its haplotype count.
    """

    __slots__ = (
        "motif", "reference_region", "sample_count", "called_allele_count",
        "all_allele_histogram", "short_allele_histogram", "hemizygous_allele_histogram", "biallelic_histogram",
        "all_allele_outlier_sample_ids", "short_allele_outlier_sample_ids", "hemizygous_allele_outlier_sample_ids",
        "motif_sample_counts", "motif_pair_sample_counts", "motif_occurrence_distributions", "allele_sequence_counts",
    )

    def __init__(self, motif, reference_region, track_outlier_sample_ids):
        self.motif = motif
        self.reference_region = reference_region
        self.sample_count = 0
        self.called_allele_count = 0
        self.all_allele_histogram = {}
        self.short_allele_histogram = {}
        self.hemizygous_allele_histogram = {}
        self.biallelic_histogram = {}
        self.all_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.short_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.hemizygous_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.motif_sample_counts = None
        self.motif_pair_sample_counts = None
        self.motif_occurrence_distributions = None
        self.allele_sequence_counts = None


def get_sort_key_of_row(line):
    """Returns get_genomic_sort_key of an output TSV line from its LocusId, VariantId and ReferenceRegion columns."""
    locus_id, variant_id, reference_region = line.split("\t", 3)[:3]
    return get_genomic_sort_key(locus_id, variant_id, reference_region)


def add_outlier_sample_ids(recorded_sample_ids, allele, sample_ids, n_outlier_sample_ids):
    """Records sample IDs at an allele size, keeping at most n_outlier_sample_ids.

    Once n samples have a size, the outlier walk stops there (see get_outlier_sample_ids_string), so more would
    never be written.
    """
    recorded = recorded_sample_ids.setdefault(allele, [])
    for sample_id in sample_ids:
        if len(recorded) >= n_outlier_sample_ids:
            break
        if sample_id not in recorded:
            recorded.append(sample_id)


def summarize_expansion_hunter_record(variant_record, expected_allele_count, sample_id, is_male, settings, stats):
    """Turns one sample's repeat record from an ExpansionHunter JSON into a RepeatSummary.

    Args:
        variant_record (dict): the record under LocusResults.<locus>.Variants.
        expected_allele_count (int): the locus's AlleleCount, the number of alleles ExpansionHunter genotyped it with
            (1 at chrY, and at chrX outside the pseudoautosomal regions in a male), or None if the JSON doesn't have
            it, in which case the number is worked out from the chromosome and is_male.
        sample_id (str): the sample's ID.
        is_male (bool): whether the sample's SampleParameters.Sex is male.
        settings (MergeSettings): the run's settings.
        stats (collections.Counter): per-input counts to update.

    Returns:
        RepeatSummary: the sample's contribution to the repeat's row.
    """
    summary = RepeatSummary(variant_record["RepeatUnit"], variant_record["ReferenceRegion"],
                            settings.track_outlier_sample_ids)
    summary.sample_count = 1
    stats[REPEAT_RECORDS] += 1

    alleles = parse_genotype(variant_record.get("Genotype"))
    if alleles is not None:
        stats[WITH_GENOTYPE] += 1
        if expected_allele_count is not None:
            expected_hemizygous = expected_allele_count == 1
        else:
            chrom = get_chromosome(summary.reference_region)
            expected_hemizygous = chrom in CHR_Y_NAMES or (is_male and chrom in CHR_X_NAMES)
        if expected_hemizygous != (len(alleles) == 1):
            stats[WITH_ALLELE_COUNT_NOT_MATCHING_SEX] += 1
        summary.called_allele_count = len(alleles)
        for allele in alleles:
            add_count(summary.all_allele_histogram, allele, 1)
        add_count(summary.short_allele_histogram, min(alleles), 1)
        if expected_hemizygous and len(alleles) == 1:
            add_count(summary.hemizygous_allele_histogram, alleles[0], 1)
        add_count(summary.biallelic_histogram, (min(alleles), max(alleles)), 1)
        if settings.track_outlier_sample_ids:
            for histogram, outlier_sample_ids in [
                (summary.all_allele_histogram, summary.all_allele_outlier_sample_ids),
                (summary.short_allele_histogram, summary.short_allele_outlier_sample_ids),
                (summary.hemizygous_allele_histogram, summary.hemizygous_allele_outlier_sample_ids),
            ]:
                for allele in histogram:
                    outlier_sample_ids[allele] = [sample_id]

    motif_composition = variant_record.get("MotifComposition")
    if motif_composition is not None:
        stats[WITH_MOTIF_COMPOSITION] += 1
        if settings.track_motif_composition:
            add_expansion_hunter_motif_composition(summary, motif_composition,
                                                   settings.track_motif_occurrence_distributions)

    consensus_sequences = variant_record.get("ConsensusSequences")
    if consensus_sequences:
        stats[WITH_CONSENSUS_SEQUENCES] += 1
        if settings.track_allele_sequences and alleles is not None \
                and not add_consensus_sequences(summary, alleles, consensus_sequences):
            stats[WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE] += 1
    return summary


def add_expansion_hunter_motif_composition(summary, motif_composition, track_motif_occurrence_distributions):
    """Records which motifs and motif pairs one sample's locus-level MotifComposition record has, translating its
    motif IDs to motif sequences.

    ExpansionHunter lists only the motifs and motif pairs it counted at least once, so each one listed counts as
    observed in the sample.

    Args:
        summary (RepeatSummary): the sample's summary of the repeat.
        motif_composition (dict): the record, e.g. {"Motifs": {"1:GCG": {"count": 561, "reads": 64}},
            "MotifPairs": {"[1][1]": {"count": 398, "reads": 63}}}. An empty record (a locus with no reads) adds
            nothing.
        track_motif_occurrence_distributions (bool): whether to also record each motif's number of occurrences.
    """
    summary.motif_sample_counts, summary.motif_pair_sample_counts = {}, {}
    if track_motif_occurrence_distributions:
        summary.motif_occurrence_distributions = {}
    motif_by_id = {}
    for key, counts in motif_composition.get("Motifs", {}).items():
        motif_id, motif = key.split(":", 1)
        motif_by_id[motif_id] = motif
        summary.motif_sample_counts[motif] = 1
        if track_motif_occurrence_distributions:
            summary.motif_occurrence_distributions[motif] = {counts["count"]: 1}
    for key in motif_composition.get("MotifPairs", {}):
        match = MOTIF_PAIR_KEY_REGEX.fullmatch(key)
        if match is None:
            raise ValueError(f"Unexpected key in an ExpansionHunter MotifComposition.MotifPairs record: {key}")
        summary.motif_pair_sample_counts[(motif_by_id[match.group(1)], motif_by_id[match.group(2)])] = 1


def add_consensus_sequences(summary, alleles, consensus_sequences):
    """Adds one call's ConsensusSequences as haplotype counts.

    ExpansionHunter writes one consensus sequence per distinct allele size, shorter allele first, so a homozygous
    diploid call's one sequence counts for both of its alleles.

    Returns:
        bool: False if the number of sequences doesn't match the number of distinct allele sizes (nothing is
            added then).
    """
    allele_sizes = sorted(set(alleles))
    if len(consensus_sequences) != len(allele_sizes):
        return False
    if summary.allele_sequence_counts is None:
        summary.allele_sequence_counts = {}
    for allele_size, sequence in zip(allele_sizes, consensus_sequences):
        add_count(summary.allele_sequence_counts, sequence, alleles.count(allele_size))
    return True


def merge_repeat_summaries(target, source, n_outlier_sample_ids, description):
    """Adds source's totals to target (two summaries of the same repeat).

    The two have the same sort key, so the same LocusId, VariantId and coordinates; a repeat whose coordinates differ
    between inputs never gets here (see the module docstring).

    Raises:
        ValueError: if the two have a different motif, or their ReferenceRegions name the chromosome differently
            (e.g. chr4 and 4).
    """
    if target.motif != source.motif or target.reference_region != source.reference_region:
        raise ValueError(
            f"{description} is at {source.reference_region} with motif {source.motif} in one input but at "
            f"{target.reference_region} with motif {target.motif} in another. Were the inputs made with different "
            f"catalogs?")
    target.sample_count += source.sample_count
    target.called_allele_count += source.called_allele_count
    for target_histogram, source_histogram, target_outliers, source_outliers in [
        (target.all_allele_histogram, source.all_allele_histogram,
         target.all_allele_outlier_sample_ids, source.all_allele_outlier_sample_ids),
        (target.short_allele_histogram, source.short_allele_histogram,
         target.short_allele_outlier_sample_ids, source.short_allele_outlier_sample_ids),
        (target.hemizygous_allele_histogram, source.hemizygous_allele_histogram,
         target.hemizygous_allele_outlier_sample_ids, source.hemizygous_allele_outlier_sample_ids),
    ]:
        for allele, count in source_histogram.items():
            add_count(target_histogram, allele, count)
        if target_outliers is not None:
            for allele, sample_ids in source_outliers.items():
                add_outlier_sample_ids(target_outliers, allele, sample_ids, n_outlier_sample_ids)
    for genotype, count in source.biallelic_histogram.items():
        add_count(target.biallelic_histogram, genotype, count)
    if source.motif_sample_counts is not None:
        if target.motif_sample_counts is None:
            target.motif_sample_counts, target.motif_pair_sample_counts = {}, {}
        for motif, sample_count in source.motif_sample_counts.items():
            add_count(target.motif_sample_counts, motif, sample_count)
        for pair, sample_count in source.motif_pair_sample_counts.items():
            add_count(target.motif_pair_sample_counts, pair, sample_count)
    if source.motif_occurrence_distributions is not None:
        if target.motif_occurrence_distributions is None:
            target.motif_occurrence_distributions = {}
        for motif, distribution in source.motif_occurrence_distributions.items():
            target_distribution = target.motif_occurrence_distributions.setdefault(motif, {})
            for occurrences, sample_count in distribution.items():
                add_count(target_distribution, occurrences, sample_count)
    if source.allele_sequence_counts is not None:
        if target.allele_sequence_counts is None:
            target.allele_sequence_counts = {}
        for sequence, count in source.allele_sequence_counts.items():
            add_count(target.allele_sequence_counts, sequence, count)


def format_biallelic_histogram(biallelic_histogram):
    return ",".join(f"{short_allele}/{long_allele}:{count}"
                    for (short_allele, long_allele), count in sorted(biallelic_histogram.items()))


def parse_biallelic_histogram(biallelic_histogram_string):
    """Parses "15/15:58,15/16:3" into {(15, 15): 58, (15, 16): 3}."""
    biallelic_histogram = {}
    for entry in filter(None, biallelic_histogram_string.split(",")):
        genotype, count = entry.split(":")
        short_allele, long_allele = genotype.split("/")
        biallelic_histogram[(int(short_allele), int(long_allele))] = int(count)
    return biallelic_histogram


def parse_outlier_sample_ids(outlier_sample_ids_string):
    """Parses "20x:sampleA,20x:sampleB,18x:sampleC" into {20: ["sampleA", "sampleB"], 18: ["sampleC"]}."""
    outlier_sample_ids = {}
    for entry in filter(None, outlier_sample_ids_string.split(",")):
        allele, sample_id = entry.split("x:", 1)
        outlier_sample_ids.setdefault(int(allele), []).append(sample_id)
    return outlier_sample_ids


def get_outlier_sample_ids_string(outlier_sample_ids, histogram, n_outlier_sample_ids):
    """Lists the samples with the largest alleles.

    Walks allele sizes from largest to smallest and stops at the first size that n_outlier_sample_ids or more samples
    have, or once more than n_outlier_sample_ids samples have been listed. When TSVs from earlier runs are combined,
    outlier_sample_ids only holds the samples that each TSV listed, and a TSV lists nothing at or below its own first
    common allele size, so a size that many samples share can arrive with only a few listed samples. It is recognized
    from the histogram, which counts every allele: the alleles of the unlisted samples number at least
    n_outlier_sample_ids. As long as each TSV was made with an n_outlier_sample_ids at least as large, its walk
    reaches at least as far down as the combined walk does, so combining TSVs lists the same samples as combining the
    JSON files they came from.

    Args:
        outlier_sample_ids (dict): allele size -> sample IDs recorded at that size.
        histogram (dict): allele size -> number of alleles of that size across all samples.
        n_outlier_sample_ids (int): number of samples with the largest alleles to list.

    Returns:
        str: comma-separated "{allele size}x:{sample ID}" entries, largest allele size first.
    """
    entries = []
    for allele in sorted(histogram, reverse=True):
        sample_ids = outlier_sample_ids.get(allele, [])
        if len(sample_ids) >= n_outlier_sample_ids:
            break
        # Each listed sample has at most 2 alleles of this size, so when every sample with it is listed,
        # histogram[allele] - len(sample_ids) <= len(sample_ids) < n_outlier_sample_ids.
        if histogram[allele] - len(sample_ids) >= n_outlier_sample_ids:
            break
        entries.extend(f"{allele}x:{sample_id}" for sample_id in sorted(sample_ids))
        if len(entries) > n_outlier_sample_ids:
            break
    return ",".join(entries)


def get_motif_ids(motif_sample_counts):
    """Numbers motifs by their number of samples, highest first, with ties broken by motif sequence.

    Returns:
        dict: motif -> ID (1 for the first), in ID order.
    """
    ordered_motifs = sorted(motif_sample_counts, key=lambda motif: (-motif_sample_counts[motif], motif))
    return {motif: index + 1 for index, motif in enumerate(ordered_motifs)}


def format_motif_composition(motif_sample_counts, motif_pair_sample_counts, motif_ids):
    """Formats motif and motif pair sample counts as the MotifComposition and MotifPairs values.

    Args:
        motif_sample_counts (dict): motif -> number of samples in which it was observed.
        motif_pair_sample_counts (dict): (motif, next motif) -> number of samples in which the pair was observed.
        motif_ids (dict): get_motif_ids(motif_sample_counts).

    Returns:
        tuple: (MotifComposition value, MotifPairs value), e.g. ("1:GCG:23,2:GCA:20", "[1][1]:23,[1][2]:18"), or two
            empty strings if there are no counts.
    """
    motifs_value = ",".join(
        f"{motif_id}:{motif}:{motif_sample_counts[motif]}" for motif, motif_id in motif_ids.items())
    ordered_pairs = sorted(motif_pair_sample_counts.items(),
                           key=lambda item: (motif_ids[item[0][0]], motif_ids[item[0][1]]))
    pairs_value = ",".join(
        f"[{motif_ids[motif]}][{motif_ids[next_motif]}]:{sample_count}"
        for (motif, next_motif), sample_count in ordered_pairs)
    return motifs_value, pairs_value


def parse_motif_composition_columns(motifs_value, pairs_value, tsv_path):
    """Parses MotifComposition and MotifPairs values back into sample counts keyed by motif sequence.

    Returns:
        tuple: (motif -> number of samples, (motif, next motif) -> number of samples)

    Raises:
        ValueError: if a value isn't in this script's format, as in TSVs made by its earlier versions.
    """
    motif_sample_counts, motif_pair_sample_counts, motif_by_id = {}, {}, {}
    for entry in filter(None, motifs_value.split(",")):
        fields = entry.split(":")
        if len(fields) != 3:
            raise ValueError(f"{tsv_path} has a MotifComposition entry, {entry}, that isn't \"<motif ID>:<motif>:"
                             f"<number of samples>\". {EARLIER_MOTIF_COMPOSITION_FORMAT_EXPLANATION}")
        motif_id, motif, sample_count = fields
        motif_by_id[motif_id] = motif
        motif_sample_counts[motif] = int(sample_count)
    for entry in filter(None, pairs_value.split(",")):
        fields = entry.split(":")
        match = MOTIF_PAIR_KEY_REGEX.fullmatch(fields[0])
        if len(fields) != 2 or match is None:
            raise ValueError(f"{tsv_path} has a MotifPairs entry, {entry}, that isn't \"[<motif ID>][<motif ID>]:"
                             f"<number of samples>\". {EARLIER_MOTIF_COMPOSITION_FORMAT_EXPLANATION}")
        motif_pair_sample_counts[(motif_by_id[match.group(1)], motif_by_id[match.group(2)])] = int(fields[1])
    return motif_sample_counts, motif_pair_sample_counts


def format_motif_occurrence_distributions(motif_occurrence_distributions, motif_ids):
    """Formats each motif's {number of occurrences in a sample: number of samples} as the MotifOccurrenceDistributions
    value, e.g. "1:GCG:10x:5,11x:12,235x:6|2:GCA:1x:18,2x:2", with the motifs in the order and with the IDs of
    motif_ids (get_motif_ids)."""
    return "|".join(f"{motif_id}:{motif}:{format_histogram(motif_occurrence_distributions[motif])}"
                    for motif, motif_id in motif_ids.items())


def parse_motif_occurrence_distributions(motif_occurrence_distributions_value):
    """Parses "1:GCG:10x:5,11x:12|2:GCA:1x:18" into {"GCG": {10: 5, 11: 12}, "GCA": {1: 18}}."""
    motif_occurrence_distributions = {}
    for entry in filter(None, motif_occurrence_distributions_value.split("|")):
        _, motif, distribution = entry.split(":", 2)
        motif_occurrence_distributions[motif] = parse_histogram(distribution)
    return motif_occurrence_distributions


def format_allele_sequences(allele_sequence_counts):
    """Formats consensus sequence counts as "SEQUENCE:count,...", most common first, ties by sequence."""
    if not allele_sequence_counts:
        return ""
    return ",".join(f"{sequence}:{count}" for sequence, count in sorted(
        allele_sequence_counts.items(), key=lambda item: (-item[1], item[0])))


def parse_allele_sequences(allele_sequences_value):
    """Parses "SEQUENCE:count,..." into {sequence: count}."""
    allele_sequence_counts = {}
    for entry in filter(None, allele_sequences_value.split(",")):
        sequence, count = entry.rsplit(":", 1)
        add_count(allele_sequence_counts, sequence, int(count))
    return allele_sequence_counts


def json_tail_command_line_has_motif_composition(input_spec):
    """Returns whether the RunInfo.CommandLine of an ExpansionHunter JSON from an older version, whose RunInfo is one
    of its last records, has --output-motif-composition.

    The whole file is decompressed to reach RunInfo, but it isn't parsed.
    """
    return run_info_has_motif_composition(
        parse_records_from_run_info_on(read_tail_from_run_info(input_spec.path)).get("RunInfo"))


def iterate_expansion_hunter_records(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each repeat record of an ExpansionHunter output JSON, in file order.

    The JSON is streamed, so the file is never fully loaded. stats records whether its RunInfo.CommandLine has
    --output-motif-composition (COMMAND_LINE_HAS_MOTIF_COMPOSITION), which in JSONs from older ExpansionHunter versions
    is only known once the whole file has been read, and whether it was read to its end (READ_TO_THE_END).
    """
    is_male = input_spec.sex.lower() == "male"
    if input_spec.command_line_has_motif_composition:
        stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] = 1

    def record_motif_composition_of_run_info_after_the_loci(run_info):
        if run_info_has_motif_composition(run_info):
            stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] = 1

    # The command line is only unknown in JSONs from older ExpansionHunter versions, whose RunInfo comes after the loci
    for locus_id, variant_id, locus_record, variant_record in iterate_repeat_records_of_json(
            input_spec.path, input_spec.command_line_has_motif_composition is None,
            record_motif_composition_of_run_info_after_the_loci):
        if settings.locus_ids_to_include is not None and locus_id not in settings.locus_ids_to_include:
            continue
        summary = summarize_expansion_hunter_record(variant_record, locus_record.get("AlleleCount"),
                                                    input_spec.sample_id, is_male, settings, stats)
        if summary.called_allele_count and sum(variant.get("VariantType") == "Repeat"
                                               for variant in locus_record.get("Variants", {}).values()) == 1:
            stats[WITH_GENOTYPE_AT_LOCUS_WITH_ONE_REPEAT] += 1
        yield get_genomic_sort_key(locus_id, variant_id, summary.reference_region), summary
    stats[READ_TO_THE_END] = 1


def write_sorted_tsv_from_expansion_hunter_json(input_spec, output_tsv_path, settings, stats, temp_dir):
    """Sorts a seeking or streaming mode JSON (loci sorted by LocusId) into a single-sample TSV in genomic order, with
    memory use of about one chromosome's rows (see write_lines_in_genomic_order)."""
    with gzip.open(output_tsv_path, "wt", compresslevel=TEMPORARY_FILE_COMPRESSION_LEVEL) as f:
        f.write("\t".join(get_output_columns(settings)) + "\n")
        write_lines_in_genomic_order(
            ((key, format_row(key, summary, settings) + "\n")
             for key, summary in iterate_expansion_hunter_records(input_spec, settings, stats)),
            f, temp_dir, get_sort_key_of_row)


def read_tsv_header(tsv_path):
    with open_text(tsv_path) as f:
        return f.readline().rstrip("\n").split("\t")


def has_motif_composition_columns(tsv_columns):
    return all(column in tsv_columns for column in MOTIF_COMPOSITION_COLUMNS)


def iterate_multisample_tsv(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each row of a TSV written by this script.

    Raises:
        ValueError: if its rows are not sorted.
    """
    return iterate_checking_that_keys_increase(parse_multisample_tsv_rows(input_spec, settings, stats),
                                               input_spec.path)


def parse_multisample_tsv_rows(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each row of a TSV written by this script, in file order, leaving out the
    loci not in settings.locus_ids_to_include."""
    for row in iterate_tsv_rows(input_spec.path):
        if settings.locus_ids_to_include is not None and row["LocusId"] not in settings.locus_ids_to_include:
            continue
        stats[TSV_ROWS] += 1
        summary = RepeatSummary(row["Motif"], row["ReferenceRegion"], settings.track_outlier_sample_ids)
        summary.sample_count = int(row["SampleCount"])
        summary.called_allele_count = int(row["CalledAlleleCount"])
        summary.all_allele_histogram = parse_histogram(row["AllAlleleHistogram"])
        summary.short_allele_histogram = parse_histogram(row["ShortAlleleHistogram"])
        summary.hemizygous_allele_histogram = parse_histogram(row["HemizygousAlleleHistogram"])
        summary.biallelic_histogram = parse_biallelic_histogram(row["BiallelicHistogram"])
        if settings.track_outlier_sample_ids:
            summary.all_allele_outlier_sample_ids = parse_outlier_sample_ids(row["OutlierSampleIds_AllAlleles"])
            summary.short_allele_outlier_sample_ids = parse_outlier_sample_ids(row["OutlierSampleIds_ShortAlleles"])
            summary.hemizygous_allele_outlier_sample_ids = parse_outlier_sample_ids(
                row["OutlierSampleIds_HemizygousAlleles"])
        # A TSV made without the motif composition or AlleleSequences columns adds nothing to them
        # (check_optional_columns warns about it). One with motif composition columns has the
        # MotifOccurrenceDistributions column too when it is tracked (check_tsv_input_columns checks).
        if settings.track_motif_composition:
            summary.motif_sample_counts, summary.motif_pair_sample_counts = parse_motif_composition_columns(
                row.get("MotifComposition", ""), row.get("MotifPairs", ""), input_spec.path)
            if settings.track_motif_occurrence_distributions:
                summary.motif_occurrence_distributions = parse_motif_occurrence_distributions(
                    row.get(MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN, ""))
        if settings.track_allele_sequences:
            summary.allele_sequence_counts = parse_allele_sequences(row.get(ALLELE_SEQUENCES_COLUMN, ""))
        yield get_genomic_sort_key(row["LocusId"], row["VariantId"], row["ReferenceRegion"]), summary


def get_output_columns(settings):
    columns = list(BASE_COLUMNS)
    if settings.track_outlier_sample_ids:
        columns += [MAX_OUTLIER_SAMPLE_IDS_COLUMN] + OUTLIER_COLUMNS
    if settings.track_motif_composition:
        columns += MOTIF_COMPOSITION_COLUMNS
        if settings.track_motif_occurrence_distributions:
            columns.append(MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN)
    if settings.track_allele_sequences:
        columns.append(ALLELE_SEQUENCES_COLUMN)
    return columns


def format_row(key, summary, settings):
    """Returns the output TSV line (without the newline) for one repeat."""
    locus_id, variant_id = key[3], key[4]
    row = [
        locus_id, variant_id, summary.reference_region, summary.motif, str(summary.sample_count),
        str(summary.called_allele_count),
        format_histogram(summary.all_allele_histogram), format_histogram(summary.short_allele_histogram),
        format_histogram(summary.hemizygous_allele_histogram), format_biallelic_histogram(summary.biallelic_histogram),
    ]
    if settings.track_outlier_sample_ids:
        row += [
            str(settings.n_outlier_sample_ids),
            get_outlier_sample_ids_string(summary.all_allele_outlier_sample_ids, summary.all_allele_histogram,
                                          settings.n_outlier_sample_ids),
            get_outlier_sample_ids_string(summary.short_allele_outlier_sample_ids, summary.short_allele_histogram,
                                          settings.n_outlier_sample_ids),
            get_outlier_sample_ids_string(summary.hemizygous_allele_outlier_sample_ids,
                                          summary.hemizygous_allele_histogram, settings.n_outlier_sample_ids),
        ]
    if settings.track_motif_composition:
        motif_ids = get_motif_ids(summary.motif_sample_counts or {})
        row += format_motif_composition(summary.motif_sample_counts or {}, summary.motif_pair_sample_counts or {},
                                        motif_ids)
        if settings.track_motif_occurrence_distributions:
            row.append(format_motif_occurrence_distributions(summary.motif_occurrence_distributions or {}, motif_ids))
    if settings.track_allele_sequences:
        row.append(format_allele_sequences(summary.allele_sequence_counts))
    return "\t".join(row)


def merge_inputs(input_specs, output_tsv_path, settings, temp_dir, show_progress_bar=False,
                 compression_level=OUTPUT_COMPRESSION_LEVEL):
    """Merges the inputs one repeat at a time into a TSV with the columns get_output_columns(settings) gives, leaving
    out the repeats of the loci in settings.locus_ids_to_exclude.

    Seeking and streaming mode JSONs are first sorted into temporary TSVs in temp_dir. With settings.max_rows_to_write,
    the merge stops once that many rows are written.

    Args:
        input_specs (list): InputSpec of each input.
        output_tsv_path (str): where to write the merged TSV (gzipped).
        settings (MergeSettings): the run's settings.
        temp_dir (str): directory for temporary files.
        show_progress_bar (bool): whether to show a progress bar of rows written.
        compression_level (int): gzip compression level of the output.

    Returns:
        tuple: (list of per-input collections.Counter stats, in input_specs order; number of rows written; largest
            SampleCount of the rows written; set of the LocusIds whose repeats were left out)
    """
    stats_by_input = [collections.Counter() for _ in input_specs]
    left_out_locus_ids = set()
    sorted_tsv_paths = []
    iterators = []
    try:
        for input_spec, stats in zip(input_specs, stats_by_input):
            if input_spec.kind == "tsv":
                iterators.append(iterate_multisample_tsv(input_spec, settings, stats))
            elif input_spec.loci_sorted_by_locus_id:
                file_descriptor, sorted_tsv_path = tempfile.mkstemp(suffix=".tsv.gz", dir=temp_dir)
                os.close(file_descriptor)
                sorted_tsv_paths.append(sorted_tsv_path)
                write_sorted_tsv_from_expansion_hunter_json(input_spec, sorted_tsv_path, settings, stats, temp_dir)
                iterators.append(iterate_multisample_tsv(tsv_input_spec(sorted_tsv_path), settings,
                                                         collections.Counter()))
            else:
                iterators.append(iterate_in_genomic_order(
                    iterate_expansion_hunter_records(input_spec, settings, stats), input_spec.path, temp_dir))
        merged = merge_records_by_key(iterators, lambda summary, other_summary, description: merge_repeat_summaries(
            summary, other_summary, settings.n_outlier_sample_ids, description))
        if show_progress_bar:
            merged = tqdm.tqdm(merged, unit=" rows", unit_scale=True, desc="Merging")

        num_rows = 0
        largest_sample_count = 0
        with gzip.open(output_tsv_path, "wt", compresslevel=compression_level) as f:
            f.write("\t".join(get_output_columns(settings)) + "\n")
            for key, summary in merged:
                # key[3] is the LocusId
                if key[3] in settings.locus_ids_to_exclude:
                    left_out_locus_ids.add(key[3])
                    continue
                f.write(format_row(key, summary, settings) + "\n")
                num_rows += 1
                largest_sample_count = max(largest_sample_count, summary.sample_count)
                if num_rows == settings.max_rows_to_write:
                    break
    finally:
        # Closing the iterators closes their files, which a merge stopped by -n leaves part-read.
        for iterator in iterators:
            iterator.close()
        remove_files(sorted_tsv_paths)
    return stats_by_input, num_rows, largest_sample_count, left_out_locus_ids


def merge_input_shard(task):
    """Merges one shard of the inputs into a temporary TSV; the worker-process entry point for --threads.

    Args:
        task (tuple): (shard index, list of InputSpec, output TSV path, MergeSettings, temp dir)

    Returns:
        tuple: (shard index, list of per-input stats, set of the LocusIds whose repeats were left out)
    """
    shard_index, input_specs, output_tsv_path, settings, temp_dir = task
    stats_by_input, _, _, left_out_locus_ids = merge_inputs(input_specs, output_tsv_path, settings, temp_dir,
                                                            compression_level=TEMPORARY_FILE_COMPRESSION_LEVEL)
    return shard_index, stats_by_input, left_out_locus_ids


def json_has_motif_composition_in_first_records(json_path, max_records_to_check):
    """Returns whether an ExpansionHunter JSON has MotifComposition in its first max_records_to_check repeat records.

    False doesn't rule out motif composition further on (ExpansionHunter's loci-with-non-ref-motifs and
    loci-with-known-motifs modes put it on few loci); main then checks a JSON's command line too.
    """
    with open_binary(json_path) as f:
        return any("MotifComposition" in variant_record for _, _, _, variant_record in itertools.islice(
            iterate_repeat_records(f), max_records_to_check))


def merge_shards(input_specs_by_shard, temp_dir, settings, threads, show_progress_bar):
    """Merges each shard of the inputs into a temporary TSV, in up to `threads` worker processes.

    A single shard is merged in this process straight into its final form (output compression level).

    Returns:
        tuple: (list of the shards' TSV paths, list of each shard's per-input stats, for a single shard the number of
            rows written and their largest SampleCount, otherwise None and None, and the set of the LocusIds whose
            repeats were left out)
    """
    num_workers = min(threads, len(input_specs_by_shard))
    print(f"Merging {format_count_with_noun(sum(map(len, input_specs_by_shard)), 'input file')}"
          + (f" using {num_workers:,d} threads" if num_workers >= 2 else ""))
    if len(input_specs_by_shard) == 1:
        merged_path = os.path.join(temp_dir, "merged.tsv.gz")
        stats_by_input, num_rows, largest_sample_count, left_out_locus_ids = merge_inputs(
            input_specs_by_shard[0], merged_path, settings, temp_dir, show_progress_bar)
        return [merged_path], [stats_by_input], num_rows, largest_sample_count, left_out_locus_ids

    shard_tsv_paths = get_shard_tsv_paths(temp_dir, len(input_specs_by_shard))
    tasks = [(index, shard, shard_tsv_path, settings, temp_dir)
             for index, (shard, shard_tsv_path) in enumerate(zip(input_specs_by_shard, shard_tsv_paths))]
    stats_by_shard = [None] * len(tasks)
    left_out_locus_ids = set()
    for shard_index, stats_by_input, shard_left_out_locus_ids in map_in_worker_processes(
            merge_input_shard, tasks, num_workers, show_progress_bar, "shards"):
        stats_by_shard[shard_index] = stats_by_input
        left_out_locus_ids |= shard_left_out_locus_ids
    return shard_tsv_paths, stats_by_shard, None, None, left_out_locus_ids


def merge_tsvs(tsv_paths, output_tsv_path, settings, temp_dir, show_progress_bar):
    """Merges TSVs written by this script into one, in rounds if there are more than one merge can read at once.

    The input TSVs are deleted once merged.

    Returns:
        tuple: (number of rows written, largest SampleCount of those rows)
    """
    # the shards' TSVs have no repeats left to leave out
    def merge_tsv_files(paths, output_path, is_final_merge):
        _, num_rows, largest_sample_count, _ = merge_inputs(
            [tsv_input_spec(path) for path in paths], output_path, settings, temp_dir,
            show_progress_bar=show_progress_bar and is_final_merge,
            compression_level=OUTPUT_COMPRESSION_LEVEL if is_final_merge else TEMPORARY_FILE_COMPRESSION_LEVEL)
        return num_rows, largest_sample_count

    return merge_files_in_rounds(tsv_paths, output_tsv_path, temp_dir, merge_tsv_files)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output-prefix",
                        help="Writes <prefix>.tsv.gz (a local path). Default: the input file name without its suffix "
                             "for one input, otherwise combined.<N>_samples.<R>_rows in the current directory, N being "
                             "the number of samples and R the number of rows written")
    parser.add_argument("-n", "--max-rows-to-write", type=int,
                        help="Write only the first N rows: the first N repeats in genomic order. Useful for testing.")
    parser.add_argument("-L", "--loci", action="append", metavar="LOCUS_ID",
                        help="Keep only the repeats of this LocusId (every repeat of a locus with several, such as "
                             "HTT). Can be given more than once, or as a comma-separated list: -L HTT -L ATXN1, or "
                             "-L HTT,ATXN1")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of worker processes. Each merges a share of the inputs into a temporary TSV "
                             "next to the output, and those are then merged (default: %(default)s)")
    parser.add_argument("--add-outlier-sample-id-columns", action="store_true",
                        help=f"Add the {MAX_OUTLIER_SAMPLE_IDS_COLUMN}, OutlierSampleIds_AllAlleles, "
                             f"OutlierSampleIds_ShortAlleles and OutlierSampleIds_HemizygousAlleles columns")
    parser.add_argument("--n-outlier-sample-ids", type=int, default=10,
                        help="Number of samples with the largest alleles to list in each outlier column (default: "
                             "%(default)s)")
    parser.add_argument("--skip-motif-composition-columns", action="store_true",
                        help="Leave out the MotifComposition and MotifPairs columns even if the inputs "
                             "have motif composition")
    parser.add_argument("--output-motif-occurrence-distributions", action="store_true",
                        help=f"Add the {MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN} column: for each motif, the number of "
                             f"samples with each number of occurrences of it, such as \"1:CAG:10x:23,11x:23|"
                             f"2:CAA:1x:40\". TSV inputs with motif composition columns must have it too.")
    parser.add_argument("--add-consensus-sequence-column", action="store_true",
                        help="Add the AlleleSequences column, from the inputs' ConsensusSequences")
    parser.add_argument("--allow-overlapping-samples", action="store_true",
                        help="Allow a SampleId to appear in more than one JSON input. Intended to be used when inputs "
                             "contain data for non-overlapping sets of loci (such as per-chromosome runs of the same "
                             "samples)")
    parser.add_argument("--samples-to-exclude",
                        help="A file of samples to leave out: a TSV whose header has a SampleId column, such as the "
                             "samples_to_exclude.tsv that start_multisample_qc_dashboard.py exports, or one SampleId "
                             "per line. "
                             "Only JSON inputs can be left out.")
    parser.add_argument("--loci-to-exclude",
                        help="A file of loci to leave out of the output, all repeats of each: a TSV whose header has a "
                             "LocusId column, such as the loci_to_exclude.tsv that start_multisample_qc_dashboard.py "
                             "exports, or one LocusId per line. Loci are matched by LocusId only.")
    parser.add_argument("--input-paths-file",
                        help="A file listing more input paths, one per line, for cohorts too large to list on the "
                             "command line")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show progress bars")
    parser.add_argument("input_paths", nargs="*",
                        help="ExpansionHunter output JSON files (.json, .json.gz, .json.bgz) and/or TSVs written by "
                             "earlier runs of this script (.tsv, .tsv.gz, .tsv.bgz), local or gs:// paths")
    args = parser.parse_args(argv)

    args.input_paths = read_input_paths(parser, args.input_paths, args.input_paths_file, JSON_SUFFIXES + TSV_SUFFIXES)
    for option, path in [("--samples-to-exclude", args.samples_to_exclude),
                         ("--loci-to-exclude", args.loci_to_exclude)]:
        if path and not file_exists(path):
            parser.error(f"{option} {path} not found")
    if args.samples_to_exclude and any(path.endswith(TSV_SUFFIXES) for path in args.input_paths):
        parser.error("--samples-to-exclude can only leave out JSON inputs, since a TSV from an earlier run has its "
                     "samples merged already. Pass it to the runs that make those TSVs.")
    if args.output_prefix:
        check_output_prefix(parser, args.output_prefix)
    if args.max_rows_to_write is not None and args.max_rows_to_write < 1:
        parser.error("-n must be at least 1")
    if args.loci is not None:
        args.loci = frozenset(locus_id for value in args.loci for locus_id in value.split(",") if locus_id)
        if not args.loci:
            parser.error("-L was given without a LocusId")
    if args.output_motif_occurrence_distributions and args.skip_motif_composition_columns:
        parser.error("--output-motif-occurrence-distributions can't be combined with --skip-motif-composition-columns")
    if args.n_outlier_sample_ids < 2:
        parser.error("--n-outlier-sample-ids must be at least 2: the outlier walk stops at the first allele size that "
                     "this many samples have, so with 1 the outlier columns would always be empty")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if not args.output_prefix and len(args.input_paths) == 1:
        # Next to a local input; in the current directory for a gs:// one. (Without -o, several inputs get a name
        # from their number of samples and rows once they have been merged: see get_default_combined_output_prefix.)
        input_path = args.input_paths[0]
        args.output_prefix = strip_suffix(os.path.basename(input_path) if is_gcs_path(input_path) else input_path,
                                          JSON_SUFFIXES + TSV_SUFFIXES)
    if args.output_prefix and output_would_overwrite_an_input(args.output_prefix + ".tsv.gz", args.input_paths):
        parser.error(f"The output, {args.output_prefix}.tsv.gz, would overwrite an input; pass a different -o")
    return args


def get_default_combined_output_prefix(input_specs, largest_sample_count, num_rows):
    """Returns the output prefix for several inputs when -o isn't given: combined.<N>_samples.<R>_rows, in the current
    directory, R being the number of rows written.

    A TSV input doesn't record which samples went into it, so N is the larger of the number of distinct SampleIds of
    the JSON inputs and the largest SampleCount among the output's rows (the number of samples at the most complete
    repeat). That counts both batches of different samples and batches of the same samples split by chromosome. With
    -n, the largest SampleCount is taken over the rows written.

    Args:
        input_specs (list): InputSpec of each input.
        largest_sample_count (int): largest SampleCount among the output's rows.
        num_rows (int): number of rows written to the output.
    """
    num_samples = max(len({input_spec.sample_id for input_spec in input_specs if input_spec.kind == "json"}),
                      largest_sample_count)
    return (f"combined.{num_samples}_sample{'' if num_samples == 1 else 's'}"
            f".{num_rows}_row{'' if num_rows == 1 else 's'}")


def read_outlier_limit_of_tsv(tsv_path, header):
    """Returns the MaxOutlierSampleIds of a TSV's first row, or None if it has no rows."""
    with open_text(tsv_path) as f:
        next(f)
        first_row = f.readline()
    if not first_row:
        return None
    return int(first_row.rstrip("\n").split("\t")[header.index(MAX_OUTLIER_SAMPLE_IDS_COLUMN)])


def check_tsv_input_columns(tsv_path, header, args):
    """Checks a TSV input has the columns of a TSV written by this script and, if this run adds outlier columns, has
    them too, with an outlier limit at least as large as this run's, and if this run adds the
    MotifOccurrenceDistributions column, has it too when it has motif composition columns.

    Args:
        tsv_path (str): the TSV.
        header (list): its columns.
        args (argparse.Namespace): the command line arguments.

    Raises:
        ValueError: if its columns are not ones this script writes, it lacks the outlier columns or the
            MotifOccurrenceDistributions column this run needs, or its outlier columns were made with a smaller
            --n-outlier-sample-ids.
    """
    headers_this_script_writes = [
        get_output_columns(MergeSettings(
            n_outlier_sample_ids=None, track_outlier_sample_ids=track_outlier_sample_ids,
            track_motif_composition=track_motif_composition,
            track_motif_occurrence_distributions=track_motif_occurrence_distributions,
            track_allele_sequences=track_allele_sequences, max_rows_to_write=None))
        for track_outlier_sample_ids in (False, True) for track_motif_composition in (False, True)
        for track_motif_occurrence_distributions in (False, True) for track_allele_sequences in (False, True)]
    if header not in headers_this_script_writes:
        raise ValueError(f"{tsv_path} doesn't have the columns of a TSV written by this script: {', '.join(header)}")
    if args.output_motif_occurrence_distributions and has_motif_composition_columns(header) \
            and MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN not in header:
        raise ValueError(f"{tsv_path} has motif composition columns but no {MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN} "
                         f"column, so its samples' numbers of motif occurrences are unknown. Remake it with "
                         f"--output-motif-occurrence-distributions.")
    if args.add_outlier_sample_id_columns:
        # The header is one this script writes, so it has either all the outlier columns or none of them
        if MAX_OUTLIER_SAMPLE_IDS_COLUMN not in header:
            raise ValueError(f"{tsv_path} has no outlier sample ID columns. It was made without "
                             f"--add-outlier-sample-id-columns, so its outlier samples are unknown.")
        outlier_limit = read_outlier_limit_of_tsv(tsv_path, header)
        if outlier_limit is not None and outlier_limit < args.n_outlier_sample_ids:
            raise ValueError(
                f"{tsv_path} was made with --n-outlier-sample-ids {outlier_limit}, smaller than this run's "
                f"{args.n_outlier_sample_ids}, so it left out samples this run would list. Remake it with "
                f"--n-outlier-sample-ids {args.n_outlier_sample_ids} or more.")


def collect_input_specs(args):
    """Reads each JSON's SampleId and Sex and checks each TSV's columns.

    Returns:
        list: InputSpec of each input.

    Raises:
        ValueError: if two JSON inputs have the same SampleId (without --allow-overlapping-samples), or a TSV lacks
            the outlier columns or the MotifOccurrenceDistributions column this run needs.
    """
    json_paths = [path for path in args.input_paths if path.endswith(JSON_SUFFIXES)]
    header_by_json_path = dict(zip(json_paths, map_in_threads(read_expansion_hunter_json_header, json_paths)))
    input_specs = []
    json_path_by_sample_id = {}
    for path in args.input_paths:
        if path.endswith(JSON_SUFFIXES):
            sample_id, sex, loci_sorted_by_locus_id, command_line = header_by_json_path[path]
            # None if RunInfo comes after the loci, so its command line isn't known yet
            command_line_has_motif_composition = (
                None if command_line is None else MOTIF_COMPOSITION_OPTION in command_line)
            if sample_id in json_path_by_sample_id and not args.allow_overlapping_samples:
                raise ValueError(
                    f"{path} and {json_path_by_sample_id[sample_id]} have the same SampleId, {sample_id}. If the "
                    f"inputs cover different loci, such as per-chromosome runs of the same sample, pass "
                    f"--allow-overlapping-samples.")
            json_path_by_sample_id[sample_id] = path
            input_specs.append(InputSpec(
                path, "json", sample_id, sex, loci_sorted_by_locus_id, command_line_has_motif_composition))
        else:
            header = read_tsv_header(path)
            check_tsv_input_columns(path, header, args)
            input_specs.append(tsv_input_spec(path, header))
    return input_specs


def check_optional_columns(input_specs, stats_by_input, args):
    """Decides whether the output gets motif composition columns, and checks the optional data is in every input.

    Returns:
        tuple: (whether to write the motif composition columns, list of warnings to print). The warnings are about
            inputs without motif composition when others have it (without --skip-motif-composition-columns), inputs
            without consensus sequences when --add-consensus-sequence-column was given, and JSONs whose motif
            composition couldn't be checked because of -n.
    """
    inputs_with_motif_composition, inputs_without_motif_composition, inputs_without_consensus_sequences = [], [], []
    inputs_with_unknown_motif_composition = []
    for input_spec, stats in zip(input_specs, stats_by_input):
        if input_spec.kind == "json":
            # A run with --output-motif-composition loci-with-non-ref-motifs (or loci-with-known-motifs) can have no
            # MotifComposition record at all, so the command line counts too.
            has_motif_composition = stats[WITH_MOTIF_COMPOSITION] > 0 or stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] > 0
            # ExpansionHunter writes ConsensusSequences only for called repeats at loci with one repeat, so a JSON
            # without any such call (or with only such loci left out by -L or -n) has none to contribute rather than
            # missing them.
            if stats[WITH_CONSENSUS_SEQUENCES] == 0 and stats[WITH_GENOTYPE_AT_LOCUS_WITH_ONE_REPEAT] > 0:
                inputs_without_consensus_sequences.append(input_spec.path)
            # With -n, a JSON from an older ExpansionHunter version may not have been read as far as its RunInfo, so
            # one without MotifComposition records so far can't be told either way.
            if not has_motif_composition and input_spec.command_line_has_motif_composition is None \
                    and not stats[READ_TO_THE_END]:
                inputs_with_unknown_motif_composition.append(input_spec.path)
                continue
        else:
            has_motif_composition = has_motif_composition_columns(input_spec.tsv_columns)
            if args.add_consensus_sequence_column and ALLELE_SEQUENCES_COLUMN not in input_spec.tsv_columns:
                inputs_without_consensus_sequences.append(input_spec.path)
        (inputs_with_motif_composition if has_motif_composition else inputs_without_motif_composition).append(
            input_spec.path)

    warnings = []
    if not args.skip_motif_composition_columns and inputs_with_motif_composition and inputs_without_motif_composition:
        warnings.append(
            f"Motif composition is missing from "
            f"{format_count_with_noun(len(inputs_without_motif_composition), 'input')} (e.g. "
            f"{inputs_without_motif_composition[0]}) but present in "
            f"{format_count_with_noun(len(inputs_with_motif_composition), 'other input')} (e.g. "
            f"{inputs_with_motif_composition[0]}), so the motif composition columns count only the samples in the "
            f"inputs that have it. Run ExpansionHunter with --output-motif-composition on the others to count every "
            f"sample, or pass --skip-motif-composition-columns to leave the columns out. (A JSON counts as having "
            f"motif composition if any of its repeat records has a MotifComposition field or its RunInfo.CommandLine "
            f"has --output-motif-composition; a TSV, if it has the motif composition columns.)")
    if args.add_consensus_sequence_column and inputs_without_consensus_sequences:
        warnings.append(
            f"--add-consensus-sequence-column was given, but consensus sequences are missing from "
            f"{format_count_with_noun(len(inputs_without_consensus_sequences), 'input')} (e.g. "
            f"{inputs_without_consensus_sequences[0]}), so the "
            f"{ALLELE_SEQUENCES_COLUMN} column counts only the samples in the other inputs. (A JSON has none if it has "
            f"called repeats but no ConsensusSequences, as when ExpansionHunter is run with "
            f"--dont-output-consensus-sequences; a TSV, if it was made without --add-consensus-sequence-column.)")
    if not args.skip_motif_composition_columns and inputs_with_unknown_motif_composition:
        warnings.append(
            f"Motif composition couldn't be checked for "
            f"{format_count_with_noun(len(inputs_with_unknown_motif_composition), 'JSON input')} from older "
            f"ExpansionHunter versions (e.g. {inputs_with_unknown_motif_composition[0]}): the rows read have no "
            f"MotifComposition, and -n stopped before RunInfo, which these versions write at the end. The motif "
            f"composition columns may be missing, or count only some of the samples.")
    return not args.skip_motif_composition_columns and bool(inputs_with_motif_composition), warnings


def print_input_stats(input_specs, stats_by_input):
    for input_spec, stats in zip(input_specs, stats_by_input):
        if input_spec.kind == "tsv":
            print(f"  {os.path.basename(input_spec.path)}: {format_count_with_noun(stats[TSV_ROWS], 'row')}")
            continue
        details = ", ".join(f"{count:,d} {description}" for description, count in stats.items()
                            if description not in (REPEAT_RECORDS, COMMAND_LINE_HAS_MOTIF_COMPOSITION, READ_TO_THE_END,
                                                   WITH_GENOTYPE_AT_LOCUS_WITH_ONE_REPEAT))
        print(f"  {os.path.basename(input_spec.path)} (sample {input_spec.sample_id}, "
              f"sex {input_spec.sex or 'unknown'}): "
              f"{format_count_with_noun(stats[REPEAT_RECORDS], 'repeat record')}" + (f", {details}" if details else ""))


def main(argv=None):
    set_up_command_line_run()
    args = parse_args(argv)
    json_count = sum(path.endswith(JSON_SUFFIXES) for path in args.input_paths)
    tsv_count = len(args.input_paths) - json_count
    if json_count and tsv_count:
        print(f"Combining {len(args.input_paths):,d} inputs: {format_count_with_noun(tsv_count, 'TSV')} and "
              f"{format_count_with_noun(json_count, 'ExpansionHunter JSON file')}")
    elif json_count:
        print(f"Combining {format_count_with_noun(json_count, 'ExpansionHunter JSON file')}")
    else:
        print(f"Combining {format_count_with_noun(tsv_count, 'TSV file')}")
    input_specs = collect_input_specs(args)
    if args.samples_to_exclude:
        input_specs = leave_out_samples(input_specs, read_ids_to_exclude(args.samples_to_exclude, "SampleId"),
                                        args.samples_to_exclude)
    locus_ids_to_exclude = (read_ids_to_exclude(args.loci_to_exclude, "LocusId") if args.loci_to_exclude
                            else frozenset())
    samples_without_sex = sorted({input_spec.sample_id for input_spec in input_specs
                                  if input_spec.kind == "json" and input_spec.sex.lower() not in ("male", "female")})
    if samples_without_sex:
        print(f"WARNING: no Sex in SampleParameters for {format_count_with_noun(len(samples_without_sex), 'sample')} "
              f"(e.g. {samples_without_sex[0]}). They are counted as female, so only their chrY calls are hemizygous")

    # Whether the output gets motif composition columns is only certain once every input has been read. A TSV's header
    # tells, and so does RunInfo.CommandLine in every JSON whose command line was read with its header, which costs
    # nothing. For the other JSONs (from older ExpansionHunter versions), the first records tell in ExpansionHunter's
    # all-loci mode. Its other modes can put MotifComposition on few loci or none, so if no JSON's command line was
    # read with its header, the smallest JSON's is read, which takes decompressing the whole file, so -n skips it
    # (check_optional_columns then warns about the JSONs whose motif composition it can't tell). Inputs usually agree
    # on motif composition (check_optional_columns warns when they don't), so that one JSON usually speaks for all. The
    # merge is rerun if MotifComposition turns up late in JSONs whose command line wasn't checked or abbreviates
    # --output-motif-composition. With -n, at most N records of each JSON are checked first, since the merge itself
    # reads about that many.
    # Each gs:// file's size takes a request to Google Cloud Storage, so it is looked up once.
    file_size_by_path = dict(zip(args.input_paths, map_in_threads(get_file_size, args.input_paths)))
    track_motif_composition = False
    if not args.skip_motif_composition_columns:
        json_specs = [input_spec for input_spec in input_specs if input_spec.kind == "json"]
        every_json_command_line_is_unknown = bool(json_specs) and all(
            input_spec.command_line_has_motif_composition is None for input_spec in json_specs)
        max_records_to_check = min(MOTIF_COMPOSITION_CHECK_RECORDS, args.max_rows_to_write or math.inf)
        track_motif_composition = (
            any(has_motif_composition_columns(input_spec.tsv_columns) for input_spec in input_specs
                if input_spec.kind == "tsv")
            or any(input_spec.command_line_has_motif_composition for input_spec in json_specs)
            or any(json_has_motif_composition_in_first_records(input_spec.path, max_records_to_check)
                   for input_spec in json_specs if input_spec.command_line_has_motif_composition is None)
            or (every_json_command_line_is_unknown and args.max_rows_to_write is None
                and json_tail_command_line_has_motif_composition(
                    min(json_specs, key=lambda input_spec: file_size_by_path[input_spec.path]))))
    settings = MergeSettings(
        n_outlier_sample_ids=args.n_outlier_sample_ids,
        track_outlier_sample_ids=args.add_outlier_sample_id_columns,
        track_motif_composition=track_motif_composition,
        track_motif_occurrence_distributions=args.output_motif_occurrence_distributions,
        track_allele_sequences=args.add_consensus_sequence_column,
        max_rows_to_write=args.max_rows_to_write,
        locus_ids_to_exclude=locus_ids_to_exclude,
        locus_ids_to_include=args.loci,
    )
    if args.loci is not None:
        print(f"Keeping only the repeats of the {format_count_with_noun(len(args.loci), 'LocusId')} given with -L, so "
              f"the input counts below cover only their records")
    if args.max_rows_to_write is not None:
        print(f"Writing only the first {format_count_with_noun(args.max_rows_to_write, 'row')} (-n), so the input "
              f"counts below cover only the records read for them")
    input_specs_by_shard = split_inputs_into_shards(input_specs, args.threads, file_size_by_path)
    # Without -o, several inputs are named after their number of samples and rows, which are only known once merged,
    # so until then the temporary files go where that output will: in the current directory, named "combined".
    with create_temporary_directory_next_to(args.output_prefix or "combined") as temp_dir:
        while True:
            shard_tsv_paths, stats_by_shard, num_rows, largest_sample_count, left_out_locus_ids = merge_shards(
                input_specs_by_shard, temp_dir, settings, args.threads, args.show_progress_bar)
            ordered_specs = [spec for shard in input_specs_by_shard for spec in shard]
            ordered_stats = [stats for shard_stats in stats_by_shard for stats in shard_stats]
            write_motif_composition_columns, optional_column_warnings = check_optional_columns(
                ordered_specs, ordered_stats, args)
            if not write_motif_composition_columns or settings.track_motif_composition:
                break
            print("Motif composition only appears after the repeat records checked before the merge, so merging again "
                  "with the motif composition columns")
            remove_files(shard_tsv_paths)
            settings = settings._replace(track_motif_composition=True)

        print("Inputs:")
        print_input_stats(ordered_specs, ordered_stats)
        for description in [WITH_ALLELE_COUNT_NOT_MATCHING_SEX, WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE]:
            total = sum(stats[description] for stats in ordered_stats)
            if total:
                print(f"WARNING: {format_count_with_noun(total, 'genotype')} {description}")
        # Printed once the merge is final, so a rerun of the merge doesn't print them twice
        for warning in optional_column_warnings:
            print(f"WARNING: {warning}")
        if args.loci_to_exclude:
            print(f"Left out the repeats of {len(left_out_locus_ids):,d} of the {len(locus_ids_to_exclude):,d} "
                  f"{'locus' if len(locus_ids_to_exclude) == 1 else 'loci'} listed in {args.loci_to_exclude}")
            # -n stops the merge before it reaches most repeats, and -L skips the other loci, so only a full merge can
            # tell which aren't in it
            missing_locus_ids = sorted(locus_ids_to_exclude - left_out_locus_ids)
            if missing_locus_ids and args.max_rows_to_write is None and args.loci is None:
                print(f"WARNING: {len(missing_locus_ids):,d} {'locus' if len(missing_locus_ids) == 1 else 'loci'} "
                      f"listed in {args.loci_to_exclude} {'is' if len(missing_locus_ids) == 1 else 'are'} not in the "
                      f"inputs (e.g. {missing_locus_ids[0]})")

        if len(shard_tsv_paths) == 1:
            merged_path = shard_tsv_paths[0]
        else:
            merged_path = os.path.join(temp_dir, "merged.tsv.gz")
            num_rows, largest_sample_count = merge_tsvs(
                shard_tsv_paths, merged_path, settings, temp_dir, args.show_progress_bar)
        output_tsv_path = (args.output_prefix or get_default_combined_output_prefix(
            ordered_specs, largest_sample_count, num_rows)) + ".tsv.gz"
        # parse_args checked -o and a single input's name; the combined name is only known now. It is a TSV input's
        # own name when that input is an earlier default-named output with the same numbers of samples and rows, e.g.
        # a new sample at loci the input already has, when no row of the input had every one of its samples.
        if output_would_overwrite_an_input(output_tsv_path, args.input_paths):
            raise ValueError(f"The output, {output_tsv_path}, would overwrite an input; pass -o")
        os.replace(merged_path, output_tsv_path)

    if args.loci is not None and args.max_rows_to_write is None:
        # With -L the output is small, so reading it back is quick
        missing_locus_ids = sorted(args.loci - {row["LocusId"] for row in iterate_tsv_rows(output_tsv_path)})
        if missing_locus_ids:
            print(f"WARNING: {format_count_with_noun(len(missing_locus_ids), 'LocusId')} given with -L "
                  f"{'has' if len(missing_locus_ids) == 1 else 'have'} no rows in the output (e.g. "
                  f"{missing_locus_ids[0]})")
    if args.output_motif_occurrence_distributions and not settings.track_motif_composition:
        print(f"WARNING: --output-motif-occurrence-distributions was given, but no input has motif composition, so "
              f"there is no {MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN} column")
    added_columns = [name for name, added in [
        ("outlier sample ID", settings.track_outlier_sample_ids),
        ("motif composition", settings.track_motif_composition),
        (MOTIF_OCCURRENCE_DISTRIBUTIONS_COLUMN,
         settings.track_motif_composition and settings.track_motif_occurrence_distributions),
        (ALLELE_SEQUENCES_COLUMN, settings.track_allele_sequences)] if added]
    print(f"Wrote {format_count_with_noun(num_rows, 'row')} to {output_tsv_path}"
          + (f" (with {', '.join(added_columns)} columns)" if added_columns else ""))


if __name__ == "__main__":
    main()
