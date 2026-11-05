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

    - ExpansionHunter's optimized-streaming mode already writes the loci nearly in genomic order.
    - ExpansionHunter's seeking and streaming modes write the loci sorted by LocusId. This script pre-sorts their JSON output 
      into a temporary single-sample TSV, one chromosome at a time, before processing it.


The output has one row per repeat: for loci with a single repeat (every locus in the TRExplorer catalog) that is one
row per LocusId, and a locus with several repeats (for example HTT and HTT_CCG) gets one
row per repeat. Rows are sorted by position. Columns include:

    LocusId, VariantId, ReferenceRegion, Motif   from the ExpansionHunter output
    LocusStructure               from --variant-catalog if given; otherwise "(<Motif>)*" for loci with one repeat
                                 and empty for loci with compound repeat definitions.
    SampleCount                  samples whose ExpansionHunter output has a record for this locus, regardless of whether 
                                 a genotype was called or not. ExpansionHunter leaves a locus out of a sample's output 
                                 if its --skip-hom-ref or --skip-missing-genotypes options were set, and also skips loci 
                                 that are too wide (more than 2x the sample's read length) or whose motif is too long 
                                 (more than half the read length). A sample without a record for the repeat is not 
                                 counted.
    CalledAlleleCount            called alleles: 2 per diploid genotype, 1 per hemizygous one, 0 per no-call
    AllAlleleHistogram           every called allele, in repeat units: "15x:120,16x:3"
    ShortAlleleHistogram         the shorter allele of each called sample (the only allele of a hemizygous call)
    HemizygousAlleleHistogram    the allele of each hemizygous call: a chrY repeat in any sample, or a chrX repeat in
                                 a male sample (based on SampleParameters.Sex in ExpansionHunter's JSON ouptut)
    BiallelicHistogram           genotype counts, "short/long:count", sorted: "15/15:58,15/16:3". A hemizygous call's
                                 single allele is written twice ("8" counts as "8/8"). Missing genotypes are not counted.

With --add-outlier-sample-id-columns, the MaxOutlierSampleIds column is set to the value of the --n-outlier-sample-ids argument, 
while the OutlierSampleIds_AllAlleles, OutlierSampleIds_ShortAlleles and OutlierSampleIds_HemizygousAlleles columns are set to a list of the
samples with the largest alleles (eg. "20x:sampleA,18x:sampleB"). Allele sizes are taken from largest to smallest, stopping at the 
first allele size that is shared by --n-outlier-sample-ids or more samples, or after more than that many samples have been listed. 

When the inputs have ExpansionHunter's motif composition (enabled by ExpansionHunter's --output-motif-composition option), these columns are also added:

    MotifComposition    "1:GCG:561:64,2:GCA:82:50" (representing motif ID, motif, total number of occurrences across reads, number of reads where the motif was observed)
    MotifPairs          "[1][1]:398:63,[1][2]:82:50" (representing the IDs of a motif and the motif after it, total number of occurrences across reads, number of reads where the motif pair was observed)

Both of these columns sum the locus-level Motifs and MotifPairs records of ExpansionHunter's MotifComposition output field over all samples (ExpansionHunter's
per-allele Allele1/Allele2 output blocks are not used). ExpansionHunter numbers motifs separately in each sample, so motif 1 can be GCG in one sample and GCA in
another: motifs and motif pairs are summed by motif sequence, and both columns renumber motifs by their total occurrences, highest first, with ties broken by motif sequence. 
Either every input must have motif composition or none.

With --add-consensus-sequence-columns, an AlleleSequences column is also added: each distinct consensus allele sequence from
ExpansionHunter's ConsensusSequences field is listed, along with its total haplotype count, most common first, "TGCTGCTGA:86,TGCTGCTGC:56". A
homozygous diploid call's one consensus sequence counts twice, and a hemizygous call's once. Calls without
ConsensusSequences add nothing, so the counts can add up to less than CalledAlleleCount. 

JSON inputs are checked to make sure that no sample id is present in more than one input file. 

Inputs (and --variant-catalog) can be local paths or gs:// paths in Google Cloud Storage, which are streamed rather
than downloaded. Reading gs:// paths needs the google-cloud-storage package and application default credentials
(gcloud auth application-default login). The output is written locally.

Example:
    python3 generate_multisample_tsv.py --threads 4 --show-progress-bar --add-outlier-sample-id-columns \\
        -o batch1 sample1.json.gz sample2.json.gz sample3.json.gz
    python3 generate_multisample_tsv.py --add-outlier-sample-id-columns -o cohort \\
        batch1.tsv.gz batch2.tsv.gz sample4.json.gz
"""

import collections
import functools
import gzip
import heapq
import io
import itertools
import json
import math
import multiprocessing
import operator
import os
import re
import resource
import sys
import tempfile

import configargparse
import ijson
import tqdm

JSON_SUFFIXES = (".json", ".json.gz", ".json.bgz")
TSV_SUFFIXES = (".tsv", ".tsv.gz", ".tsv.bgz")

BASE_COLUMNS = [
    "LocusId",
    "VariantId",
    "ReferenceRegion",
    "Motif",
    "LocusStructure",
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
ALLELE_SEQUENCES_COLUMN = "AlleleSequences"

MOTIF_PAIR_KEY_REGEX = re.compile(r"\[(\d+)\]\[(\d+)\]")
SINGLE_REPEAT_LOCUS_STRUCTURE_REGEX = re.compile(r"\([A-Za-z]+\)\*")
# How far (in bp) a repeat may come after repeats that start later on the same chromosome. ExpansionHunter's
# optimized-streaming mode writes a locus once the reads have passed it, so a locus can come up to about a flank
# length (1,000 bp) after loci that start later; on a full TRExplorer v2.1 catalog the largest such distance is 955 bp.
# The buffer holds every repeat in the window, so this is about 10x the largest distance seen and no more.
REORDER_WINDOW_BP = 10_000
# ExpansionHunter genotypes chrY with one allele in every sample and chrX with one allele in males
# (determineExpectedAlleleCount in ehunter/core/LocusStats.cpp).
CHR_X_NAMES = {"chrX", "X"}
CHR_Y_NAMES = {"chrY", "Y"}
# Repeat records of each JSON checked for MotifComposition before the merge, to decide whether to track it.
# ExpansionHunter's all-loci mode puts a MotifComposition record on every eligible repeat, so one shows up within the
# first few records. Its other modes put it on few loci, so for those the check relies on RunInfo.CommandLine.
MOTIF_COMPOSITION_CHECK_RECORDS = 10_000
# gzip's default level 9 is several times slower than these for little gain on these tables.
OUTPUT_COMPRESSION_LEVEL = 6
TEMPORARY_FILE_COMPRESSION_LEVEL = 1
# Open files a merge needs besides its inputs (output, bucket files of a JSON being sorted, Python's own).
OPEN_FILES_NEEDED_BESIDES_INPUTS = 64
# Each input a merge reads costs about 1 MB of buffers (plus GCS_READ_CHUNK_BYTES for a gs:// input), so larger
# cohorts are merged in groups of at most this many inputs, and the groups' TSVs merged after.
MAX_INPUTS_PER_MERGE = 256
# Bytes fetched per request when reading a gs:// file. Each open gs:// file holds one such chunk in memory.
GCS_READ_CHUNK_BYTES = 2 * 1024 * 1024
# RunInfo's CommandLine tells whether ExpansionHunter was run with --output-motif-composition. In JSONs from older
# ExpansionHunter versions RunInfo is among the last records. It holds the whole command line, which a --locus list of
# thousands of loci can make megabytes long, so while such a JSON is streamed everything from this key on is kept, up
# to MAX_RUN_INFO_BYTES (well past any command line the operating system accepts).
RUN_INFO_KEY = b'"RunInfo"'
MAX_RUN_INFO_BYTES = 64 * 1024 * 1024
# Per-input stats keys. The JSON ones are printed as "<count> <key>" after "<count> repeat records".
REPEAT_RECORDS = "repeat records"
WITH_GENOTYPE = "with genotype"
WITH_MOTIF_COMPOSITION = "with motif composition"
WITH_CONSENSUS_SEQUENCES = "with consensus sequences"
WITH_ALLELE_COUNT_NOT_MATCHING_SEX = "with an allele count that doesn't match the sample's sex"
WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE = "whose number of consensus sequences doesn't match the genotype"
TSV_ROWS = "rows"
# These two describe the input rather than counting records, so they aren't printed. The first is set from
# RunInfo.CommandLine, the second once a JSON has been read to its end (with -n it may not be).
COMMAND_LINE_HAS_MOTIF_COMPOSITION = "ExpansionHunter command line has --output-motif-composition"
READ_TO_THE_END = "read to the end"

# What every merge needs to know about the run, passed to worker processes. max_rows_to_write is -n (None for all).
MergeSettings = collections.namedtuple("MergeSettings", [
    "n_outlier_sample_ids",
    "track_outlier_sample_ids",
    "track_motif_composition",
    "track_allele_sequences",
    "catalog_locus_structures",
    "max_rows_to_write",
])

# One input: kind is "json" or "tsv". sample_id, sex, loci_sorted_by_locus_id (whether the JSON was written by
# ExpansionHunter's seeking or streaming mode) and command_line_has_motif_composition (whether its RunInfo.CommandLine
# has --output-motif-composition, or None if RunInfo couldn't be read before the loci are streamed, as in older
# optimized-streaming output, whose RunInfo is its last record) are set for JSON inputs only.
InputSpec = collections.namedtuple(
    "InputSpec", ["path", "kind", "sample_id", "sex", "loci_sorted_by_locus_id", "command_line_has_motif_composition"])


def tsv_input_spec(path):
    return InputSpec(path, "tsv", None, None, False, None)


class RepeatSummary:
    """Totals for one repeat over some samples: one sample's record, one TSV row, or a merge of those.

    Histograms map allele size (repeat units) to allele count, and the biallelic histogram maps a genotype
    (short allele, long allele) to its number of samples. Outlier sample IDs map allele size to the IDs of up
    to n_outlier_sample_ids samples with that allele, and are None when outlier columns weren't requested. The
    optional counts stay None unless the repeat has data for them: motif counts map a motif to [occurrences, reads],
    motif pair counts map (motif, next motif) to [occurrences, reads], and allele sequence counts map a consensus
    sequence to its haplotype count.
    """

    __slots__ = (
        "motif", "reference_region", "locus_structure", "sample_count", "called_allele_count",
        "all_allele_histogram", "short_allele_histogram", "hemizygous_allele_histogram", "biallelic_histogram",
        "all_allele_outlier_sample_ids", "short_allele_outlier_sample_ids", "hemizygous_allele_outlier_sample_ids",
        "motif_counts", "motif_pair_counts", "allele_sequence_counts",
    )

    def __init__(self, motif, reference_region, locus_structure, track_outlier_sample_ids):
        self.motif = motif
        self.reference_region = reference_region
        self.locus_structure = locus_structure
        self.sample_count = 0
        self.called_allele_count = 0
        self.all_allele_histogram = {}
        self.short_allele_histogram = {}
        self.hemizygous_allele_histogram = {}
        self.biallelic_histogram = {}
        self.all_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.short_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.hemizygous_allele_outlier_sample_ids = {} if track_outlier_sample_ids else None
        self.motif_counts = None
        self.motif_pair_counts = None
        self.allele_sequence_counts = None


def is_gcs_path(path):
    return path.startswith("gs://")


def get_gcs_client():
    """Returns this process's Google Cloud Storage client. Keyed by process ID so that a worker process forked from the
    main one makes its own client rather than sharing the main process's connections."""
    return get_gcs_client_for_process(os.getpid())


@functools.lru_cache(maxsize=None)
def get_gcs_client_for_process(process_id):
    try:
        from google.cloud import storage
    except ImportError as e:
        raise ImportError("Reading gs:// paths needs the google-cloud-storage package (pip3 install "
                          "google-cloud-storage)") from e
    # Reading files needs no project; without project=None the client fails when no default project is configured.
    return storage.Client(project=None)


def get_gcs_blob(path):
    """Returns the google.cloud.storage Blob of a gs://<bucket>/<name> path, or None if there is no such file."""
    bucket_name, _, blob_name = path[len("gs://"):].partition("/")
    return get_gcs_client().bucket(bucket_name).get_blob(blob_name)


def file_exists(path):
    return get_gcs_blob(path) is not None if is_gcs_path(path) else os.path.isfile(path)


def get_file_size(path):
    if not is_gcs_path(path):
        return os.path.getsize(path)
    blob = get_gcs_blob(path)
    if blob is None:
        raise FileNotFoundError(path)
    return blob.size


class GzipFileThatClosesItsSource(gzip.GzipFile):
    """A GzipFile that decompresses an open file and closes that file when it is closed (GzipFile leaves it open)."""

    def __init__(self, source):
        super().__init__(fileobj=source, mode="rb")
        self.source = source

    def close(self):
        try:
            super().close()
        finally:
            self.source.close()


def open_binary(path):
    """Opens a local or gs:// file for reading bytes, decompressing it if it ends with .gz or .bgz."""
    is_compressed = path.endswith((".gz", ".bgz"))
    if not is_gcs_path(path):
        return gzip.open(path, "rb") if is_compressed else open(path, "rb")
    blob = get_gcs_blob(path)
    if blob is None:
        raise FileNotFoundError(path)
    reader = blob.open("rb", chunk_size=GCS_READ_CHUNK_BYTES)
    return GzipFileThatClosesItsSource(reader) if is_compressed else reader


def open_text(path):
    """Opens a local or gs:// file for reading text, decompressing it if it ends with .gz or .bgz."""
    return io.TextIOWrapper(open_binary(path), encoding="utf-8")


def strip_suffix(path, suffixes):
    """Returns path without whichever of suffixes it ends with (or unchanged if none)."""
    for suffix in sorted(suffixes, key=len, reverse=True):
        if path.endswith(suffix):
            return path[:-len(suffix)]
    return path


def remove_files(paths):
    for path in paths:
        if os.path.isfile(path):
            os.remove(path)


def parse_genotype(genotype):
    """Parses an ExpansionHunter repeat genotype such as "15/17" or "15" (hemizygous).

    Returns:
        list: the allele sizes in repeat units, or None if there is no genotype.
    """
    if not genotype:
        return None
    alleles = genotype.split("/")
    if not all(allele.isdigit() for allele in alleles):
        return None
    return [int(allele) for allele in alleles]


def get_chromosome(reference_region):
    return reference_region.rsplit(":", 1)[0]


def get_genomic_sort_key(locus_id, variant_id, reference_region):
    """Sorts by chromosome (1-22, X, Y, M, then others by name), start, end, LocusId and VariantId."""
    chrom, interval = reference_region.rsplit(":", 1)
    start, end = interval.split("-")
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        chrom_key = (0, int(name), "")
    elif name.upper() in ("X", "Y", "M", "MT"):
        chrom_key = (1, ("X", "Y", "M", "MT").index(name.upper()), "")
    else:
        chrom_key = (2, 0, name)
    return chrom_key, int(start), int(end), locus_id, variant_id


def get_sort_key_of_row(line):
    """Returns get_genomic_sort_key of an output TSV line from its LocusId, VariantId and ReferenceRegion columns."""
    locus_id, variant_id, reference_region = line.split("\t", 3)[:3]
    return get_genomic_sort_key(locus_id, variant_id, reference_region)


def add_count(counts_by_key, key, count):
    counts_by_key[key] = counts_by_key.get(key, 0) + count


def add_occurrences_and_reads(counts_by_key, key, occurrences, reads):
    entry = counts_by_key.setdefault(key, [0, 0])
    entry[0] += occurrences
    entry[1] += reads


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


def summarize_expansion_hunter_record(variant_record, sample_id, is_male, locus_structure, settings, stats):
    """Turns one sample's repeat record from an ExpansionHunter JSON into a RepeatSummary.

    Args:
        variant_record (dict): the record under LocusResults.<locus>.Variants.
        sample_id (str): the sample's ID.
        is_male (bool): whether the sample's SampleParameters.Sex is male.
        locus_structure (str): the LocusStructure to report for the repeat.
        settings (MergeSettings): the run's settings.
        stats (collections.Counter): per-input counts to update.

    Returns:
        RepeatSummary: the sample's contribution to the repeat's row.
    """
    summary = RepeatSummary(variant_record["RepeatUnit"], variant_record["ReferenceRegion"], locus_structure,
                            settings.track_outlier_sample_ids)
    summary.sample_count = 1
    stats[REPEAT_RECORDS] += 1

    alleles = parse_genotype(variant_record.get("Genotype"))
    if alleles is not None:
        stats[WITH_GENOTYPE] += 1
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
            add_expansion_hunter_motif_composition(summary, motif_composition)

    consensus_sequences = variant_record.get("ConsensusSequences")
    if consensus_sequences:
        stats[WITH_CONSENSUS_SEQUENCES] += 1
        if settings.track_allele_sequences and alleles is not None \
                and not add_consensus_sequences(summary, alleles, consensus_sequences):
            stats[WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE] += 1
    return summary


def add_expansion_hunter_motif_composition(summary, motif_composition):
    """Adds one sample's locus-level MotifComposition record, translating its motif IDs to motif sequences.

    Args:
        summary (RepeatSummary): the repeat's summary.
        motif_composition (dict): the record, e.g. {"Motifs": {"1:GCG": {"count": 561, "reads": 64}},
            "MotifPairs": {"[1][1]": {"count": 398, "reads": 63}}}. An empty record (a locus with no reads) adds
            nothing.
    """
    if summary.motif_counts is None:
        summary.motif_counts, summary.motif_pair_counts = {}, {}
    motif_by_id = {}
    for key, counts in motif_composition.get("Motifs", {}).items():
        motif_id, motif = key.split(":", 1)
        motif_by_id[motif_id] = motif
        add_occurrences_and_reads(summary.motif_counts, motif, counts["count"], counts["reads"])
    for key, counts in motif_composition.get("MotifPairs", {}).items():
        match = MOTIF_PAIR_KEY_REGEX.fullmatch(key)
        if match is None:
            raise ValueError(f"Unexpected key in an ExpansionHunter MotifComposition.MotifPairs record: {key}")
        pair = (motif_by_id[match.group(1)], motif_by_id[match.group(2)])
        add_occurrences_and_reads(summary.motif_pair_counts, pair, counts["count"], counts["reads"])


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
    if not target.locus_structure:
        target.locus_structure = source.locus_structure
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
    if source.motif_counts is not None:
        if target.motif_counts is None:
            target.motif_counts, target.motif_pair_counts = {}, {}
        for motif, (occurrences, reads) in source.motif_counts.items():
            add_occurrences_and_reads(target.motif_counts, motif, occurrences, reads)
        for pair, (occurrences, reads) in source.motif_pair_counts.items():
            add_occurrences_and_reads(target.motif_pair_counts, pair, occurrences, reads)
    if source.allele_sequence_counts is not None:
        if target.allele_sequence_counts is None:
            target.allele_sequence_counts = {}
        for sequence, count in source.allele_sequence_counts.items():
            add_count(target.allele_sequence_counts, sequence, count)


def format_histogram(histogram):
    return ",".join(f"{allele}x:{count}" for allele, count in sorted(histogram.items()))


def parse_histogram(histogram_string):
    """Parses "15x:120,16x:3" into {15: 120, 16: 3}."""
    histogram = {}
    for entry in filter(None, histogram_string.split(",")):
        allele, count = entry.split("x:")
        histogram[int(allele)] = int(count)
    return histogram


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
    """Lists the samples with the largest alleles, the same way combine_single_sample_LPS_to_allele_histograms.py does.

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


def format_motif_composition(motif_counts, motif_pair_counts):
    """Formats summed motif and motif pair counts as the MotifComposition and MotifPairs values.

    Motifs are numbered by total occurrences, highest first, with ties broken by motif sequence, and the motif pairs
    use the same numbers.

    Returns:
        tuple: (MotifComposition value, MotifPairs value), e.g. ("1:GCG:561:64,2:GCA:82:50",
            "[1][1]:398:63,[1][2]:82:50"), or two empty strings if there are no counts.
    """
    if not motif_counts:
        return "", ""
    ordered_motifs = sorted(motif_counts.items(), key=lambda item: (-item[1][0], item[0]))
    motif_ids = {motif: index + 1 for index, (motif, _) in enumerate(ordered_motifs)}
    motifs_value = ",".join(
        f"{motif_ids[motif]}:{motif}:{occurrences}:{reads}" for motif, (occurrences, reads) in ordered_motifs)
    ordered_pairs = sorted(motif_pair_counts.items(), key=lambda item: (motif_ids[item[0][0]], motif_ids[item[0][1]]))
    pairs_value = ",".join(
        f"[{motif_ids[motif]}][{motif_ids[next_motif]}]:{occurrences}:{reads}"
        for (motif, next_motif), (occurrences, reads) in ordered_pairs)
    return motifs_value, pairs_value


def parse_motif_composition_columns(motifs_value, pairs_value):
    """Parses MotifComposition and MotifPairs values back into counts keyed by motif sequence.

    Returns:
        tuple: (motif -> [occurrences, reads], (motif, next motif) -> [occurrences, reads])
    """
    motif_counts, motif_pair_counts, motif_by_id = {}, {}, {}
    for entry in filter(None, motifs_value.split(",")):
        motif_id, motif, occurrences, reads = entry.split(":")
        motif_by_id[motif_id] = motif
        add_occurrences_and_reads(motif_counts, motif, int(occurrences), int(reads))
    for entry in filter(None, pairs_value.split(",")):
        key, occurrences, reads = entry.split(":")
        match = MOTIF_PAIR_KEY_REGEX.fullmatch(key)
        if match is None:
            raise ValueError(f"Unexpected entry in a TSV's MotifPairs column: {entry}")
        pair = (motif_by_id[match.group(1)], motif_by_id[match.group(2)])
        add_occurrences_and_reads(motif_pair_counts, pair, int(occurrences), int(reads))
    return motif_counts, motif_pair_counts


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


def get_locus_structure(locus_id, motif, num_repeats_in_locus, catalog_locus_structures):
    """Returns the LocusStructure to report: the catalog's if known, "(<motif>)*" for a single repeat, else ""."""
    if locus_id in catalog_locus_structures:
        return catalog_locus_structures[locus_id]
    return f"({motif})*" if num_repeats_in_locus == 1 else ""


def read_records_before_locus_results(json_path):
    """Returns {name: record} for the top-level records of an ExpansionHunter JSON that come before LocusResults (all
    of them if it has no LocusResults)."""
    records = {}
    with open_binary(json_path) as f:
        events = ijson.parse(f, use_float=True)
        for prefix, event, name in events:
            if prefix != "" or event != "map_key":
                continue
            if name == "LocusResults":
                break
            builder = ijson.ObjectBuilder()
            depth = 0
            for _, value_event, value in events:
                builder.event(value_event, value)
                if value_event in ("start_map", "start_array"):
                    depth += 1
                elif value_event in ("end_map", "end_array"):
                    depth -= 1
                if depth == 0:
                    break
            records[name] = builder.value
    return records


def read_expansion_hunter_json_header(json_path):
    """Returns (SampleId, Sex, loci_sorted_by_locus_id, command_line_has_motif_composition) for an ExpansionHunter
    output JSON, where command_line_has_motif_composition says whether RunInfo.CommandLine has
    --output-motif-composition.

    ExpansionHunter writes SampleParameters, RunInfo, LocusResults and RunStats, in that order, in every mode. Its
    seeking and streaming modes write the loci sorted by LocusId, and its optimized-streaming mode nearly in genomic
    order. Older versions wrote RunInfo after the loci (with RunStats' fields in it): optimized-streaming output as
    SampleParameters, LocusResults, RunInfo, and seeking and streaming output with the keys sorted, as LocusResults,
    RunInfo, SampleParameters. Only optimized-streaming output starts with SampleParameters. Getting SampleParameters
    from older seeking or streaming output means reading the whole file, which passes RunInfo on the way; for older
    optimized-streaming output, command_line_has_motif_composition is returned as None (RunInfo is read from the end
    of the file while the loci are streamed).

    SampleId falls back to the file name, and Sex to "" if the JSON doesn't have it.
    """
    records_before_loci = read_records_before_locus_results(json_path)
    sample_parameters = records_before_loci.get("SampleParameters")
    run_info = records_before_loci.get("RunInfo")
    if run_info is not None:
        loci_sorted_by_locus_id = run_info.get("AnalysisMode") != "optimized-streaming"
    else:
        loci_sorted_by_locus_id = "SampleParameters" not in records_before_loci
    if sample_parameters is None:
        sample_parameters = {}
        with open_binary(json_path) as raw_file:
            f = RunInfoKeepingReader(raw_file)
            for sample_parameters in ijson.items(f, "SampleParameters"):
                break
        if run_info is None:
            run_info = parse_run_info(f.tail)
    command_line_has_motif_composition = None if run_info is None else run_info_has_motif_composition(run_info)
    sample_id = sample_parameters.get("SampleId") or os.path.basename(strip_suffix(json_path, JSON_SUFFIXES))
    return sample_id, sample_parameters.get("Sex", ""), loci_sorted_by_locus_id, command_line_has_motif_composition


class RunInfoKeepingReader:
    """Wraps a binary ExpansionHunter JSON and keeps the bytes read from it since the last RUN_INFO_KEY.

    Until the key is seen, only its length minus one bytes are kept, so a key split between two reads is still found.
    If more than MAX_RUN_INFO_BYTES follow the key, they are dropped, and tail no longer holds RunInfo.
    """

    def __init__(self, f):
        self.f = f
        self.tail = bytearray()
        self.tail_starts_with_run_info_key = False

    def read(self, size=-1):
        data = self.f.read(size)
        search_start = max(0, len(self.tail) - len(RUN_INFO_KEY) + 1)
        self.tail += data
        index = self.tail.rfind(RUN_INFO_KEY, search_start)
        if index >= 0:
            del self.tail[:index]
            self.tail_starts_with_run_info_key = True
        elif len(self.tail) > MAX_RUN_INFO_BYTES:
            self.tail_starts_with_run_info_key = False
        if not self.tail_starts_with_run_info_key:
            del self.tail[:-(len(RUN_INFO_KEY) - 1)]
        return data


def parse_run_info(json_tail):
    """Returns the RunInfo record found in the end of an ExpansionHunter JSON, or None if it has none."""
    index = json_tail.rfind(RUN_INFO_KEY)
    if index < 0:
        return None
    text = json_tail[json_tail.index(b":", index) + 1:].decode("utf-8", "replace").lstrip()
    try:
        return json.JSONDecoder().raw_decode(text)[0]
    except ValueError:
        return None


def run_info_has_motif_composition(run_info):
    """Returns whether a RunInfo record (or None) has --output-motif-composition in its CommandLine."""
    return run_info is not None and "--output-motif-composition" in run_info.get("CommandLine", "")


def json_command_line_has_motif_composition(input_spec):
    """Returns whether an ExpansionHunter JSON's RunInfo.CommandLine has --output-motif-composition.

    In JSONs from older ExpansionHunter versions RunInfo is one of the last records, so the whole file is
    decompressed, but it isn't parsed.
    """
    if input_spec.command_line_has_motif_composition is not None:
        return input_spec.command_line_has_motif_composition
    with open_binary(input_spec.path) as raw_file:
        f = RunInfoKeepingReader(raw_file)
        while f.read(1 << 20):
            pass
    return run_info_has_motif_composition(parse_run_info(f.tail))


def iterate_expansion_hunter_records(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each repeat record of an ExpansionHunter output JSON, in file order.

    The JSON is streamed, so the file is never fully loaded. stats records whether its RunInfo.CommandLine has
    --output-motif-composition (COMMAND_LINE_HAS_MOTIF_COMPOSITION), which in JSONs from older ExpansionHunter versions
    is only known once the whole file has been read, and whether it was read to its end (READ_TO_THE_END).
    """
    is_male = input_spec.sex.lower() == "male"
    command_line_has_motif_composition = input_spec.command_line_has_motif_composition
    if command_line_has_motif_composition:
        stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] = 1
    with open_binary(input_spec.path) as raw_file:
        # In JSONs from older ExpansionHunter versions RunInfo comes after the loci, so it is picked out while they
        # are read. (Wrapping a JSON whose RunInfo comes first would keep everything after it.)
        f = raw_file if command_line_has_motif_composition is not None else RunInfoKeepingReader(raw_file)
        for locus_id, locus_record in ijson.kvitems(f, "LocusResults", use_float=True):
            repeat_records = [(variant_id, variant_record)
                              for variant_id, variant_record in locus_record.get("Variants", {}).items()
                              if variant_record.get("VariantType") == "Repeat"]
            for variant_id, variant_record in repeat_records:
                locus_structure = get_locus_structure(locus_id, variant_record["RepeatUnit"], len(repeat_records),
                                                      settings.catalog_locus_structures)
                summary = summarize_expansion_hunter_record(variant_record, input_spec.sample_id, is_male,
                                                            locus_structure, settings, stats)
                yield get_genomic_sort_key(locus_id, variant_id, summary.reference_region), summary
        # kvitems only finishes once ijson has read the whole file, so an older JSON's tail holds RunInfo (its last
        # record in optimized-streaming output; in seeking and streaming output only SampleParameters comes after it).
        if command_line_has_motif_composition is None and run_info_has_motif_composition(parse_run_info(f.tail)):
            stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] = 1
    stats[READ_TO_THE_END] = 1


def iterate_expansion_hunter_json_in_genomic_order(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each repeat record of an optimized-streaming JSON, in genomic order.

    Records are put back in order through a buffer that holds the repeats of the last REORDER_WINDOW_BP on the
    chromosome.

    Raises:
        ValueError: if a repeat comes more than REORDER_WINDOW_BP after repeats that start later, or the JSON's
            chromosomes are not in the order get_genomic_sort_key sorts them in.
    """
    reorder_buffer = []
    counter = 0
    current_chrom_key = None
    largest_start = 0
    last_yielded_key = None
    for key, summary in iterate_expansion_hunter_records(input_spec, settings, stats):
        if key[0] != current_chrom_key:
            if current_chrom_key is not None and key[0] < current_chrom_key:
                raise ValueError(
                    f"{input_spec.path} lists {get_chromosome(summary.reference_region)} after a chromosome that "
                    f"sorts after it. Optimized-streaming output lists chromosomes in the order of the BAM or CRAM "
                    f"header, and this script needs 1-22, X, Y, M, then others by name.")
            while reorder_buffer:
                last_yielded_key, _, buffered_summary = heapq.heappop(reorder_buffer)
                yield last_yielded_key, buffered_summary
            current_chrom_key = key[0]
            largest_start = 0
        if last_yielded_key is not None and key < last_yielded_key:
            raise ValueError(
                f"{input_spec.path} lists {key[4]} ({summary.reference_region}) more than {REORDER_WINDOW_BP:,d} bp "
                f"after repeats that start later, which optimized-streaming output isn't expected to do.")
        heapq.heappush(reorder_buffer, (key, counter, summary))
        counter += 1
        largest_start = max(largest_start, key[1])
        while reorder_buffer and reorder_buffer[0][0][1] < largest_start - REORDER_WINDOW_BP:
            last_yielded_key, _, buffered_summary = heapq.heappop(reorder_buffer)
            yield last_yielded_key, buffered_summary
    while reorder_buffer:
        key, _, buffered_summary = heapq.heappop(reorder_buffer)
        yield key, buffered_summary


def write_sorted_tsv_from_expansion_hunter_json(input_spec, output_tsv_path, settings, stats, temp_dir):
    """Sorts a seeking or streaming mode JSON (loci sorted by LocusId) into a single-sample TSV in genomic order.

    Each repeat's row is written to a temporary file for its chromosome, and then each chromosome's rows are sorted
    in memory and written out, so memory use is about one chromosome's rows.
    """
    with tempfile.TemporaryDirectory(dir=temp_dir) as bucket_dir:
        bucket_files = {}
        try:
            for key, summary in iterate_expansion_hunter_records(input_spec, settings, stats):
                if key[0] not in bucket_files:
                    bucket_files[key[0]] = open(os.path.join(bucket_dir, f"chromosome_{len(bucket_files)}.tsv"), "w")
                bucket_files[key[0]].write(format_row(key, summary, settings) + "\n")
        finally:
            for bucket_file in bucket_files.values():
                bucket_file.close()
        with gzip.open(output_tsv_path, "wt", compresslevel=TEMPORARY_FILE_COMPRESSION_LEVEL) as f:
            f.write("\t".join(get_output_columns(settings)) + "\n")
            for chrom_key in sorted(bucket_files):
                with open(bucket_files[chrom_key].name) as bucket_file:
                    lines = bucket_file.readlines()
                lines.sort(key=get_sort_key_of_row)
                f.writelines(lines)


def read_tsv_header(tsv_path):
    with open_text(tsv_path) as f:
        return f.readline().rstrip("\n").split("\t")


def tsv_has_motif_composition_columns(tsv_path):
    header = read_tsv_header(tsv_path)
    return all(column in header for column in MOTIF_COMPOSITION_COLUMNS)


def iterate_multisample_tsv(input_spec, settings, stats):
    """Yields (sort key, RepeatSummary) for each row of a TSV written by this script.

    Raises:
        ValueError: if its rows are not sorted.
    """
    header = read_tsv_header(input_spec.path)
    column_index = {column: index for index, column in enumerate(header)}
    previous_key = None
    with open_text(input_spec.path) as f:
        next(f)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            row = {column: fields[index] for column, index in column_index.items()}
            stats[TSV_ROWS] += 1
            locus_structure = settings.catalog_locus_structures.get(row["LocusId"]) or row["LocusStructure"]
            summary = RepeatSummary(row["Motif"], row["ReferenceRegion"], locus_structure,
                                    settings.track_outlier_sample_ids)
            summary.sample_count = int(row["SampleCount"])
            summary.called_allele_count = int(row["CalledAlleleCount"])
            summary.all_allele_histogram = parse_histogram(row["AllAlleleHistogram"])
            summary.short_allele_histogram = parse_histogram(row["ShortAlleleHistogram"])
            summary.hemizygous_allele_histogram = parse_histogram(row["HemizygousAlleleHistogram"])
            summary.biallelic_histogram = parse_biallelic_histogram(row["BiallelicHistogram"])
            if settings.track_outlier_sample_ids:
                summary.all_allele_outlier_sample_ids = parse_outlier_sample_ids(row["OutlierSampleIds_AllAlleles"])
                summary.short_allele_outlier_sample_ids = parse_outlier_sample_ids(
                    row["OutlierSampleIds_ShortAlleles"])
                summary.hemizygous_allele_outlier_sample_ids = parse_outlier_sample_ids(
                    row["OutlierSampleIds_HemizygousAlleles"])
            if settings.track_motif_composition:
                summary.motif_counts, summary.motif_pair_counts = parse_motif_composition_columns(
                    row["MotifComposition"], row["MotifPairs"])
            if settings.track_allele_sequences:
                summary.allele_sequence_counts = parse_allele_sequences(row[ALLELE_SEQUENCES_COLUMN])
            key = get_genomic_sort_key(row["LocusId"], row["VariantId"], row["ReferenceRegion"])
            if previous_key is not None and key <= previous_key:
                raise ValueError(f"{input_spec.path} is not sorted: {row['LocusId']} {row['VariantId']} comes after "
                                 f"{previous_key[3]} {previous_key[4]}")
            previous_key = key
            yield key, summary


def get_output_columns(settings):
    columns = list(BASE_COLUMNS)
    if settings.track_outlier_sample_ids:
        columns += [MAX_OUTLIER_SAMPLE_IDS_COLUMN] + OUTLIER_COLUMNS
    if settings.track_motif_composition:
        columns += MOTIF_COMPOSITION_COLUMNS
    if settings.track_allele_sequences:
        columns.append(ALLELE_SEQUENCES_COLUMN)
    return columns


def format_row(key, summary, settings):
    """Returns the output TSV line (without the newline) for one repeat."""
    locus_id, variant_id = key[3], key[4]
    row = [
        locus_id, variant_id, summary.reference_region, summary.motif, summary.locus_structure,
        str(summary.sample_count), str(summary.called_allele_count),
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
        row += format_motif_composition(summary.motif_counts, summary.motif_pair_counts)
    if settings.track_allele_sequences:
        row.append(format_allele_sequences(summary.allele_sequence_counts))
    return "\t".join(row)


def merge_inputs(input_specs, output_tsv_path, settings, temp_dir, show_progress_bar=False,
                 compression_level=OUTPUT_COMPRESSION_LEVEL):
    """Merges the inputs one repeat at a time into a TSV with the columns get_output_columns(settings) gives.

    Seeking and streaming mode JSONs are first sorted into temporary TSVs in temp_dir. With settings.max_rows_to_write,
    the merge stops once that many rows are written.

    Args:
        input_specs (list): InputSpec of each input.
        output_tsv_path (str): where to write the merged TSV (gzipped).
        settings (MergeSettings): the run's settings.
        temp_dir (str): directory for temporary files.
        show_progress_bar (bool): whether to show a progress bar of input records merged.
        compression_level (int): gzip compression level of the output.

    Returns:
        tuple: (list of per-input collections.Counter stats, in input_specs order; number of rows written;
            number of rows whose LocusStructure is unknown)
    """
    stats_by_input = [collections.Counter() for _ in input_specs]
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
                iterators.append(iterate_expansion_hunter_json_in_genomic_order(input_spec, settings, stats))
        merged = heapq.merge(*iterators, key=operator.itemgetter(0))
        if show_progress_bar:
            merged = tqdm.tqdm(merged, unit=" records", unit_scale=True, desc="Merging")

        num_rows = 0
        num_rows_without_locus_structure = 0
        with gzip.open(output_tsv_path, "wt", compresslevel=compression_level) as f:
            f.write("\t".join(get_output_columns(settings)) + "\n")
            for key, records in itertools.groupby(merged, key=operator.itemgetter(0)):
                _, summary = next(records)
                for _, other_summary in records:
                    merge_repeat_summaries(summary, other_summary, settings.n_outlier_sample_ids,
                                           f"{key[4]} in locus {key[3]}")
                f.write(format_row(key, summary, settings) + "\n")
                num_rows += 1
                num_rows_without_locus_structure += not summary.locus_structure
                if num_rows == settings.max_rows_to_write:
                    break
    finally:
        # Closing the iterators closes their files, which a merge stopped by -n leaves part-read.
        for iterator in iterators:
            iterator.close()
        remove_files(sorted_tsv_paths)
    return stats_by_input, num_rows, num_rows_without_locus_structure


def merge_input_group(task):
    """Merges one group of inputs into a temporary TSV; the worker-process entry point for --threads.

    Args:
        task (tuple): (group index, list of InputSpec, output TSV path, MergeSettings, temp dir)

    Returns:
        tuple: (group index, list of per-input stats)
    """
    group_index, input_specs, output_tsv_path, settings, temp_dir = task
    print(f"Merging {len(input_specs):,d} input(s) into {output_tsv_path}", flush=True)
    stats_by_input, _, _ = merge_inputs(input_specs, output_tsv_path, settings, temp_dir,
                                        compression_level=TEMPORARY_FILE_COMPRESSION_LEVEL)
    return group_index, stats_by_input


def json_has_motif_composition_in_first_records(json_path, max_records_to_check):
    """Returns whether an ExpansionHunter JSON has MotifComposition in its first max_records_to_check repeat records.

    False doesn't rule out motif composition further on (ExpansionHunter's loci-with-non-ref-motifs and
    loci-with-known-motifs modes put it on few loci); main then checks a JSON's command line too.
    """
    num_checked = 0
    with open_binary(json_path) as f:
        for _, locus_record in ijson.kvitems(f, "LocusResults", use_float=True):
            for variant_record in locus_record.get("Variants", {}).values():
                if variant_record.get("VariantType") != "Repeat":
                    continue
                if "MotifComposition" in variant_record:
                    return True
                num_checked += 1
                if num_checked >= max_records_to_check:
                    return False
    return False


def raise_open_file_limit():
    """Raises this process's soft limit on open files as far as the system allows (worker processes inherit it)."""
    soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)
    # macOS rejects RLIM_INFINITY and values above kern.maxfilesperproc, so fall back to smaller targets.
    for target in [hard_limit, 65_536, 10_240]:
        if target == resource.RLIM_INFINITY or target <= soft_limit:
            continue
        if hard_limit != resource.RLIM_INFINITY and target > hard_limit:
            continue
        try:
            resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard_limit))
            return
        except (ValueError, OSError):
            continue


def get_max_inputs_per_merge():
    """Returns how many inputs one merge may read at once: MAX_INPUTS_PER_MERGE, or fewer if the open-file limit is
    lower."""
    soft_limit = resource.getrlimit(resource.RLIMIT_NOFILE)[0]
    return max(2, min(MAX_INPUTS_PER_MERGE, soft_limit - OPEN_FILES_NEEDED_BESIDES_INPUTS))


def split_into_groups(input_specs, num_groups, file_size_by_path):
    """Splits inputs into num_groups groups of similar total file size, with at most ceil(N / num_groups) each."""
    max_group_size = math.ceil(len(input_specs) / num_groups)
    groups = [[] for _ in range(num_groups)]
    group_sizes = [0] * num_groups
    for input_spec in sorted(input_specs, key=lambda spec: (-file_size_by_path[spec.path], spec.path)):
        open_groups = [index for index in range(num_groups) if len(groups[index]) < max_group_size]
        smallest_group = min(open_groups, key=lambda index: (group_sizes[index], index))
        groups[smallest_group].append(input_spec)
        group_sizes[smallest_group] += file_size_by_path[input_spec.path]
    return groups


def merge_groups(grouped_specs, temp_dir, settings, threads, show_progress_bar):
    """Merges each group of inputs into a temporary TSV, in up to `threads` worker processes.

    A single group is merged in this process straight into its final form (output compression level).

    Returns:
        tuple: (list of the groups' TSV paths, list of each group's per-input stats, number of rows written and
            number of rows without LocusStructure for a single group, or None and None)
    """
    if len(grouped_specs) == 1:
        merged_path = os.path.join(temp_dir, "merged.tsv.gz")
        print(f"Merging the inputs into {merged_path}")
        stats_by_input, num_rows, num_rows_without_locus_structure = merge_inputs(
            grouped_specs[0], merged_path, settings, temp_dir, show_progress_bar)
        return [merged_path], [stats_by_input], num_rows, num_rows_without_locus_structure

    group_paths = [os.path.join(temp_dir, f"group_{index + 1}.tsv.gz") for index in range(len(grouped_specs))]
    tasks = [(index, group, group_path, settings, temp_dir)
             for index, (group, group_path) in enumerate(zip(grouped_specs, group_paths))]
    num_workers = min(threads, len(tasks))
    print(f"Merging the inputs in {len(tasks):,d} groups, {num_workers:,d} at a time")
    grouped_stats = [None] * len(tasks)
    if num_workers == 1:
        results = map(merge_input_group, tasks)
        if show_progress_bar:
            results = tqdm.tqdm(results, total=len(tasks), unit=" groups")
        for group_index, stats_by_input in results:
            grouped_stats[group_index] = stats_by_input
    else:
        # Leaving the with block (normally, or because a group's merge raised) terminates the other workers.
        with multiprocessing.Pool(processes=num_workers) as pool:
            results = pool.imap_unordered(merge_input_group, tasks)
            if show_progress_bar:
                results = tqdm.tqdm(results, total=len(tasks), unit=" groups")
            for group_index, stats_by_input in results:
                grouped_stats[group_index] = stats_by_input
    return group_paths, grouped_stats, None, None


def merge_tsvs(tsv_paths, output_tsv_path, settings, temp_dir, show_progress_bar):
    """Merges TSVs written by this script into one, in rounds if there are more than one merge can read at once.

    The input TSVs are deleted once merged.

    Returns:
        tuple: (number of rows written, number of rows whose LocusStructure is unknown)
    """
    max_inputs = get_max_inputs_per_merge()
    round_number = 0
    while len(tsv_paths) > max_inputs:
        round_number += 1
        next_round_paths = []
        for start in range(0, len(tsv_paths), max_inputs):
            next_round_path = os.path.join(temp_dir, f"round_{round_number}_part_{len(next_round_paths) + 1}.tsv.gz")
            merge_inputs([tsv_input_spec(path) for path in tsv_paths[start:start + max_inputs]], next_round_path,
                         settings, temp_dir, compression_level=TEMPORARY_FILE_COMPRESSION_LEVEL)
            next_round_paths.append(next_round_path)
        remove_files(tsv_paths)
        tsv_paths = next_round_paths
    _, num_rows, num_rows_without_locus_structure = merge_inputs(
        [tsv_input_spec(path) for path in tsv_paths], output_tsv_path, settings, temp_dir, show_progress_bar)
    remove_files(tsv_paths)
    return num_rows, num_rows_without_locus_structure


def read_locus_structures_from_catalog(catalog_path):
    """Returns LocusId -> LocusStructure for the catalog loci whose structure isn't a single "(<motif>)*" repeat.

    Single-repeat structures are not kept, since "(<Motif>)*" can be rebuilt from the ExpansionHunter output.
    """
    locus_structures = {}
    with open_binary(catalog_path) as f:
        for record in ijson.items(f, "item"):
            locus_structure = record.get("LocusStructure", "")
            if not SINGLE_REPEAT_LOCUS_STRUCTURE_REGEX.fullmatch(locus_structure):
                locus_structures[record["LocusId"]] = locus_structure
    return locus_structures


def parse_args(argv=None):
    parser = configargparse.ArgumentParser(
        description=__doc__, formatter_class=configargparse.RawDescriptionHelpFormatter)
    parser.add_argument("-c", "--config-file", is_config_file=True,
                        help="Config file that sets any of the options below (e.g. 'show-progress-bar = true')")
    parser.add_argument("-o", "--output-prefix",
                        help="Writes <prefix>.tsv.gz (a local path). Default: derived from "
                             "the input file name for one input, otherwise combined.<N>_inputs.multisample")
    parser.add_argument("-n", "--max-rows-to-write", type=int,
                        help="Write only the first N rows: the first N repeats in genomic order. Don't combine such a "
                             "TSV with other inputs later. With optimized-streaming JSONs from older ExpansionHunter "
                             "versions the motif composition columns can be missing or incomplete; a warning says "
                             "when")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of worker processes. Each merges a share of the inputs into a temporary TSV "
                             "next to the output, and those are then merged (default: %(default)s)")
    parser.add_argument("--variant-catalog",
                        help="ExpansionHunter catalog JSON (.json or .json.gz) to take the LocusStructure of loci "
                             "with several repeats, or with fixed sequence next to the repeat, from")
    parser.add_argument("--add-outlier-sample-id-columns", action="store_true",
                        help=f"Add the {MAX_OUTLIER_SAMPLE_IDS_COLUMN}, OutlierSampleIds_AllAlleles, "
                             f"OutlierSampleIds_ShortAlleles and OutlierSampleIds_HemizygousAlleles columns")
    parser.add_argument("--n-outlier-sample-ids", type=int, default=10,
                        help="Number of samples with the largest alleles to list in each outlier column (default: "
                             "%(default)s)")
    parser.add_argument("--skip-motif-composition-columns", action="store_true",
                        help="Leave out the MotifComposition and MotifPairs columns even if the inputs "
                             "have motif composition")
    parser.add_argument("--add-consensus-sequence-columns", action="store_true",
                        help="Add the AlleleSequences column, from the inputs' ConsensusSequences")
    parser.add_argument("--allow-overlapping-samples", action="store_true",
                        help="Allow a SampleId to appear in more than one JSON input. Only use this when those "
                             "inputs cover different loci (such as per-chromosome runs of the same samples): a sample "
                             "counted twice at a locus has its alleles counted twice.")
    parser.add_argument("--show-progress-bar", action="store_true", help="Show progress bars")
    parser.add_argument("input_paths", nargs="+",
                        help="ExpansionHunter output JSON files (.json, .json.gz, .json.bgz) and/or TSVs written by "
                             "earlier runs of this script (.tsv, .tsv.gz, .tsv.bgz), local or gs:// paths")
    args = parser.parse_args(argv)

    for path in args.input_paths:
        if not path.endswith(JSON_SUFFIXES + TSV_SUFFIXES):
            parser.error(f"{path} must end with one of {', '.join(JSON_SUFFIXES + TSV_SUFFIXES)}")
        if not file_exists(path):
            parser.error(f"{path} not found")
    if len(set(args.input_paths)) < len(args.input_paths):
        parser.error("An input path is listed more than once")
    if args.variant_catalog and not file_exists(args.variant_catalog):
        parser.error(f"{args.variant_catalog} not found")
    if args.output_prefix and is_gcs_path(args.output_prefix):
        parser.error("--output-prefix must be a local path; copy the output to Google Cloud Storage afterwards")
    if args.max_rows_to_write is not None and args.max_rows_to_write < 1:
        parser.error("-n must be at least 1")
    if args.n_outlier_sample_ids < 2:
        parser.error("--n-outlier-sample-ids must be at least 2: the outlier walk stops at the first allele size that "
                     "this many samples have, so with 1 the outlier columns would always be empty")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if not args.output_prefix:
        if len(args.input_paths) == 1:
            # Next to a local input; in the current directory for a gs:// one.
            input_path = args.input_paths[0]
            args.output_prefix = strip_suffix(os.path.basename(input_path) if is_gcs_path(input_path) else input_path,
                                              JSON_SUFFIXES + TSV_SUFFIXES) + ".multisample"
        else:
            args.output_prefix = f"combined.{len(args.input_paths)}_inputs.multisample"
    return args


def read_outlier_limit_of_tsv(tsv_path, header):
    """Returns the MaxOutlierSampleIds of a TSV's first row, or None if it has no rows."""
    with open_text(tsv_path) as f:
        next(f)
        first_row = f.readline()
    if not first_row:
        return None
    return int(first_row.rstrip("\n").split("\t")[header.index(MAX_OUTLIER_SAMPLE_IDS_COLUMN)])


def check_tsv_input_columns(tsv_path, args):
    """Checks a TSV input has the columns of a TSV written by this script, the optional ones this run needs, and an
    outlier limit at least as large as this run's.

    Raises:
        ValueError: if its columns are not ones this script writes, a needed column is missing, or the TSV's outlier
            columns were made with a smaller --n-outlier-sample-ids.
    """
    header = read_tsv_header(tsv_path)
    headers_this_script_writes = [
        get_output_columns(MergeSettings(None, track_outlier_sample_ids, track_motif_composition, track_allele_sequences,
                                         None, None))
        for track_outlier_sample_ids in (False, True) for track_motif_composition in (False, True)
        for track_allele_sequences in (False, True)]
    if header not in headers_this_script_writes:
        raise ValueError(f"{tsv_path} doesn't have the columns of a TSV written by this script: {', '.join(header)}")
    missing_columns = []
    hint = ""
    if args.add_outlier_sample_id_columns:
        missing_outlier_columns = [column for column in [MAX_OUTLIER_SAMPLE_IDS_COLUMN] + OUTLIER_COLUMNS
                                   if column not in header]
        if missing_outlier_columns:
            missing_columns += missing_outlier_columns
            hint = " It was made without --add-outlier-sample-id-columns, so its outlier samples are unknown."
    if args.add_consensus_sequence_columns and ALLELE_SEQUENCES_COLUMN not in header:
        missing_columns.append(ALLELE_SEQUENCES_COLUMN)
        hint += " It was made without --add-consensus-sequence-columns."
    if missing_columns:
        raise ValueError(f"{tsv_path} is missing the column(s) {', '.join(missing_columns)}.{hint}")
    if args.add_outlier_sample_id_columns:
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
            a column this run needs.
    """
    input_specs = []
    json_path_by_sample_id = {}
    for path in args.input_paths:
        if path.endswith(JSON_SUFFIXES):
            sample_id, sex, loci_sorted_by_locus_id, command_line_has_motif_composition = (
                read_expansion_hunter_json_header(path))
            if sample_id in json_path_by_sample_id and not args.allow_overlapping_samples:
                raise ValueError(
                    f"{path} and {json_path_by_sample_id[sample_id]} have the same SampleId, {sample_id}. If the "
                    f"inputs cover different loci, such as per-chromosome runs of the same sample, pass "
                    f"--allow-overlapping-samples.")
            json_path_by_sample_id[sample_id] = path
            input_specs.append(InputSpec(
                path, "json", sample_id, sex, loci_sorted_by_locus_id, command_line_has_motif_composition))
        else:
            check_tsv_input_columns(path, args)
            input_specs.append(tsv_input_spec(path))
    return input_specs


def check_optional_columns(input_specs, stats_by_input, args):
    """Decides whether the output gets motif composition columns, and checks the optional data is in every input.

    Returns:
        bool: whether to write the motif composition columns.

    Raises:
        ValueError: if only some inputs have motif composition (without --skip-motif-composition-columns), or
            --add-consensus-sequence-columns was given and an input has no consensus sequences.
    """
    inputs_with_motif_composition, inputs_without_motif_composition, inputs_without_consensus_sequences = [], [], []
    inputs_with_unknown_motif_composition = []
    for input_spec, stats in zip(input_specs, stats_by_input):
        if input_spec.kind == "json":
            # A run with --output-motif-composition loci-with-non-ref-motifs (or loci-with-known-motifs) can have no
            # MotifComposition record at all, so the command line counts too.
            has_motif_composition = stats[WITH_MOTIF_COMPOSITION] > 0 or stats[COMMAND_LINE_HAS_MOTIF_COMPOSITION] > 0
            # ExpansionHunter writes ConsensusSequences only for called repeats, so a JSON without any call has none
            # to contribute rather than missing them. (TSVs without the column were rejected before the merge.)
            if stats[WITH_CONSENSUS_SEQUENCES] == 0 and stats[WITH_GENOTYPE] > 0:
                inputs_without_consensus_sequences.append(input_spec.path)
            # With -n, a JSON from an older ExpansionHunter version may not have been read as far as its RunInfo, so
            # one without MotifComposition records so far can't be told either way.
            if not has_motif_composition and input_spec.command_line_has_motif_composition is None \
                    and not stats[READ_TO_THE_END]:
                inputs_with_unknown_motif_composition.append(input_spec.path)
                continue
        else:
            has_motif_composition = tsv_has_motif_composition_columns(input_spec.path)
        (inputs_with_motif_composition if has_motif_composition else inputs_without_motif_composition).append(
            input_spec.path)

    if not args.skip_motif_composition_columns and inputs_with_motif_composition and inputs_without_motif_composition:
        raise ValueError(
            f"{len(inputs_with_motif_composition):,d} input(s) have motif composition (e.g. "
            f"{inputs_with_motif_composition[0]}) but {len(inputs_without_motif_composition):,d} don't (e.g. "
            f"{inputs_without_motif_composition[0]}), so the motif composition columns would only count some of the "
            f"samples. Run ExpansionHunter with --output-motif-composition on all of them, or pass "
            f"--skip-motif-composition-columns. (A JSON counts as having motif composition if any of its repeat "
            f"records has a MotifComposition field or its RunInfo.CommandLine has --output-motif-composition; a TSV, "
            f"if it has the motif composition columns.)")
    if args.add_consensus_sequence_columns and inputs_without_consensus_sequences:
        raise ValueError(
            f"--add-consensus-sequence-columns was given, but {len(inputs_without_consensus_sequences):,d} JSON "
            f"input(s) have called repeats without any ConsensusSequences (e.g. "
            f"{inputs_without_consensus_sequences[0]}). Was ExpansionHunter run with "
            f"--dont-output-consensus-sequences?")
    if not args.skip_motif_composition_columns and inputs_with_unknown_motif_composition:
        print(f"WARNING: {len(inputs_with_unknown_motif_composition):,d} JSON input(s) from older ExpansionHunter "
              f"versions (e.g. {inputs_with_unknown_motif_composition[0]}) have no MotifComposition in the rows read, "
              f"and -n stopped before their RunInfo, so whether they have motif composition couldn't be checked. The "
              f"motif composition columns may be missing, or count only some of the samples.")
    return not args.skip_motif_composition_columns and bool(inputs_with_motif_composition)


def print_input_stats(input_specs, stats_by_input):
    for input_spec, stats in zip(input_specs, stats_by_input):
        if input_spec.kind == "tsv":
            print(f"  {input_spec.path}: {stats[TSV_ROWS]:,d} rows")
            continue
        details = ", ".join(f"{count:,d} {description}" for description, count in stats.items()
                            if description not in (REPEAT_RECORDS, COMMAND_LINE_HAS_MOTIF_COMPOSITION, READ_TO_THE_END))
        print(f"  {input_spec.path} (sample {input_spec.sample_id}, sex {input_spec.sex or 'unknown'}): "
              f"{stats[REPEAT_RECORDS]:,d} repeat records" + (f", {details}" if details else ""))


def main(argv=None):
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)  # so status messages reach a log file as they happen
    raise_open_file_limit()
    args = parse_args(argv)
    json_count = sum(path.endswith(JSON_SUFFIXES) for path in args.input_paths)
    print(f"Combining {len(args.input_paths):,d} inputs: {json_count:,d} ExpansionHunter JSON file(s) and "
          f"{len(args.input_paths) - json_count:,d} TSV(s) from earlier runs")
    input_specs = collect_input_specs(args)
    samples_without_sex = sorted({input_spec.sample_id for input_spec in input_specs
                                  if input_spec.kind == "json" and input_spec.sex.lower() not in ("male", "female")})
    if samples_without_sex:
        print(f"WARNING: {len(samples_without_sex):,d} samples have no Sex in SampleParameters (e.g. "
              f"{samples_without_sex[0]}); they are counted as female, so only their chrY calls are hemizygous")

    catalog_locus_structures = {}
    if args.variant_catalog:
        print(f"Reading LocusStructures from {args.variant_catalog}")
        catalog_locus_structures = read_locus_structures_from_catalog(args.variant_catalog)
        print(f"  {len(catalog_locus_structures):,d} loci have a LocusStructure other than a single "
              f"\"(<motif>)*\" repeat")

    # Whether the output gets motif composition columns is only certain once every input has been read. A TSV's header
    # tells, and so do a JSON's first records in ExpansionHunter's all-loci mode. Its other modes can put
    # MotifComposition on few loci or none, so one JSON's RunInfo.CommandLine is checked too: inputs that differ in
    # motif composition are an error, so one JSON speaks for all. A JSON whose command line was read with its header
    # is used if there is one; otherwise the smallest, whose whole file has to be decompressed, which -n skips
    # (check_optional_columns then warns about the JSONs whose motif composition it can't tell). The merge is rerun if
    # MotifComposition turns up late in JSONs whose command line wasn't checked or abbreviates
    # --output-motif-composition. With -n, at most N records of each JSON are checked first, since the merge itself
    # reads about that many.
    # Each gs:// file's size takes a request to Google Cloud Storage, so it is looked up once.
    file_size_by_path = {input_spec.path: get_file_size(input_spec.path) for input_spec in input_specs}
    track_motif_composition = False
    if not args.skip_motif_composition_columns:
        tsv_paths = [input_spec.path for input_spec in input_specs if input_spec.kind == "tsv"]
        tsvs_without_motif_composition = [path for path in tsv_paths if not tsv_has_motif_composition_columns(path)]
        json_specs = [input_spec for input_spec in input_specs if input_spec.kind == "json"]
        json_spec_to_check_command_line_of = min(json_specs, default=None, key=lambda input_spec: (
            input_spec.command_line_has_motif_composition is None, file_size_by_path[input_spec.path]))
        check_command_line = json_spec_to_check_command_line_of is not None and (
            json_spec_to_check_command_line_of.command_line_has_motif_composition is not None
            or args.max_rows_to_write is None)
        max_records_to_check = min(MOTIF_COMPOSITION_CHECK_RECORDS, args.max_rows_to_write or math.inf)
        track_motif_composition = (
            len(tsvs_without_motif_composition) < len(tsv_paths)
            or any(json_has_motif_composition_in_first_records(input_spec.path, max_records_to_check)
                   for input_spec in json_specs)
            or (check_command_line and json_command_line_has_motif_composition(json_spec_to_check_command_line_of)))
        if track_motif_composition and tsvs_without_motif_composition:
            raise ValueError(
                f"Some inputs have motif composition but {tsvs_without_motif_composition[0]} has no motif composition "
                f"columns, so those columns would only count some of the samples. Pass "
                f"--skip-motif-composition-columns to leave them out.")
    settings = MergeSettings(
        n_outlier_sample_ids=args.n_outlier_sample_ids,
        track_outlier_sample_ids=args.add_outlier_sample_id_columns,
        track_motif_composition=track_motif_composition,
        track_allele_sequences=args.add_consensus_sequence_columns,
        catalog_locus_structures=catalog_locus_structures,
        max_rows_to_write=args.max_rows_to_write,
    )
    if args.max_rows_to_write is not None:
        print(f"Writing only the first {args.max_rows_to_write:,d} rows (-n), so the input counts below cover only the "
              f"records read for them")
    output_tsv_path = f"{args.output_prefix}.tsv.gz"
    num_groups = max(min(args.threads, len(input_specs)), math.ceil(len(input_specs) / get_max_inputs_per_merge()))
    grouped_specs = split_into_groups(input_specs, num_groups, file_size_by_path)
    with tempfile.TemporaryDirectory(prefix=os.path.basename(args.output_prefix) + ".tmp_",
                                     dir=os.path.dirname(os.path.abspath(output_tsv_path))) as temp_dir:
        while True:
            group_paths, grouped_stats, num_rows, num_rows_without_locus_structure = merge_groups(
                grouped_specs, temp_dir, settings, args.threads, args.show_progress_bar)
            ordered_specs = [spec for group in grouped_specs for spec in group]
            ordered_stats = [stats for group_stats in grouped_stats for stats in group_stats]
            write_motif_composition_columns = check_optional_columns(ordered_specs, ordered_stats, args)
            if not write_motif_composition_columns or settings.track_motif_composition:
                break
            print("Motif composition only appears after the repeat records checked before the merge, so merging again "
                  "with the motif composition columns")
            remove_files(group_paths)
            settings = settings._replace(track_motif_composition=True)

        print("Inputs:")
        print_input_stats(ordered_specs, ordered_stats)
        for description in [WITH_ALLELE_COUNT_NOT_MATCHING_SEX, WITH_CONSENSUS_SEQUENCES_NOT_MATCHING_GENOTYPE]:
            total = sum(stats[description] for stats in ordered_stats)
            if total:
                print(f"WARNING: {total:,d} genotypes {description}")

        if len(group_paths) == 1:
            merged_path = group_paths[0]
        else:
            merged_path = os.path.join(temp_dir, "merged.tsv.gz")
            print(f"Merging the {len(group_paths):,d} groups' TSVs into {merged_path}")
            num_rows, num_rows_without_locus_structure = merge_tsvs(
                group_paths, merged_path, settings, temp_dir, args.show_progress_bar)
        os.replace(merged_path, output_tsv_path)

    if num_rows_without_locus_structure:
        print(f"WARNING: {num_rows_without_locus_structure:,d} rows are repeats of loci with several repeats, whose "
              f"LocusStructure can't be rebuilt from the ExpansionHunter output; pass --variant-catalog to fill it in")
    added_columns = [name for name, added in [
        ("outlier sample ID", settings.track_outlier_sample_ids),
        ("motif composition", settings.track_motif_composition),
        (ALLELE_SEQUENCES_COLUMN, settings.track_allele_sequences)] if added]
    print(f"Wrote {num_rows:,d} rows to {output_tsv_path}"
          + (f" (with {', '.join(added_columns)} columns)" if added_columns else ""))


if __name__ == "__main__":
    main()
