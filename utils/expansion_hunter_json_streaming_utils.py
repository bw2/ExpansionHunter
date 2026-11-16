"""Code shared by the scripts that combine ExpansionHunter output JSON files from many samples
(generate_multisample_tsv.py and generate_multisample_qc_stats.py).

It covers:
    - reading local or gs:// files, gzipped or not, as streams
    - reading an ExpansionHunter JSON's SampleId, Sex and analysis mode from the records before its loci
    - streaming a JSON's repeat records and putting optimized-streaming output back into genomic order
    - splitting many inputs into shards that worker processes merge, and merging their outputs in rounds that stay
      under the open-file limit
"""

import collections
import concurrent.futures
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
import pickle
import resource
import sys
import tempfile

import ijson
import tqdm

JSON_SUFFIXES = (".json", ".json.gz", ".json.bgz")
# gzip's default level 9 is several times slower than these for little gain on the output tables
OUTPUT_COMPRESSION_LEVEL = 6
TEMPORARY_FILE_COMPRESSION_LEVEL = 1
# Checking that a gs:// file exists, or reading its header, takes a request or a few to Google Cloud Storage, so many
# files are checked in parallel threads
THREADS_FOR_GCS_REQUESTS = 32

# How far (in bp) a repeat may come after repeats that start later on the same chromosome. ExpansionHunter's
# optimized-streaming mode writes a locus once the reads have passed it, so a locus can come up to about a flank
# length (1,000 bp) after loci that start later; on a full TRExplorer v2.1 catalog the largest such distance is 955 bp.
# The buffer holds every repeat in the window, so this is about 10x the largest distance seen and no more.
REORDER_WINDOW_BP = 10_000
# The chromosome keys of get_genomic_sort_key in the order it sorts them, and each one's successor (None, before the
# first chromosome, maps to chr1). A JSON whose chromosomes follow each other this way is put back in order as it is
# streamed; one in another order, such as a BAM header sorted by name (chr1, chr10, ..., chr19, chr2, ...), is sorted.
STANDARD_CHROMOSOME_ORDER = [(0, number, "") for number in range(1, 23)] + [(1, 0, ""), (1, 1, ""), (1, 2, "")]
NEXT_STANDARD_CHROMOSOME_KEY = dict(zip([None] + STANDARD_CHROMOSOME_ORDER, STANDARD_CHROMOSOME_ORDER))
# ExpansionHunter genotypes chrY with one allele in every sample and chrX with one allele in males, except at loci its
# catalog marks "IsInsidePseudoautosomalRegion": true (determineExpectedAlleleCount in ehunter/core/LocusStats.cpp).
# Each locus's AlleleCount in the output JSON says which it did, so generate_multisample_tsv.py uses these names only
# when AlleleCount is missing. generate_multisample_qc_stats.py uses them to tell autosomes from sex chromosomes and
# to leave female samples' chrY repeats out of the call rates.
CHR_X_NAMES = {"chrX", "X"}
CHR_Y_NAMES = {"chrY", "Y"}
# Open files a merge needs besides its inputs (output, bucket files of a JSON being sorted, Python's own).
OPEN_FILES_NEEDED_BESIDES_INPUTS = 64
# Each input a merge reads costs about 1 MB of buffers (plus GCS_READ_CHUNK_BYTES for a gs:// input), so larger
# cohorts are merged in shards of at most this many inputs, and the shards' outputs merged after.
MAX_INPUTS_PER_MERGE = 256
# Bytes fetched per request when reading a gs:// file. Each open gs:// file holds one such chunk in memory.
GCS_READ_CHUNK_BYTES = 2 * 1024 * 1024
# RunInfo's CommandLine tells whether ExpansionHunter was run with --output-motif-composition. In JSONs from older
# ExpansionHunter versions RunInfo is among the last records. It holds the whole command line, which a --locus list of
# thousands of loci can make megabytes long, so while such a JSON is streamed everything from this key on is kept, up
# to MAX_RUN_INFO_BYTES (well past any command line the operating system accepts).
RUN_INFO_KEY = b'"RunInfo"'
MAX_RUN_INFO_BYTES = 64 * 1024 * 1024
MOTIF_COMPOSITION_OPTION = "--output-motif-composition"


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


def iterate_tsv_rows(path):
    """Yields {column: value} for each row of a TSV whose first line is its header."""
    with open_text(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def map_in_threads(function, items):
    """Returns [function(item) for item in items], computed in THREADS_FOR_GCS_REQUESTS threads (for functions that
    wait on Google Cloud Storage requests)."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=THREADS_FOR_GCS_REQUESTS) as executor:
        return list(executor.map(function, items))


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


def add_count(counts_by_key, key, count):
    counts_by_key[key] = counts_by_key.get(key, 0) + count


def format_histogram(histogram):
    """Formats {15: 120, 16: 3} (allele size in repeat units: count) as "15x:120,16x:3"."""
    return ",".join(f"{allele}x:{count}" for allele, count in sorted(histogram.items()))


def parse_histogram(histogram_string):
    """Parses "15x:120,16x:3" into {15: 120, 16: 3}."""
    histogram = {}
    for entry in filter(None, histogram_string.split(",")):
        allele, count = entry.split("x:")
        histogram[int(allele)] = int(count)
    return histogram


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


def format_count_with_noun(count, noun):
    """Returns the count and the noun, plural unless the count is 1: "1 input", "0 inputs", "1,024 inputs"."""
    return f"{count:,d} {noun}{'' if count == 1 else 's'}"


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


ExpansionHunterJsonHeader = collections.namedtuple("ExpansionHunterJsonHeader", [
    "sample_id", "sex", "loci_sorted_by_locus_id", "command_line"])


def read_expansion_hunter_json_header(json_path):
    """Returns an ExpansionHunterJsonHeader for an ExpansionHunter output JSON: (SampleId, Sex,
    loci_sorted_by_locus_id, command_line), where command_line is RunInfo.CommandLine ("" if RunInfo has none).

    ExpansionHunter writes SampleParameters, RunInfo, LocusResults and RunStats, in that order, in every mode. Its
    seeking and streaming modes write the loci sorted by LocusId, and its optimized-streaming mode nearly in genomic
    order. Older versions wrote RunInfo after the loci (with RunStats' fields in it): optimized-streaming output as
    SampleParameters, LocusResults, RunInfo, and seeking and streaming output with the keys sorted, as LocusResults,
    RunInfo, SampleParameters. Only optimized-streaming output starts with SampleParameters. Getting SampleParameters
    from older seeking or streaming output means decompressing the whole file, and both it and RunInfo are taken
    from the end of the file without parsing the loci (Illumina's ExpansionHunter, whose output has no RunInfo, needs
    them parsed); for older optimized-streaming output, command_line is returned as None (RunInfo is read from the end
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
        records_from_run_info_on = parse_records_from_run_info_on(read_tail_from_run_info(json_path))
        sample_parameters = records_from_run_info_on.get("SampleParameters")
        if run_info is None:
            run_info = records_from_run_info_on.get("RunInfo")
        if sample_parameters is None:
            sample_parameters = {}
            with open_binary(json_path) as f:
                for sample_parameters in ijson.items(f, "SampleParameters"):
                    break
    sample_id = sample_parameters.get("SampleId") or os.path.basename(strip_suffix(json_path, JSON_SUFFIXES))
    return ExpansionHunterJsonHeader(sample_id, sample_parameters.get("Sex", ""), loci_sorted_by_locus_id,
                                     None if run_info is None else run_info.get("CommandLine", ""))


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


def read_tail_from_run_info(json_path):
    """Decompresses a whole ExpansionHunter JSON without parsing it, and returns the tail RunInfoKeepingReader keeps:
    the bytes from its last RunInfo key to the end of the file (only the last few bytes if it has no RunInfo)."""
    with open_binary(json_path) as raw_file:
        f = RunInfoKeepingReader(raw_file)
        while f.read(1 << 20):
            pass
    return f.tail


def parse_records_from_run_info_on(json_tail):
    """Returns {name: record} for the top-level records from RunInfo to the end of an ExpansionHunter JSON, given the
    tail that RunInfoKeepingReader keeps, or {} if the tail doesn't start at RunInfo (the JSON has none, or more than
    MAX_RUN_INFO_BYTES follow it) or isn't the end of a complete JSON."""
    if not json_tail.startswith(RUN_INFO_KEY):
        return {}
    try:
        return json.loads(b"{" + json_tail)
    except ValueError:
        return {}


def run_info_has_motif_composition(run_info):
    """Returns whether a RunInfo record (or None) has --output-motif-composition in its CommandLine."""
    return run_info is not None and MOTIF_COMPOSITION_OPTION in run_info.get("CommandLine", "")


def iterate_repeat_records_of_json(path, run_info_comes_after_the_loci, handle_run_info_after_the_loci):
    """Streams an ExpansionHunter JSON and yields iterate_repeat_records' tuples for it.

    Args:
        path (str): the JSON, a local or gs:// path.
        run_info_comes_after_the_loci (bool): whether the JSON is from an older ExpansionHunter version that wrote
            RunInfo after the loci. Its RunInfo is then picked out while the loci are read, and passed to
            handle_run_info_after_the_loci once they all have been. (Doing that for a JSON whose RunInfo comes first
            would keep everything after it in memory.)
        handle_run_info_after_the_loci: function(RunInfo record, or None if the JSON has none).
    """
    with open_binary(path) as raw_file:
        f = RunInfoKeepingReader(raw_file) if run_info_comes_after_the_loci else raw_file
        yield from iterate_repeat_records(f)
        # iterate_repeat_records only finishes once ijson has read the whole file, so an older JSON's tail holds
        # RunInfo (its last record in optimized-streaming output; in seeking and streaming output only
        # SampleParameters comes after it)
        if run_info_comes_after_the_loci:
            handle_run_info_after_the_loci(parse_records_from_run_info_on(f.tail).get("RunInfo"))


def iterate_repeat_records(binary_file):
    """Yields (LocusId, VariantId, locus record, variant record) for each Repeat variant of an ExpansionHunter JSON,
    in file order. The locus record (LocusResults.<locus>) holds the locus-level fields, such as AlleleCount.

    The JSON is streamed with ijson, so it is never fully loaded.

    Args:
        binary_file: the JSON, open for reading bytes (decompressed).
    """
    for locus_id, locus_record in ijson.kvitems(binary_file, "LocusResults", use_float=True):
        for variant_id, variant_record in locus_record.get("Variants", {}).items():
            if variant_record.get("VariantType") == "Repeat":
                yield locus_id, variant_id, locus_record, variant_record


def iterate_in_genomic_order(keyed_items, input_path, temp_dir):
    """Puts the repeats of an optimized-streaming JSON back into genomic order.

    ExpansionHunter's optimized-streaming mode writes the chromosomes in the order of the BAM or CRAM header, and each
    chromosome's repeats nearly in genomic order. While each chromosome is followed by the next one in
    STANDARD_CHROMOSOME_ORDER, as in a header sorted that way, the repeats are put back in order through a buffer that
    holds the repeats of the last REORDER_WINDOW_BP on the chromosome. At the first chromosome that is not (such as
    chr10 after chr1, in a header sorted by name), the rest of the JSON is sorted with iterate_sorted_by_chromosome,
    which uses temporary files in temp_dir. A catalog without some chromosomes is sorted that way too, with the same
    result.

    Args:
        keyed_items: iterable of (get_genomic_sort_key of the repeat, item) in file order, where each item has a
            reference_region attribute (used in error messages).
        input_path (str): the JSON's path (used in error messages).
        temp_dir (str): directory for the temporary files.

    Yields:
        tuple: (sort key, item) in genomic order.

    Raises:
        ValueError: if a repeat comes more than REORDER_WINDOW_BP after repeats that start later, or a chromosome's
            repeats are split between the part of the JSON read in order and the part that was sorted.
    """
    keyed_items = iter(keyed_items)
    reorder_buffer = []
    counter = 0
    current_chrom_key = None
    largest_start = 0
    last_yielded_key = None
    for key, item in keyed_items:
        if key[0] != current_chrom_key:
            while reorder_buffer:
                last_yielded_key, _, buffered_item = heapq.heappop(reorder_buffer)
                yield last_yielded_key, buffered_item
            if key[0] != NEXT_STANDARD_CHROMOSOME_KEY.get(current_chrom_key):
                yield from iterate_sorted_by_chromosome(
                    itertools.chain([(key, item)], keyed_items), input_path, temp_dir, last_yielded_key)
                return
            current_chrom_key = key[0]
            largest_start = 0
        if last_yielded_key is not None and key < last_yielded_key:
            raise ValueError(
                f"{input_path} lists {key[4]} ({item.reference_region}) more than {REORDER_WINDOW_BP:,d} bp "
                f"after repeats that start later, which optimized-streaming output isn't expected to do.")
        heapq.heappush(reorder_buffer, (key, counter, item))
        counter += 1
        largest_start = max(largest_start, key[1])
        while reorder_buffer and reorder_buffer[0][0][1] < largest_start - REORDER_WINDOW_BP:
            last_yielded_key, _, buffered_item = heapq.heappop(reorder_buffer)
            yield last_yielded_key, buffered_item
    while reorder_buffer:
        key, _, buffered_item = heapq.heappop(reorder_buffer)
        yield key, buffered_item


def iterate_sorted_by_chromosome(keyed_items, input_path, temp_dir, last_yielded_key):
    """Yields (sort key, item) pairs that come in any order in genomic order: they are written to a temporary file per
    chromosome (pickled, in temp_dir), and then each chromosome's pairs are sorted in memory and yielded in turn, so
    memory use is about one chromosome's pairs.

    Args:
        keyed_items: iterable of (get_genomic_sort_key of the repeat, item).
        input_path (str): the JSON's path (used in error messages).
        temp_dir (str): directory for the temporary files.
        last_yielded_key: the key of the last pair already yielded from the same JSON, or None.

    Raises:
        ValueError: if a pair sorts before last_yielded_key.
    """
    with tempfile.TemporaryDirectory(dir=temp_dir) as bucket_dir:
        bucket_files = {}
        try:
            for key, item in keyed_items:
                if key[0] not in bucket_files:
                    bucket_files[key[0]] = open(os.path.join(bucket_dir, f"chromosome_{len(bucket_files)}.pickle"),
                                                "wb")
                pickle.dump((key, item), bucket_files[key[0]], protocol=pickle.HIGHEST_PROTOCOL)
        finally:
            for bucket_file in bucket_files.values():
                bucket_file.close()
        for chromosome_key in sorted(bucket_files):
            pairs = []
            with open(bucket_files[chromosome_key].name, "rb") as bucket_file:
                while True:
                    try:
                        pairs.append(pickle.load(bucket_file))
                    except EOFError:
                        break
            pairs.sort(key=operator.itemgetter(0))
            if last_yielded_key is not None and pairs[0][0] < last_yielded_key:
                raise ValueError(f"{input_path} lists the repeats of {get_chromosome(pairs[0][1].reference_region)} "
                                 f"in more than one place")
            yield from pairs
            last_yielded_key = pairs[-1][0]


def write_lines_in_genomic_order(keyed_lines, output_file, temp_dir, get_sort_key_of_line):
    """Writes lines that come in any order, such as the repeats of a seeking or streaming mode JSON (sorted by
    LocusId), to output_file in genomic order.

    Each line is written to a temporary file for its chromosome, and then each chromosome's lines are sorted in memory
    and written out, so memory use is about one chromosome's lines.

    Args:
        keyed_lines: iterable of (get_genomic_sort_key of the line's repeat, the line, ending with a newline).
        output_file: the open file to write the lines to.
        temp_dir (str): directory for the temporary files.
        get_sort_key_of_line: function(line) that returns the line's get_genomic_sort_key.
    """
    with tempfile.TemporaryDirectory(dir=temp_dir) as bucket_dir:
        bucket_files = {}
        try:
            for key, line in keyed_lines:
                if key[0] not in bucket_files:
                    bucket_files[key[0]] = open(os.path.join(bucket_dir, f"chromosome_{len(bucket_files)}.txt"), "w")
                bucket_files[key[0]].write(line)
        finally:
            for bucket_file in bucket_files.values():
                bucket_file.close()
        for chromosome_key in sorted(bucket_files):
            with open(bucket_files[chromosome_key].name) as bucket_file:
                lines = bucket_file.readlines()
            lines.sort(key=get_sort_key_of_line)
            output_file.writelines(lines)


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


def set_up_command_line_run():
    """Makes status messages reach a log file as they are printed, and raises the open-file limit."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(line_buffering=True)
    raise_open_file_limit()


def check_input_paths(parser, input_paths, allowed_suffixes):
    """Exits with an argparse error if an input path is listed twice, has none of the allowed suffixes, or doesn't
    exist."""
    if len(set(input_paths)) < len(input_paths):
        parser.error("An input path is listed more than once")
    for path in input_paths:
        if not path.endswith(allowed_suffixes):
            parser.error(f"{path} must end with one of {', '.join(allowed_suffixes)}")
    missing_paths = [path for path, exists in zip(input_paths, map_in_threads(file_exists, input_paths))
                     if not exists]
    if missing_paths:
        parser.error(f"{format_count_with_noun(len(missing_paths), 'input')} not found, e.g. {missing_paths[0]}")


def read_input_paths(parser, input_paths, input_paths_file, allowed_suffixes):
    """Returns the input paths given on the command line followed by those listed in the --input-paths-file (one per
    line; blank lines and lines that start with # are skipped), after checking them with check_input_paths.

    Args:
        parser (argparse.ArgumentParser): the parser whose error() reports a problem.
        input_paths (list): the input paths given on the command line.
        input_paths_file (str): the --input-paths-file, a local or gs:// path, or None.
        allowed_suffixes (tuple): the suffixes an input path may end with.

    Returns:
        list: the input paths.
    """
    input_paths = list(input_paths)
    if input_paths_file:
        if not file_exists(input_paths_file):
            parser.error(f"--input-paths-file {input_paths_file} not found")
        with open_text(input_paths_file) as f:
            input_paths += [line.strip() for line in f if line.strip() and not line.startswith("#")]
    if not input_paths:
        parser.error("No inputs: pass input paths on the command line, or with --input-paths-file")
    check_input_paths(parser, input_paths, allowed_suffixes)
    return input_paths


def read_ids_to_exclude(path, id_column):
    """Returns the IDs listed in a file: those in the id_column of a TSV whose header has that column, or else one per
    line. Blank lines and lines that start with # are skipped."""
    with open_text(path) as f:
        lines = [line.rstrip("\n") for line in f if line.strip() and not line.startswith("#")]
    if lines and id_column in lines[0].split("\t"):
        column_index = lines[0].split("\t").index(id_column)
        return frozenset(line.split("\t")[column_index] for line in lines[1:])
    return frozenset(line.strip() for line in lines)


def leave_out_samples(input_specs, sample_ids_to_exclude, samples_to_exclude_path):
    """Returns the input specs without the JSON inputs of the samples to exclude (those whose sample_id is listed), and
    prints how many were left out and how many of the listed samples aren't among the inputs.

    Raises:
        ValueError: if every input is left out.
    """
    kept_input_specs = [input_spec for input_spec in input_specs if input_spec.sample_id not in sample_ids_to_exclude]
    found_sample_ids = {input_spec.sample_id for input_spec in input_specs} & sample_ids_to_exclude
    print(f"Leaving out {format_count_with_noun(len(input_specs) - len(kept_input_specs), 'JSON input')} of the "
          f"{format_count_with_noun(len(sample_ids_to_exclude), 'sample')} listed in {samples_to_exclude_path}")
    missing_sample_ids = sorted(sample_ids_to_exclude - found_sample_ids)
    if missing_sample_ids:
        print(f"WARNING: {format_count_with_noun(len(missing_sample_ids), 'sample')} listed in "
              f"{samples_to_exclude_path} {'is' if len(missing_sample_ids) == 1 else 'are'} not among the inputs "
              f"(e.g. {missing_sample_ids[0]})")
    if not kept_input_specs:
        raise ValueError(f"Every input is a sample listed in {samples_to_exclude_path}")
    return kept_input_specs


def output_would_overwrite_an_input(output_path, input_paths):
    """Returns whether output_path is one of the local input paths (after resolving symbolic links)."""
    return os.path.realpath(output_path) in {os.path.realpath(path) for path in input_paths if not is_gcs_path(path)}


def check_output_prefix(parser, output_prefix):
    """Exits with an argparse error if --output-prefix is a gs:// path or its folder doesn't exist."""
    if is_gcs_path(output_prefix):
        parser.error("--output-prefix must be a local path; copy the output files to Google Cloud Storage afterwards")
    if not os.path.isdir(os.path.dirname(os.path.abspath(output_prefix))):
        parser.error(f"The output folder {os.path.dirname(os.path.abspath(output_prefix))} doesn't exist")


def create_temporary_directory_next_to(output_prefix):
    """Returns a TemporaryDirectory in the output's folder, named after the output, for the temporary TSVs."""
    return tempfile.TemporaryDirectory(prefix=os.path.basename(output_prefix) + ".tmp_",
                                       dir=os.path.dirname(os.path.abspath(output_prefix)))


def get_shard_tsv_paths(temp_dir, num_shards):
    return [os.path.join(temp_dir, f"shard_{index + 1}.tsv.gz") for index in range(num_shards)]


def split_inputs_into_shards(input_specs, threads, file_size_by_path):
    """Splits inputs into shards of similar total file size, one per thread, or more if a shard would otherwise have
    more inputs than one merge can read at once (get_max_inputs_per_merge()). Each shard gets at most
    ceil(N / number of shards) inputs.

    Args:
        input_specs (list): the inputs, each with a path attribute.
        threads (int): number of worker processes.
        file_size_by_path (dict): file size of each input's path.

    Returns:
        list: lists of input specs, one per shard.
    """
    num_shards = max(min(threads, len(input_specs)), math.ceil(len(input_specs) / get_max_inputs_per_merge()))
    max_shard_size = math.ceil(len(input_specs) / num_shards)
    shards = [[] for _ in range(num_shards)]
    shard_file_sizes = [0] * num_shards
    for input_spec in sorted(input_specs, key=lambda spec: (-file_size_by_path[spec.path], spec.path)):
        shards_with_room = [index for index in range(num_shards) if len(shards[index]) < max_shard_size]
        smallest_shard = min(shards_with_room, key=lambda index: (shard_file_sizes[index], index))
        shards[smallest_shard].append(input_spec)
        shard_file_sizes[smallest_shard] += file_size_by_path[input_spec.path]
    return shards


def map_in_worker_processes(function, tasks, num_workers, show_progress_bar, unit):
    """Runs function on each task, in up to num_workers worker processes (in this process if num_workers is 1).

    Leaving the pool (normally, or because a task raised) terminates the other workers.

    Args:
        function: a module-level function, so that worker processes can import it.
        tasks (list): its arguments, one per call.
        num_workers (int): number of worker processes.
        show_progress_bar (bool): whether to show a progress bar of finished tasks.
        unit (str): what a task is, for the progress bar.

    Returns:
        list: the results, in the order the tasks finished (a function whose results have to be matched to their task
            returns what identifies the task).
    """
    if num_workers <= 1:
        results = map(function, tasks)
        if show_progress_bar:
            results = tqdm.tqdm(results, total=len(tasks), unit=f" {unit}")
        return list(results)
    with multiprocessing.Pool(processes=num_workers) as pool:
        results = pool.imap_unordered(function, tasks)
        if show_progress_bar:
            results = tqdm.tqdm(results, total=len(tasks), unit=f" {unit}")
        return list(results)


def merge_files_in_rounds(paths, output_path, temp_dir, merge_files):
    """Merges files into one, first in rounds into temporary files if there are more than one merge can read at once
    (get_max_inputs_per_merge()). Each round's files, and the given ones, are deleted once merged.

    Args:
        paths (list): the files to merge.
        output_path (str): the merged file to write.
        temp_dir (str): directory for the merged files of each round.
        merge_files: function(list of paths, output path, is_final_merge) that merges files into one.

    Returns:
        what merge_files returns for the final merge.
    """
    max_inputs = get_max_inputs_per_merge()
    round_number = 0
    while len(paths) > max_inputs:
        round_number += 1
        next_round_paths = []
        for start in range(0, len(paths), max_inputs):
            next_round_path = os.path.join(temp_dir, f"round_{round_number}_part_{len(next_round_paths) + 1}.tsv.gz")
            merge_files(paths[start:start + max_inputs], next_round_path, False)
            next_round_paths.append(next_round_path)
        remove_files(paths)
        paths = next_round_paths
    result = merge_files(paths, output_path, True)
    remove_files(paths)
    return result


def merge_records_by_key(iterators, merge_into_first):
    """Merges iterators of (sort key, record), each sorted by key, merging the records that share a key.

    Args:
        iterators (list): the iterators to merge. Their keys are get_genomic_sort_key's.
        merge_into_first: function(first record, other record, description) that adds the other record's totals to
            the first's. The description names the repeat, for error messages.

    Yields:
        tuple: (sort key, merged record), in key order.
    """
    merged = heapq.merge(*iterators, key=operator.itemgetter(0))
    for key, records in itertools.groupby(merged, key=operator.itemgetter(0)):
        _, record = next(records)
        for _, other_record in records:
            merge_into_first(record, other_record, f"{key[4]} in locus {key[3]}")
        yield key, record


def iterate_checking_that_keys_increase(keyed_items, path):
    """Yields the (sort key, item) pairs of a file, in order.

    Raises:
        ValueError: if a key isn't larger than the one before it, meaning the file at path isn't sorted.
    """
    previous_key = None
    for key, item in keyed_items:
        if previous_key is not None and key <= previous_key:
            raise ValueError(f"{path} is not sorted: {key[3]} {key[4]} comes after {previous_key[3]} "
                             f"{previous_key[4]}")
        previous_key = key
        yield key, item
