# Usage

Expansion Hunter requires the following inputs:
1. A BAM or a CRAM file containing aligned reads from a PCR-free WGS sample.
    1. The BAM or CRAM file must be sorted, and must be indexed for every [analysis mode](#analysis-modes) except plain `streaming`.
    2. The BAM or CRAM file may be a local filesystem path or [URL](#url-support).
4. A FASTA file with a reference genome assembly (which must be the same as the one used to align the reads)
5. A [variant catalog file](04_VariantCatalogFiles.md).

Expansion Hunter outputs a VCF file and a JSON file with variant genotypes and
other useful information. It can also optionally output a BAMlet containing
alignments of reads that overlap or are located in close proximity to each
variant, but only if the `--enable-bamlet-output` flag is specified (see
[Optional arguments](#optional-arguments) below). The VCF and JSON files
are largely equivalent, but the JSON file may be easier to parse
programmatically. Here is a template with the names of the required parameters.

```bash
ExpansionHunter --reads <aligned reads BAM/CRAM file/URL> \
                --reference <reference genome FASTA file> \
                --catalog <JSON file specifying variants to genotype> \
                --output-prefix <Prefix for the output files>
```

## Optional arguments

In addition to the required program options listed above, there are a number of
optional arguments.

* `--sex <arg>` Specifies sex of the sample; can be either `male` or `female`
  (default). This parameter only affects repeats on sex chromosomes.
* `--threads <int>` Specifies how many threads to can be used accelerate analysis
   of large variant catalogs. Set to 1 by default. Typically seeking mode can
   benefit from relatively high thread counts, while for streaming mode
   there is limited benefit beyond about 16 threads.
* `--min-locus-coverage <int>` Specifies minimum read coverage depth at loci
   on diploid chromosomes required to attempt genotyping. Automatically reduced
   to half for loci on haploid chromosomes. The locus will be skipped if the
   coverage falls below this value. Set to 10 by default.
* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.
* `--max-depth <int>` In `low-mem-streaming` and `optimized-streaming` modes,
   this sets a limit on the number of reads processed per locus using reservoir sampling.
   The intention is to bound the memory usage and runtime at extremely high-coverage loci
   (e.g. centromeric/satellite repeats) where millions of reads can pile up and slow down processing.
   Set to 100 by default; set to 0 to disable the cap.
* `--reads-index <BAM/CRAM index file/URL>` Specifies the BAM/CRAM index file
  path or URL explicitly, instead of auto-detecting it from the `--reads` path.
  This is useful when the index file is in a different location than the reads
  file, or when using cloud URLs where auto-detection may not work.
* `--analysis-mode <mode>` Specify analysis mode, which can be `seeking`,
  `streaming`, `low-mem-streaming`, or `optimized-streaming`. The default mode
  is `seeking`. See further description of analysis modes below.
* `--dont-output-quality-metrics` Disable per-allele quality metrics computation. By
  default, ExpansionHunter computes quality metrics (QD, strand bias, flank depth,
  etc.) for each allele and outputs them in the JSON file. Use this flag to skip
  this computation if the metrics are not needed.
* `--enable-bamlet-output` Output a BAM file containing realigned reads
  that overlap or are located in close proximity to each variant. The file is
  written to `<output-prefix>_realigned.bam`.
* `--copy-catalog-fields` Copy extra annotation fields from the input variant
  catalog to the output JSON. This allows custom fields like `Gene`, `Diseases`,
  `PathogenicMin`, etc. to be preserved in the output, making it easier to
  annotate results without a separate join step.
* `--resume` Make an interrupted run restartable. See [Resuming an interrupted
  run](#resuming-an-interrupted-run) below.


Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.

### URL support

The aligned reads input BAM or CRAM file may be a local filesystem path or URL.
Supported protocols for URL input include ftp, https, s3, and gs (Google Cloud
Storage). S3 and GCS bucket access can be configured using the URL syntax and
environment variables supported by samtools/htslib.

### Analysis modes

#### Seeking mode

In seeking mode, alignment file indexing is used to seek specific read sets for the
analysis of each variant. Seeking mode is recommended for analysis of small catalogs.
This mode requires that the input BAM or CRAM file is already sorted and indexed.

#### Streaming mode

In streaming mode, the alignment file is read in a single pass and all variants are
analyzed during this reading operation. Streaming mode is recommended for the analysis
of large catalogs, but does require more memory as a funciton of catalog size. This mode
does not require that the BAM or CRAM file is sorted or indexed.

#### Low-mem-streaming mode

Changes how data is read from the input BAM or CRAM file in order to keep memory usage (typically < 10 GB) and is independent of catalog size.
The output stays nearly identical to `streaming` mode.

#### Optimized-streaming mode

`optimized-streaming` mode uses a fast heuristic genotyper to identify loci that can be quickly genotyped using spanning reads. It then runs the full graph-based genotyper
only on the subset of loci that appear to have larger expansions. This significantly speeds up analysis of large
catalogs (> ~10k loci) since the majority of loci can be genotyped using only spanning reads. Memory usage is similar to `low-mem-streaming` mode.

### Resuming an interrupted run

A run that is killed part way through (an out-of-memory kill, a preempted machine, Ctrl-C, a crash)
normally has to be started over from the first locus. Adding `--resume` makes it restartable:

```bash
ExpansionHunter --reads sample.cram --reference reference.fa --catalog catalog.json \
  --output-prefix sample --analysis-mode optimized-streaming --threads 8 --resume
```

As each locus is genotyped it is recorded in three checkpoint files next to the output:

| File | Contents |
|------|----------|
| `<output-prefix>.json[.gz].unfinished` | each finished locus's JSON record |
| `<output-prefix>.vcf[.gz].unfinished`  | each finished locus's VCF lines |
| `<output-prefix>.processed_loci.unfinished` | one finished locus per line, with how many records it wrote to each file above |

All three are deleted once the final `.json` and `.vcf` files have been written, so a completed run
leaves nothing extra behind. Re-running the same command with `--resume` picks up the checkpoint
files, skips the loci they already contain, and genotypes only the rest. The final output is
identical to what an uninterrupted run would have produced.

Things worth knowing:

* `--resume` only works with `--analysis-mode low-mem-streaming` and `optimized-streaming`, which
  write their output as they go. `seeking` and `streaming` write everything only after the last
  locus, so there is nothing to resume from; passing `--resume` there logs a warning and changes
  nothing.
* The resumed run must use the same catalog, reads, sample and output-affecting options as the
  interrupted one. If they differ, ExpansionHunter stops with an error naming the option that
  changed, rather than silently mixing results from two different runs. Delete the `.unfinished`
  files (or use a different `--output-prefix`) to start over.
* `--threads` may differ between the two runs, with one exception. A `--threads 1` run can be resumed
  with any thread count, and a `--threads > 1` run with any count above 1; every finished locus is
  reused. Only resuming a `--threads > 1` run with `--threads 1` is refused, and only when that run
  finished contigs out of order (it usually has), since one coordinate sweep cannot append to that. The
  run then stops with an error and leaves the checkpoint untouched, so re-running with `--threads > 1`
  still recovers everything.
* Resuming skips the genotyping work, not the read scanning: the BAM/CRAM is still read from the
  beginning. On large catalogs genotyping dominates, so this is still a large saving.
* Checkpointing costs one extra write of each record, and needs extra disk while the run is in
  progress: the checkpoint files (compressed if `-z` is used), plus the per-slice temp files, which are
  **never compressed**. At `--threads > 1` those temps are the per-contig files the merge already used
  before `--resume` existed; at `--threads 1`, `--resume` adds one covering the whole output. So with
  `-z` the peak extra disk is roughly one uncompressed copy of the output plus one compressed copy.
* Every `LocusId` in the catalog must be unique. Resume tells loci apart by their id, so a catalog with
  duplicates is rejected with `--resume` (it already produces a colliding JSON record without it).
* Keep `-z` the same across the interrupted run and its resume. The checkpoint file names carry the
  output's `.gz` suffix; changing `-z` is detected and reported rather than silently starting over.
* `--resume` cannot be combined with `--enable-bamlet-output`. The bamlet is rewritten from scratch
  on every run and carries no record of which loci it covers, so a resumed run would replace it with
  one holding only the loci that run genotyped.

#### Known limitations of `low-mem-streaming` and `optimized-streaming`

These two newer modes ignore `OfftargetRegions` entries in the variant catalog. This can affect loci that do
explicitly list off-target regions in the catalog, such as **C9ORF72**, **FMR1**. For these loci,
`--analysis-mode seeking` or `--analysis-mode streaming` are recommended.