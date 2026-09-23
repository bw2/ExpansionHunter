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
   Set to 150 by default (500 when `--output-motif-composition` is used); set to 0 to disable the cap.
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
* `--output-motif-composition <all-loci|loci-with-non-ref-motifs>` Add a `MotifComposition` record to
  the JSON output: counts of each motif and of each pair of adjacent motifs, for the locus and, when the two alleles differ by at least 2 repeat units, for each allele
  (see [Motif composition](05_OutputJsonFiles.md#motif-composition)). `all-loci` adds it to every
  genotype whose motif size is 2 bp or longer and at most a third of the read length; `loci-with-non-ref-motifs`
  adds it only where the reads show a motif not present in the reference repeat sequence. This option works with the
  `optimized-streaming` and `low-mem-streaming` analysis modes. When this flag is used and `--max-depth` is not specified, the
  default `--max-depth` is raised to 500 in order to better capture rare motifs at high-coverage loci.
* `--resume` Make an interrupted run restartable. Works with the `optimized-streaming` and
  `low-mem-streaming` analysis modes. See [Resuming an interrupted run](#resuming-an-interrupted-run)
  below.


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
normally has to be started over from the first locus. Adding `--resume` makes it restartable.
`--resume` works with `--analysis-mode optimized-streaming` and `--analysis-mode low-mem-streaming`:

```bash
ExpansionHunter --reads sample.cram --reference reference.fa --catalog catalog.json \
  --output-prefix sample --analysis-mode optimized-streaming --threads 8 --resume
```

Resuming requires the following intermediate files to be kept until ExpansionHunter is run again
(for example, copied back if the rerun happens on a different machine):

```
<output-prefix>.processed_loci.txt
<output-prefix>.contig<N>.json
<output-prefix>.contig<N>.vcf
```

* `<output-prefix>.processed_loci.txt` lists the loci that are finished.
* `<output-prefix>.contig<N>.json` and `<output-prefix>.contig<N>.vcf` hold the results of the
  finished loci on one chromosome.

If `<output-prefix>.processed_loci.txt` is missing, the run starts over from the first locus. If a
chromosome's temp files are missing, that chromosome's loci are genotyped again; other chromosomes
are unaffected.

The temp files and the list are deleted once the final `.json` and `.vcf` files have been written, so
a completed run leaves nothing extra behind. The final output is identical to what an uninterrupted
run would have produced.

Things worth knowing:

* `--resume` does not work with `seeking` or `streaming` mode. Those modes write everything only
  after the last locus, so there is nothing to resume from; passing `--resume` with them is an
  error, and the run stops before doing anything.
* The resumed run must use the same catalog, reads, reference, sample, ExpansionHunter build and
  output-affecting options as the interrupted one. Otherwise, ExpansionHunter stops with an error.
  Delete `<output-prefix>.processed_loci.txt` (or use a different `--output-prefix`) to start over.
* Resuming skips the genotyping work, not the read scanning: the BAM/CRAM is still read from the
  beginning. On large catalogs genotyping dominates, so this is still a large saving.
* While the run is in progress, the temp files take roughly one uncompressed copy of the output on
  disk (they are never compressed, even with `-z`). At `--threads > 1` this is no different from a
  run without `--resume`. At `--threads 1` it is extra, because without `--resume` that mode writes
  the final files directly and skips the merge.
* Every `LocusId` in the catalog must be unique. Resume tells loci apart by their id, so a catalog with
  duplicates is rejected with `--resume` (it already produces a colliding JSON record without it).
* `--resume` cannot be combined with `--enable-bamlet-output`. The bamlet is rewritten from scratch
  on every run and carries no record of which loci it covers, so a resumed run would replace it with
  one holding only the loci that run genotyped.

#### Known limitations of `low-mem-streaming` and `optimized-streaming`

These two newer modes ignore `OfftargetRegions` entries in the variant catalog. This can affect loci that do
explicitly list off-target regions in the catalog, such as **C9ORF72**, **FMR1**. For these loci,
`--analysis-mode seeking` or `--analysis-mode streaming` are recommended.

Their VCF output is not always sorted by position within a chromosome. Loci are written in the order
they finish. `bcftools index` and `tabix` require sorted input, so sort the VCF before indexing it:

```bash
bcftools sort -Oz -o sample.sorted.vcf.gz sample.vcf
bcftools index -t sample.sorted.vcf.gz
```

The same commands work on the `-z` output (`sample.vcf.gz`).
