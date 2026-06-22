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
   limits the number of reads processed per locus using reservoir sampling, to
   bound memory and runtime at extremely high-coverage loci (e.g.
   centromeric/satellite repeats) where millions of reads can otherwise pile up.
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

#### Known limitations of `low-mem-streaming` and `optimized-streaming`

These two newer modes ignore `OfftargetRegions` entries in the variant catalog. This can affect loci that do
explicitly list off-target regions in the catalog, such as **C9ORF72**, **FMR1**. For these loci,
`--analysis-mode seeking` or `--analysis-mode streaming` are recommended.