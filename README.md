[![Build](https://github.com/bw2/ExpansionHunter/actions/workflows/build.yml/badge.svg)](https://github.com/bw2/ExpansionHunter/actions/workflows/build.yml)

### ExpansionHunter fork - under active development

This modified version of ExpansionHunter introduces the following new features:
- **New analysis modes**:
  - `--analysis-mode low-mem-streaming` is like `streaming` mode and produces nearly identical output, but uses much less memory (< 10Gb). It achieves this by having ExpansionHunter read the BAM/CRAM file in two sequential passes.  The first pass caches all mate pairs that aligned far away (> ~2kb) from each other, while the second pass genotypes all loci. Since all far-away reads that may be needed for genotyping a locus are already cached in memory, the second pass can ingest the locally-aligned reads around a locus, genotype the locus, and then immediately discard these reads from memory before moving on to the next locus. This avoids having to keep reads from all loci in memory before genotyping begins. Since reading through a file sequentially is a relatively fast operation (taking minutes rather than hours), this two-pass approach isn't significantly slower than the single-pass approach of the original `streaming` mode.  
  - `--analysis-mode optimized-streaming` significantly speeds up analysis of large catalogs (> ~10k loci) by using simple heuristics to detect which loci can be confidently genotyped using only spanning reads, and then quickly computing their genotypes without running the full computationally-expensive ExpansionHunter genotyping algorithm. Since any given individual in the population will have no more than ~5k to 10k large expansions relative to the reference genome (see [[Weisburd 2023](https://pubmed.ncbi.nlm.nih.gov/37214979/)]), while genome-wide catalogs can have hundreds of thousands or millions of TR loci, this quick heuristic-based genotyping can be used for the majority of loci, yielding a 3x or more speedup depending on the catalog. The memory usage of this mode is also low (< 10Gb) and independent of catalog size, similar to `low-mem-streaming` mode.
  - `June 11, 2026`: these two new modes now fully support multi-threading via `--threads N`
- **Integrated read visualizations**: REViewer functionality is now built directly into ExpansionHunter, outputting SVG read pileup images without needing a separate post-processing step (see [VariantCatalog docs](docs/04_VariantCatalogFiles.md)).
  - `--plot-all` generates read visualizations for every locus
  - `--disable-all-plots` disables all image generation (overrides catalog settings)
  - `PlotReadVisualization` field in the variant catalog enables conditional image generation based on genotype thresholds (e.g., only visualize when long allele >= 400 repeats)
- **Consensus allele sequences**: Consensus nucleotide sequences are now reported for each allele. This is a simplistic first implementation that just collapses confidently-placed (ie. darker-colored) reads within the REViewer visualization and takes the most common base at each position. Insertions and deletions within the reads are not incorporated into the consensus sequence. Also, any positions not covered by confidently-placed reads are reported as N's (see [Consensus Sequences docs](docs/05_OutputJsonFiles.md#consensus-sequences)).
  - `--dont-output-consensus-sequences` disables consensus sequence computation if not needed
- **Per-allele quality metrics**: New `AlleleQualityMetrics` in JSON output provides detailed quality information for each allele (see [AlleleQualityMetrics docs](docs/07_AlleleQualityMetrics.md)).
  - Metrics include QD (quality by depth), strand bias, flank depth, insertion/deletion rates, and more
  - `--dont-output-quality-metrics` disables quality metrics computation if not needed
- **Misc. new convenience features and options**:
  - supports gzip-compressed input catalogs, and provides a `-z` option to compress the output files
  - **Converts N chars error to a warning**: changes the `Flanks can contain at most 5 characters N but found x Ns` error to a warning, allowing ExpansionHunter to run to completion without terminating on these errors
  - `--start-with`, `--n-loci`, and `--sort-catalog-by` options allow processing a fixed number of loci from the input catalog
  - `--locus` for filtering the input catalog to specific LocusId(s)
  - `--reads-index` explicitly specifies the BAM/CRAM index file path or URL, useful when the index is in a different location than the reads file or when auto-detection doesn't work with cloud URLs
  - `--region` for filtering the input catalog to a specific genomic region
  - `--skip-hom-ref` skips output of loci where all variants are homozygous reference, reducing output file size
  - `--skip-missing-genotypes` skips output of loci with missing genotypes (eg. due to low coverage)
  - `--copy-catalog-fields` copies extra annotation fields (e.g., Gene, Diseases) from the input catalog to the output JSON
  - `--enable-bamlet-output` writes a "bamlet" BAM file containing the realigned reads for each locus
  - `--quick-heuristic-genotyping-only` modifies `optimized-streaming` mode so that it only genotypes loci that can be confidently genotyped using spanning reads, while skipping full genotyping completely. Provided mainly for benchmarking or debugging purposes.
  - `--cache-mates` enables a cross-locus read cache in `seeking` analysis mode to make it run faster on catalogs where many loci have the same motif (eg. if you have a catalog of only/mostly `CGG` and `CCG` repeats). Since in-repeat reads from all these loci will typically mismap to the same few places in the genome, caching the reads in-memory can subsantially reduce disk access latency. For large catalogs (> 10k loci), it is still better to use `low-mem-streaming` or `optimized-streaming`. 
  - `--max-depth` (default `100`) limits the number of reads processed per locus in `low-mem-streaming` and `optimized-streaming` modes using reservoir sampling. This bounds memory and runtime at pathological high-coverage loci (e.g. centromeric/satellite repeats) where millions of reads can otherwise pile up. The cap is per-locus and scales with the locus window width, so it limits all loci to the same depth rather than the same absolute read count; loci below the cap are unaffected and the retained sample is deterministic and identical across `--threads`. Set to `0` to disable the cap.
- **Input BAM or FASTA can be read directly from cloud buckets**: allows direct access to remote BAM/CRAM or reference FASTA files in Google Cloud Storage or S3 via functionality provided by htslib 
  - for access to private buckets, set environment variable:  
    `export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`
  - for access to requester-pays buckets, also set environment variable  
    `export GCS_REQUESTER_PAYS_PROJECT=<your gcloud project>`


Thank you to [@maarten-k](https://github.com/maarten-k) for testing out early versions and introducing substantial optimizations to the build process.

### Citation
If you use this modified version of ExpansionHunter, please cite:
```
Insights from a genome-wide truth set of tandem repeat variation
Ben Weisburd, Grace Tiao, Heidi L. Rehm
bioRxiv 2023.05.05.539588; doi: https://doi.org/10.1101/2023.05.05.539588
```

---


# Expansion Hunter: a tool for estimating repeat sizes

There are a number of regions in the human genome consisting of repetitions of
short unit sequence (commonly a trimer). Such repeat regions can expand to a
size much larger than the read length and thereby cause a disease.
[Fragile X Syndrome](https://en.wikipedia.org/wiki/Fragile_X_syndrome),
[ALS](https://en.wikipedia.org/wiki/Amyotrophic_lateral_sclerosis), and
[Huntington's Disease](https://en.wikipedia.org/wiki/Huntington%27s_disease)
are well known examples.

Expansion Hunter aims to estimate sizes of such repeats by performing a targeted
search through a BAM/CRAM file for reads that span, flank, and are fully
contained in each repeat.

Linux and macOS operating systems are currently supported.

## License

Expansion Hunter is provided under the terms and conditions of the
[Apache License Version 2.0](LICENSE.txt). It relies on several third party
packages provided under other open source licenses, please see
[COPYRIGHT.txt](COPYRIGHT.txt) for additional details.

## Documentation

Installation instructions, usage guide, and description of file formats are
contained in the [docs folder](docs/01_Introduction.md).

## Companion tools and resources

- [A genome-wide STR catalog](https://github.com/Illumina/RepeatCatalogs)
  containing polymorphic repeats with similar properties to known pathogenic and
  functional STRs
- [REViewer](https://github.com/Illumina/REViewer), a tool for visualizing
  alignments of reads in regions containing tandem repeats

## Method

The method is described in the following papers:

- Egor Dolzhenko, Joke van Vugt, Richard Shaw, Mitch Bekritsky, and others,
  [Detection of long repeat expansions from PCR-free whole-genome sequence data](http://genome.cshlp.org/content/27/11/1895),
  Genome Research 2017

- Egor Dolzhenko, Viraj Deshpande, Felix Schlesinger, Peter Krusche, Roman Petrovski, and others,
[ExpansionHunter: A sequence-graph based tool to analyze variation in short tandem repeat regions](https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btz431/5499079),
Bioinformatics 2019
