[![Build](https://github.com/bw2/ExpansionHunter/actions/workflows/build.yml/badge.svg)](https://github.com/bw2/ExpansionHunter/actions/workflows/build.yml)
[![Docker Build](https://github.com/bw2/ExpansionHunter/actions/workflows/docker.yml/badge.svg)](https://github.com/bw2/ExpansionHunter/actions/workflows/docker.yml)

### ExpansionHunter fork - under active development

This modified version of ExpansionHunter introduces the following new features:
- **New analysis mode**:
  - `--analysis-mode optimized-streaming` is like `streaming` mode, but uses much less memory (< 10Gb on 30x coverage BAMs) and significantly speeds up analysis of large catalogs (> ~5k loci). It introduces 2 key modifications. 1) It uses simple heuristics to detect which loci can be confidently genotyped using only spanning reads. For these loci, it quickly computes their genotypes without running the full computationally-expensive ExpansionHunter genotyping algorithm. Since any given individual has fewer than ~5k to 10k large expansions in their genome (see [[Weisburd 2023](https://pubmed.ncbi.nlm.nih.gov/37214979/)]), this quick heuristic-based genotyping can take care of 90+% of loci in genome-wide catalogs, while maintaining overall accuracy within ~2% of the original algorithm (see <a href="https://broadinstitute.github.io/str-truth-set-v2/tool_comparison_viewer.html#?s=HG002&st=illumina&t=EHv5-bw2-optimized&c=31x&m=2to6bp&g=all&seqs=HG002&seqst=illumina&seqt=EHv5-bw2-optimized&seqc=31x&seqmet=edit_distance_bp&seqm=2to6bp&seqg=all">tool benchmark</a>). Adding `--genotyping-approach only-full` turns the quick heuristic off and runs the full genotyping algorithm on every locus (this was previously provided as a separate `low-mem-streaming` analysis mode). 2) It keeps memory usage low and independent of catalog size by having ExpansionHunter read the BAM/CRAM file in two sequential passes. The first pass caches all mate pairs that aligned far away (> ~2kb) from each other, while the second pass genotypes all loci. Since all far-away reads that may be needed for genotyping a locus are already cached in memory, the second pass can ingest the locally-aligned reads, genotype the locus, and then immediately discard these reads from memory before moving on to the next locus. This avoids having to keep reads from all loci in memory before genotyping begins. Since reading through a file sequentially is a relatively fast operation (taking minutes rather than hours), this two-pass approach isn't much slower than the single-pass approach used by the original `streaming` mode. 
  - `June 21, 2026`: reservoir sampling further reduces processing time by downsampling reads at loci with extremely high read depth (see `--max-depth`).
- **Integrated read visualizations**: REViewer functionality is now built directly into ExpansionHunter, outputting SVG read pileup images without needing a separate post-processing step (see [VariantCatalog docs](docs/04_VariantCatalogFiles.md)).
  - `--plot-all` generates read visualizations for every locus
  - `--disable-all-plots` disables all image generation (overrides catalog settings)
  - `PlotReadVisualization` field in the variant catalog enables conditional image generation based on genotype thresholds (e.g., only visualize when long allele >= 400 repeats)
- **Per-allele quality metrics**: New `AlleleQualityMetrics` in JSON output provides detailed quality information for each allele (see [AlleleQualityMetrics docs](docs/07_AlleleQualityMetrics.md)).
  - Metrics include QD (quality by depth), strand bias, flank depth, insertion/deletion rates, and more
  - These metrics were used to train an **Allele Quality Score** model on a genome-wide tandem repeat truth set derived from HG002 and CHM1_CHM13 T2T assemblies. This model adds 3 probabilities to the output -  `pOk`, `pTooShort`, and `pTooLong` - predicting whether the called allele size is close to, shorter than, or longer than the true allele size. It also outputs the `PredictedLengthCorrectionFactor` - dividing ExpansionHunter's allele size by this factor should bring it closer to the true allele size, especially for the larger alleles that have 0 spanning reads and where `pOk` < 0.5 (see [AlleleQualityMetrics docs](docs/07_AlleleQualityMetrics.md#genotype-quality-model-predictions)).
- **Consensus allele sequences**: Consensus nucleotide sequences are now reported for each allele. This is a simplistic first implementation that just collapses confidently-placed (ie. darker-colored) reads within the REViewer visualization and takes the most common base at each position. Insertions and deletions within the reads are not incorporated into the consensus sequence. Also, any positions not covered by confidently-placed reads are reported as N's (see [Consensus Sequences docs](docs/05_OutputJsonFiles.md#consensus-sequences)).
  - `--dont-output-consensus-sequences` disables consensus sequence computation if not needed
  - `--dont-output-quality-metrics` disables quality metrics computation if not needed
- **Motif composition**: `--output-motif-composition all-loci` (or `loci-with-non-ref-motifs`, or `loci-with-known-motifs`) outputs how many times different motifs (including the catalog motif) were observed in the reads aligned to each locus. These counts are written to a new "MotifComposition" section in the output JSON and are reported per locus. When the two allele sizes in a genotype differ by 2 or more repeats and reads can be unambiguously assigned to one or the other allele based on repeat length, per-allele counts are also reported. Optionally, a "KnownMotifs" list can be added to the locus definition in the input catalog to guide motif discovery (see [Motif composition docs](docs/05_OutputJsonFiles.md#motif-composition)). Enabling motif composition output increases runtime by ~5% and doesn't affect memory usage. The output JSON file size increases by 5% to 10%.
- **Misc. new convenience features and options**:
  - supports gzip-compressed input catalogs, and provides a `-z` option to compress the output files
  - **Converts N chars error to a warning**: changes `Flanks can contain at most 5 characters N but found x Ns` from an error to a warning. ExpansionHunter now just prints this warning and skips the offending locus instead of exiting.
  - `--start-with`, `--n-loci`, and `--sort-catalog-by` options allow processing a fixed number of loci from the input catalog
  - `--locus` for filtering the input catalog to specific LocusId(s)
  - `--reads-index` explicitly specifies the BAM/CRAM index file path or URL, useful when the index is in a different location than the reads file or when auto-detection doesn't work with cloud URLs
  - `--region` for filtering the input catalog to a specific genomic interval
  - `--skip-hom-ref` doesn't output results for loci where the genotype is homozygous reference, reducing output file size
  - `--skip-missing-genotypes` doesn't output results for loci with a missing genotype (eg. due to low coverage)
  - `--copy-catalog-fields` copies extra annotation fields (e.g., Gene, Diseases) from the input catalog to the output JSON. By default, ExpansionHunter simply ignores custom fields in the input catalog. 
  - `--enable-bamlet-output` writes a "bamlet" BAM file containing the realigned reads for each locus
  - `--genotyping-approach <auto|only-quick|only-full>` selects which genotyper(s) `optimized-streaming` mode runs on each locus. `auto` (the default) uses the quick spanning-read heuristic where it can confidently genotype a locus, while using the full genotyping algorithm elsewhere. `only-quick` skips full genotyping completely and is provided mainly for benchmarking or debugging purposes. `only-full` runs the full genotyping algorithm on every locus (previously `only-full` was implemented as a separate `--analysis-mode` called `low-mem-streaming`).
  - `--cache-mates` enables a cross-locus read cache in `seeking` analysis mode to make it run faster on catalogs where many loci have the same motif (eg. if your catalog contains mostly `CGG` and `CCG` repeats). Since in-repeat reads from all these loci will typically mismap to the same few places in the genome, caching the reads in-memory can subsantially reduce disk access latency. For large catalogs (> ~5k loci), it is still better to use `optimized-streaming`. 
  - `--max-depth` (default `150`) limits the number of reads processed per locus in `optimized-streaming` mode using reservoir sampling. This bounds memory and runtime at extremely high-coverage loci (e.g. centromeric/satellite repeats) where millions of reads can pile up and slow down processing. Set to `0` to disable the cap.
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
