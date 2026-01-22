### ExpansionHunter fork - under active development

This modified version of ExpansionHunter introduces the following new features:
- **Converts N chars error to warning**: changes the `Flanks can contain at most 5 characters N but found x Ns` error to a warning, allowing ExpansionHunter to run to completion without terminating on these errors
- **New analysis modes**:
  - `--analysis-mode low-mem-streaming` is like `streaming` mode and produces nearly identical output, but uses much less memory.
  - `--analysis-mode optimized-streaming` significantly speeds up analysis of large catalogs (> ~10k loci) by uses simple heuristics to detect which loci are almost certainly homozygous reference, and avoids running the full computationally-expensive genotyping algorithm on them. Its memory usage is also low (< 10Gb) and independent of catalog size, similar to `low-mem-streaming` mode.
- **Integrated read visualizations**: REViewer functionality is now built directly into ExpansionHunter, outputting SVG read pileup images without needing a separate post-processing step (see [VariantCatalog docs](docs/04_VariantCatalogFiles.md)).
  - `--plot-all` generates read visualizations for every locus
  - `--disable-all-plots` disables all image generation (overrides catalog settings)
  - `PlotReadVisualization` field in the variant catalog enables conditional image generation based on genotype thresholds (e.g., only visualize when long allele >= 400 repeats)
- **Consensus allele sequences**: Consensus nucleotide sequences are now reported for each allele. This is a simplistic first implementation that just collapses confidently-placed (ie. darker-colored) reads within the REViewer visualization and takes the most common base at each position. Insertions and deletions within the reads are not incorporated into the consensus sequence. Also, any positions not covered by confidently-placed reads are reported as N's (see [Consensus Sequences docs](docs/05_OutputJsonFiles.md#consensus-sequences)).
- **Per-allele quality metrics**: New `AlleleQualityMetrics` in JSON output provides detailed quality information for each allele (see [AlleleQualityMetrics docs](docs/07_AlleleQualityMetrics.md)).
  - Metrics include QD (quality by depth), strand bias, flank depth, insertion/deletion rates, and more
  - `--disable-quality-metrics` disables quality metrics computation if not needed
- **Misc. new convenience features and options**:
  - supports gzip-compressed input catalogs, and provides a `-z` option to compress the output files
  - `--start-with`, `--n-loci`, and `--sort-catalog-by` options allow processing a fixed number of loci from the input catalog
  - `--locus` for filtering the input catalog to specific LocusId(s) 
  - `--region` for filtering the input catalog to a specific genomic region
  - `--skip-hom-ref` skips output of loci where all variants are homozygous reference, reducing output file size
  - `--copy-catalog-fields` copies extra annotation fields (e.g., Gene, Diseases) from the input catalog to the output JSON
- **Input BAM or FASTA can be read directly from cloud buckets**: allows direct access to remote BAM/CRAM or reference FASTA files in Google Cloud Storage or S3 via functionality provided by htslib 
  - for access to private buckets, set environment variable:  
    `export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`
  - for access to requester-pays buckets, also set environment variable  
    `export GCS_REQUESTER_PAYS_PROJECT=<your gcloud project>`
- **Cache to speed up seeking mode**: `--cache-mates` option makes `--analysis-mode seeking` run 2x to 3x faster without changing the output
  - for large catalogs, it is better to use the new "low-mem-streaming" analysis mode. However, if you do want to split a larger variant catalog into multiple shards and then process them using "seeking" mode with `--cache-mates`, it's important to presort the catalog by normalized motif (the alphabetically-first cyclic shift of a motif - ie. AGC rather than CAG). This ensures that loci with the same motif will be processed in the same shard, increasing cache hit rates and therefore speed due to this optimization.

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
