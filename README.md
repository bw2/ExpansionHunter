### ExpansionHunter-dev

This modified version of ExpansionHunter introduces the following new features:

- `--analysis-mode low-mem-streaming` which is like `streaming` mode and produces essentially the same output, but uses much less memory 
- support for gzip-compressed input catalogs, and provides a `-z` option to compress the output files
- `--start-with`, `--n-loci`, and `--sort-catalog-by` options to process only a fixed number of loci from the input catalog
- `--locus` to only process catalog loci with the given LocusId(s) 
- `--region` to only process catalog loci within a specific genomic region
- changes the `Flanks can contain at most 5 characters N but found x Ns` error to a warning, allowing ExpansionHunter to run to completion without terminating on these errors
- allows direct access to remote BAM/CRAM or reference FASTA files in Google Cloud Storage or S3
  - for access to private buckets, set environment variable:  
    `export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)`
  - for access to requester-pays buckets, also set environment variable  
    `export GCS_REQUESTER_PAYS_PROJECT=<your gcloud project>`
- optimization of the default "seeking" analysis mode to make it 1.5x to 3x faster without changing the output
  - it works by introducing an in-memory read cache that reduces the number of disk operations needed to retrieve mismapped mate pairs
  - by default, the cache is reset for each locus, leading to a modest speedup with negligible memory overhead
  - the new `--cache-mates` option activates reuse of the cache across loci, leading to a more significant speed increase, though at a cost of increased memory usage (typically in the range of 1-2GB of memory usage for catalogs with 100s to 1000s of loci)
  - for large catalogs, it is better to use the new "low-mem-streaming" analysis mode. However, if you do want to split a larger variant catalog into multiple shards and use "seeeking" mode, it's important to presort the catalog by normalized motif (which is the cyclic shift of a motif that is alphabetically first - ie. AGC rather than CAG). This ensures that loci with the same normalized motif will be processed in the same shard, increasing cache hit rates and therefore speed due to this optimization.


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
