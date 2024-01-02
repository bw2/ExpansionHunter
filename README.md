## Note 

This modified version of ExpansionHunter introduces the following new features:
- support for gzip-compressed input catalogs
- changes the `Flanks can contain at most 5 characters N but found x Ns` error to a warning. 
  - This allows ExpansionHunter to run to completion without exiting on these loci and makes it easier to process large catalogs without having to find and exclude these loci first.
- optimization of ExpansionHunter's "seeking" analysis mode that yields a 1.5x to 3x speed increase without changing the output.
  - it works by introducing an in-memory read cache that reduces the number of disk accesses required to retrieve mismapped mate pairs.
  - by default, the cache is reset after each locus, leading to a modest speedup with negligible memory overhead.
  - the new `--cache-mates` option activates reuse of the cache across loci, leading to a more significant speed increase, though at a cost of increased memory usage (typically in the range of 1-2GB of memory usage for catalogs with 100s to 1000s of loci). 
  - if/when spliting a large variant catalog into multiple shards, it's important to presort the loci by their normalized motif (which is the cyclic shift of a motif that is alphabetically first - ie. AGC rather than CAG). 
    This ensures that loci with the same normalized motif will be processed in the same shard, increasing cache hit rates and therefore speed for this optimization.
  

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
