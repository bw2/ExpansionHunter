## Note about this optimized version of ExpansionHunter 

This fork of ExpansionHunter introduces read caching to speed up ExpansionHunter's default **seeking** mode by 2-3x without changing the output. 

To take advantage of these optimizations:
```
1) Install this optimized version
2) Run it with the new --cache-mates option 
3) If you are splitting your variant catalog into multiple shards, it's important to presort the catalog
   by the normalized motif so that loci with the same motif are grouped into the same shard. 
   This significantly improves performance by increasing cache hit rates. 
```

NOTE: if you want to process a very large catalog and have a machine with a lot of memory, you may find that it's faster / cheaper to use ExpansionHunter's **streaming** mode than this optimized **seeking** mode.

This fork includes several additional modifications:

It changes the `Flanks can contain at most 5 characters N but found x Ns` error to a warning that ExpansionHunter will print before continuing to the next locus. This allows ExpansionHunter to run to completion without exiting on these loci and makes it easier to process large catalogs without having to find and exclude these loci first.

Also, it adds an option to generate a table with per-locus runtimes. Since ExpansionHunter takes longer to process some loci, this makes it easier to find the small subset of outliers and exclude them from a large variant catalog:
```
--record-timing  write out a .tsv file with information on how long each locus takes to process. Processing time can be longer for some loci compared to others, and this allows the few slowest outliers to be excluded from a large variant catalog. 
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
[PolyForm Strict License 1.0.0](LICENSE.txt). It relies on several third party
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
