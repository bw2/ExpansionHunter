# `offtarget_regions_fixtures`: synthetic fixture with off-target read pairs

A tiny, deterministic ExpansionHunter fixture for checking that `optimized-streaming` mode genotypes a
locus that lists `OfftargetRegions` exactly as `seeking` mode does. Used by
`ehunter/tests/offtarget_regions_tests.py`.

Regenerate with:

```
python3 ehunter/tests/data/offtarget_regions_fixtures/make_fixture.py
```

It is seeded (`random.Random(20261015)`), so every run writes the same bytes.

## Files

| File | Description |
|------|-------------|
| `reference.fa` (+ `.fai`) | 3 contigs, `chr1..chr3`, random sequence with embedded repeat sequences. |
| `variant_catalog.json` | 2 loci: `C9LIKE` (with `OfftargetRegions`) and `CHR2_CAG` (without). |
| `reads.bam` (+ `.bai`) | Coordinate-sorted and indexed, 495 read pairs, 150 bp reads, 400 bp fragments. |
| `make_fixture.py` | The generator. |

## Loci and reads

| Contig | Content |
|--------|---------|
| chr1 | `C9LIKE`: `(GGCCCC)*`, 5 repeat units in the reference, `VariantType` `RareRepeat`, `OfftargetRegions` `chr2:1000-1500` and `chr3:1000-1500`. |
| chr2 | A 6-unit `GGCCCC` sequence inside the first off-target region, and `CHR2_CAG`: `(CAG)*`, 15 units, no off-target regions. |
| chr3 | A 6-unit `GGGGCC` sequence (the reverse complement) inside the second off-target region. No loci. |

`C9LIKE` is heterozygous: one allele has the reference length (90 spanning read pairs) and the other has 200
repeat units (1200 bp), tiled with fragments every 12 bp:

- reads in a flank, or crossing an edge of the repeat, are mapped on chr1 (soft-clipped past the reference repeat);
- a read inside the repeat whose mate is in a flank is unmapped and placed at its mate, as BWA does;
- read pairs with both reads inside the repeat are mapped (MAPQ 0) in the off-target regions: 27 pairs with both
  reads on chr2, 27 with both reads on chr3 (on the opposite strand, since chr3 holds the reverse complement),
  and 13 split between chr2 and chr3 (far-apart mates).

`CHR2_CAG` has 180 spanning read pairs of the reference allele.

## Expected results (verified)

- `seeking` calls `C9LIKE` 5/101. With `OfftargetRegions` removed from the catalog it calls 5/42, so the
  off-target read pairs change the call.
- `optimized-streaming` writes the same `C9LIKE` JSON record and VCF record as `seeking` with every
  `--genotyping-approach` (`auto`, `only-full`, `only-quick`) and every `--threads` value.
- `CHR2_CAG` is called 15/15 in every mode.
- With `--output-motif-composition all-loci`, each of the 134 reads in the off-target pairs adds to the number of
  reads with `GGCCCC` at `C9LIKE` (313 with the off-target regions, 179 without them).
