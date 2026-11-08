# `pseudoautosomal_region_fixtures`: synthetic fixture for `IsInsidePseudoautosomalRegion`

A tiny, deterministic ExpansionHunter fixture for checking that catalog loci marked
`"IsInsidePseudoautosomalRegion": true` are genotyped with two alleles in a male sample, and that nothing else
changes. Used by `ehunter/tests/pseudoautosomal_region_tests.py`.

Regenerate with:

```
python3 ehunter/tests/data/pseudoautosomal_region_fixtures/make_fixture.py
```

It is seeded (`random.Random(20261108)`), so every run writes the same bytes.

## Files

| File | Description |
|------|-------------|
| `reference.fa` (+ `.fai`) | 2 contigs: `chr1` (6 kb) and `chrX` (58 kb), random sequence with 12-unit CAG repeats. |
| `variant_catalog.json` | 14 `(CAG)*` loci (below); the 12 `PAR_CAG_*` loci have `"IsInsidePseudoautosomalRegion": true`. |
| `reads.bam` (+ `.bai`) | Coordinate-sorted and indexed, 1,400 read pairs, 150 bp reads, 400 bp fragments, mapping quality 60. |
| `make_fixture.py` | The generator. |

## Loci and reads

| Locus | Repeat (0-based, end exclusive) | Marked |
|-------|------|------|
| `PAR_CAG_01` ... `PAR_CAG_12` | `chrX:2000-2036`, `chrX:6000-6036`, ..., `chrX:46000-46036` (4 kb apart) | yes |
| `NONPAR_CAG` | `chrX:54000-54036` | no |
| `AUTOSOMAL_CAG` | `chr1:3000-3036` | no |

Every locus is heterozygous: 50 spanning read pairs from a 12-unit (reference length) allele and 50 from a 16-unit
allele, whose reads are aligned with a 12 bp insertion after the 6th repeat unit. Read 1 spans the repeat with at
least 20 bp of flank on each side; read 2 lies in the right flank.

## Expected results (verified)

- `--sex female`: `12/16` at all 14 loci, with or without the field.
- `--sex male`, catalog without the field: one allele at the 13 chrX loci, `12/16` at `AUTOSOMAL_CAG`.
- `--sex male`: `12/16` at the 12 marked loci in every analysis mode; `NONPAR_CAG` and `AUTOSOMAL_CAG` have the same
  JSON and VCF records as without the field.
