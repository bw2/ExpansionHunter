# `chrparallel` — synthetic multi-contig EH test fixture

A tiny, fully deterministic ExpansionHunter (EH) fixture for validating that
`low-mem-streaming` and `optimized-streaming` produce **byte-identical** output
regardless of `--threads`, and that the per-contig / cross-contig / mate-cache
machinery handles a set of edge cases.

Regenerate with:

```
python3 ehunter/tests/data/chrparallel/make_fixture.py
```

It is seeded (`random.Random(20260619)`) with no timestamps, so every run emits
the same bytes.

## Files

| File | Description |
|------|-------------|
| `reference.fa` (+ `.fai`) | 5 contigs `chr1..chr5` (this exact order; EH derives contig index from FASTA/BAM header order), 4.2–5.5 kb each, random ACGT with embedded perfect repeat tracts. |
| `variant_catalog.json` | 4 loci, position-sortable. |
| `reads.bam` (+ `.bai`) | Coordinate-sorted, indexed; proper FR pairs (619 pairs / 1238 records). |
| `reads.cram` (+ `.crai`) | The same reads as CRAM against `reference.fa`. |
| `make_fixture.py` | The generator (kept in-tree so the data is reproducible / auditable). |

## Contig → locus → scenario

| Contig | Locus | Motif / structure | Read len | Scenario exercised |
|--------|-------|-------------------|----------|--------------------|
| chr1 | `CHR1_CAG` | `(CAG)*`, 15 copies | 100 bp | (f) mixed read length; (b) same-contig **far pair** (one end ~1500 bp upstream of the locus, mate distance ≥ 1000 bp) |
| chr2 | `CHR2_CCG` | `(CCG)*`, 12 copies | 150 bp | (f) mixed read length; (a) **cross-contig mate** (spanning reads on chr2 whose mate maps to chr5) |
| chr3 | `CHR3_GAA` | `(GAA)*`, 10 copies | — | (d) **zero coverage** — the locus is in the catalog but no reads cover it; emitted as a no-call |
| chr4 | `CHR4_MULTI` | `(A\|T)CCCC(ATTCT)*` (SmallVariant + Repeat) | 150 bp | (c) **multi-variant locus** — two variant entries in one locus, routed through full graph genotyping (the optimized fast path only handles single-repeat loci) |
| chr5 | *(none)* | — | 150 bp | (e) **no-loci contig that is the physical home of cross-contig mates** from (a) |

Scenario (c) note: the task asked for a locus whose *variantId string-order
differs from its positional order*. Variant IDs here are auto-derived from the
ReferenceRegion (`CHR4_MULTI_chr4:1400-1401`, `CHR4_MULTI_chr4:1405-1475`), so
for this locus string order happens to equal positional order. What the locus
**does** exercise — and the substantive point for cross-thread determinism — is
a multi-variant (SmallVariant + Repeat) locus flowing through full genotyping
with multiple variant records, plus the `unordered_map`-keyed per-locus cache
iteration in the streaming analyzers (the byte-identical-across-threads check
below is what proves that iteration order does not leak into the output).

## How the reads are built (two non-obvious EH requirements)

Both were discovered by iterating against the EH binary; they are documented in
`make_fixture.py` and matter for anyone editing this fixture:

1. **Reverse-strand SEQ is stored in forward-reference orientation.**
   pysam/htslib store `query_sequence` verbatim and convey strand only via the
   `is_reverse` flag. If you reverse-complement the bytes yourself for a
   reverse-strand read, the stored SEQ ends up double-reversed and EH's graph
   aligner rejects that mate, so the pair is never classified and the locus
   silently fails to genotype in full-genotyping (seeking / low-mem) modes.

2. **Coverage must clear `--min-locus-coverage` (default 10).**
   The full graph genotyper's `isLowDepth()` filter discards a locus whose
   estimated diploid depth is below the threshold and emits empty findings (no
   genotype). With the default `--region-extension-length` of 1000, the graph
   flank nodes are 1000 bp each, so `depth ≈ readLen × pairs / (2000 − readLen)`.
   The generator sizes each genotyped locus (`pairs_for_read_len`) to clear
   `depth ≈ 14`. The optimized-streaming fast path does **not** apply this
   filter, which is why under-covered loci can still genotype there.

## Expected results (verified)

Run from the repo root with the built binary
`build/ehunter-prefix/src/ehunter-build/ExpansionHunter`:

```
ExpansionHunter --reads   ehunter/tests/data/chrparallel/reads.bam \
                --reference ehunter/tests/data/chrparallel/reference.fa \
                --variant-catalog ehunter/tests/data/chrparallel/variant_catalog.json \
                --analysis-mode <low-mem-streaming|optimized-streaming> \
                --sort-catalog-by position --threads <N> \
                --output-prefix <prefix>
```

- exit code 0 in both modes; 4 `LocusResults` spanning **4 distinct contigs** (chr1–chr4).
- **Genotypes (identical in both modes):** `CHR1_CAG` = 15/15, `CHR2_CCG` = 12/12,
  `CHR4_MULTI` = SNV 0/0 + `(ATTCT)*` 14/14.
- `CHR3_GAA` is emitted as a **zero-coverage no-call** (`Coverage: 0.0`, no `Genotype`).
- `optimized-streaming` reports **2 (50.0%) loci via fast genotyping** (CHR1_CAG,
  CHR2_CCG), 1 via full genotyping (CHR4_MULTI), 1 zero-coverage no-call.
- `LocusResults` and `SampleParameters` are **byte-identical across `--threads 1`,
  `2`, and `4`** in both modes (the primary purpose of this fixture). The
  `RunInfo` record (`Started`/`Completed`/`Runtime`/`Threads`) is expected to
  differ, since it reports true wall-clock timing and the actual `--threads`
  value for each run.
- The **CRAM** input produces the same LocusIds and the same genotypes as the BAM.
