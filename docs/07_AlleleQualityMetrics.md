# Allele Quality Metrics

ExpansionHunter computes per-allele quality metrics for repeat variants. These
metrics provide detailed information about the quality of each genotype call,
enabling downstream filtering and quality assessment.

## Enabling/Disabling

Quality metrics are **enabled by default**. To disable them:

```bash
ExpansionHunter --reads input.bam \
                --reference ref.fa \
                --catalog catalog.json \
                --output-prefix output \
                --dont-output-quality-metrics
```

## Output Structure

The `AlleleQualityMetrics` field appears within a repeat variant record in the
JSON output when the metrics could be computed for that variant. It is absent when
the read-realignment workflow that produces the metrics did not run for the locus
— for example loci that also contain a `SmallVariant`, no-call records, or runs with `--dont-output-quality-metrics`. For
heterozygous calls, there are two entries in the `Alleles` array; for homozygous
or hemizygous calls, there is one entry.

```json
"AlleleQualityMetrics": {
    "VariantId": "HTT",
    "Alleles": [
        {
            "AlleleNumber": 1,
            "AlleleSize": 19,
            "Depth": 15.0,
            "QD": 0.92,
            "MeanInsertedBasesWithinRepeats": 0.0,
            "MeanDeletedBasesWithinRepeats": 0.0,
            "ReadRepeatPurity": 1.0,
            "StrandBiasBinomialPhred": 0.0,
            "LeftFlankNormalizedDepth": 1.2,
            "RightFlankNormalizedDepth": 1.1,
            "HighQualityUnambiguousReads": 12,
            "ConfidenceIntervalDividedByAlleleSize": 0.0,
            "PredictedLengthCorrectionFactor": 1.0,
            "pOk": 0.94,
            "pTooShort": 0.04,
            "pTooLong": 0.02
        },
        {
            "AlleleNumber": 2,
            "AlleleSize": 45,
            "Depth": 12.0,
            "QD": 0.88,
            ...
        }
    ]
}
```



## Metric Definitions

### AlleleNumber

1-based allele index. For heterozygous calls:
- `1` = shorter allele (or only allele in homozygous/hemizygous cases)
- `2` = longer allele

### AlleleSize

Size of the allele in repeat units, matching the value in `Genotype`.

### Depth

Allele-specific read depth at the repeat: the number of read bases
matched to the repeat node divided by the allele's length in base pairs.

### QD (Quality by Depth)

The mean per-read alignment quality over the reads overlapping the repeat.

**Formula:** `QD = sum(readQuality) / numReadsOverlappingRepeat`

Where the denominator is the number of reads overlapping the repeat node, while `readQuality` for each read = `matchedBases / totalBases`

A value close to 1.0 indicates high-quality alignments with few mismatches,
insertions, or deletions. Lower values may indicate mapping issues or sequence
errors. `QD` is not emitted on quick-path (`QuickGenotype: true`) calls (see below).

### MeanInsertedBasesWithinRepeats

Average number of inserted bases within the repeat region per read.

**Formula:** `MeanInsertedBasesWithinRepeats = totalInsertedBasesInRepeat / numReadsOverlappingRepeat`

Elevated values may indicate:
- Somatic mosaicism
- PCR stutter artifacts
- Alignment errors

### MeanDeletedBasesWithinRepeats

Average number of deleted bases within the repeat node per read.

**Formula:** `MeanDeletedBasesWithinRepeats = totalDeletedBasesInRepeat / numReadsOverlappingRepeat`

Similar interpretation to `MeanInsertedBasesWithinRepeats`.

### ReadRepeatPurity

Fraction of the repeat-region read bases (the read bases falling within the repeat region, among the
reads supporting this allele) that match a perfect repeat sequence (eg. CAG.CAG.CAG.CAG) of the
catalog motif — a measure of how clean the observed repeat is, in `[0, 1]` (`1.0` = perfectly pure).

**Formula:** `ReadRepeatPurity = matchedRepeatRegionReadBases / totalRepeatRegionReadBases`

Each read contributes only the portion that overlaps the repeat region. This is derived from the
read-to-graph alignment (match vs mismatch/insertion against the repeat node), except in
optimized-streaming mode when `QuickGenotype` is true, where it is computed from spanning reads only.
Low values flag interrupted or compound repeats. The field is omitted when no repeat-region read bases
were observed for the allele. See also the per-locus `ReferenceRepeatPurity` (the same measure applied to
the reference repeat region), described in [Output JSON Files](05_OutputJsonFiles.md).

### StrandBiasBinomialPhred

Phred-scaled p-value from a binomial test for strand bias.

**Formula:** `StrandBiasBinomialPhred = -10 * log10(binomial_pvalue)`

Tests whether the forward/reverse strand distribution deviates significantly
from the expected 50/50 ratio. Higher values indicate greater strand bias.

| Value | Interpretation |
|-------|----------------|
| 0     | No bias (equal strand distribution) |
| < 10  | Minimal bias |
| 10-20 | Moderate bias |
| > 20  | Strong bias (potential mapping artifact) |

### LeftFlankNormalizedDepth and RightFlankNormalizedDepth

Coverage in the 20bp flanking regions normalized by allele depth.

**Formula:** `FlankNormalizedDepth = (flankCoveredBases / 20) / Depth`

Values close to 1.0 indicate consistent coverage across the locus. Very low
values may indicate:
- Reads not extending into the flanks
- Alignment truncation
- Potential false positive expansions

Very high values may indicate:
- Repetitive flanking sequence causing multi-mapping

### HighQualityUnambiguousReads

Count of reads that:
1. Have high alignment quality (match rate ≥ 0.9)
2. Unambiguously support only this allele (not consistent with other haplotypes)

For **homozygous** (e.g., 5/5) and **hemizygous** (e.g., male X chromosome) genotypes,
all reads are considered unambiguous by definition since there is only one distinct
allele they could support. The high-quality threshold still applies.

For **heterozygous** genotypes with different allele sizes, a read is unambiguous
only if it can be assigned to exactly one of the two alleles based on its alignment.

Low counts relative to `Depth` may indicate ambiguous read assignments (for
heterozygous calls) or lower-quality alignments.

### ConfidenceIntervalDividedByAlleleSize

Width of the confidence interval normalized by allele size.

**Formula:** `(CI_upper - CI_lower) / (AlleleSize + 1)`

Smaller values indicate more precise genotype estimates. Values close to 0
indicate high confidence; larger values indicate uncertainty.

| Value | Interpretation |
|-------|----------------|
| 0     | Exact genotype (typically for small alleles) |
| < 0.1 | High confidence |
| 0.1-0.5 | Moderate confidence |
| > 0.5 | Low confidence (large uncertainty relative to size) |

## Genotype-quality model predictions

In addition to the read-derived metrics above, ExpansionHunter now attaches
per-allele predictions from a genotype-quality model — a small
gradient-boosted-tree model (either the default one, or
supplied with `--genotype-quality-model`). The model takes the allele's quality
metrics and read profile as input and predicts (a) how far the call is likely to be
from the true allele size and (b) the direction of any error.

These four fields are added to each allele object:

| Field | Meaning |
|-------|---------|
| `PredictedLengthCorrectionFactor` | Estimated `called_size / true_size`. |
| `pOk` | Probability the call is close to the true size. |
| `pTooShort` | Probability the call under-estimates the true size. |
| `pTooLong` | Probability the call over-estimates the true size. |

`pOk + pTooShort + pTooLong` sum to 1 (up to rounding).

### PredictedLengthCorrectionFactor

The model's estimate of the called size divided by the true size (`called / true`):

- `≈ 1.0` — the call is consistent with the predicted true size.
- `< 1.0` — the call is likely **too short** (the true allele is larger).
- `> 1.0` — the call is likely **too long** (the true allele is smaller).

A size-corrected estimate can be computed as  `AlleleSize / PredictedLengthCorrectionFactor`.

### pOk, pTooShort, pTooLong

A calibrated 3-way probability over whether the call is within tolerance (`pOk`), an
under-call (`pTooShort`), or an over-call (`pTooLong`). `pOk` is the primary
confidence signal: high `pOk` means the call is probably already correct, while low
`pOk` flags a call the model expects to be off (and the two directional probabilities
say which way).

**Suggested use (apply gate).** The length correction is most useful exactly where the
model is unsure the call is right. A reasonable policy is to apply
`PredictedLengthCorrectionFactor` only when `pOk < 0.5`, leaving confident calls unchanged:

```python
def corrected_size(allele):
    eh = allele["AlleleSize"]
    if allele.get("pOk", 1.0) < 0.5 and "PredictedLengthCorrectionFactor" in allele:
        return eh / allele["PredictedLengthCorrectionFactor"]
    return eh
```

Applying the correction to high-`pOk` calls is generally counter-productive (it moves
already-correct calls), so gating on `pOk` is recommended rather than correcting every
allele.

## Quick-path (`QuickGenotype: true`) rows

In `--analysis-mode optimized-streaming` some loci are genotyped by a fast heuristic
instead of the full genotyper. Those variant records carry
`"QuickGenotype": true` (see [Output JSON files](05_OutputJsonFiles.md)). On these
rows the allele quality metrics differ from those computed using the full genotyper:

- **Omitted fields:** `QD`, `LeftFlankNormalizedDepth`, and
  `RightFlankNormalizedDepth` are not computed and are absent from each allele
  object. A consumer that uses these should treat their absence as "not available" rather than as a zero
  value.
- **Approximated fields:** `Depth`, `HighQualityUnambiguousReads`,
  `StrandBiasBinomialPhred`, `MeanInsertedBasesWithinRepeats`, and
  `MeanDeletedBasesWithinRepeats` are derived from the high-quality spanning-read
  counts rather than from full graph realignment. In particular,
  `HighQualityUnambiguousReads` is just the per-allele spanning-read depth on the
  quick path.

## Notes

- For homozygous calls, reads from both haplotypes are combined into a single
  allele's metrics
