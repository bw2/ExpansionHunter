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

The `AlleleQualityMetrics` field appears within each repeat variant record in the
JSON output. For heterozygous calls, there are two entries in the `Alleles` array;
for homozygous or hemizygous calls, there is one entry.

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
            "StrandBiasBinomialPhred": 0.0,
            "LeftFlankNormalizedDepth": 1.2,
            "RightFlankNormalizedDepth": 1.1,
            "HighQualityUnambiguousReads": 12,
            "ConfidenceIntervalDividedByAlleleSize": 0.0
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

## Fast-path (`QuickGenotype`) rows

In `--analysis-mode optimized-streaming` some loci are genotyped by a fast heuristic
(`processLocusFast`) instead of the full genotyper. Those variant records carry
`"QuickGenotype": true` (see [Output JSON files](05_OutputJsonFiles.md)). On these
rows the allele quality metrics differ from the full genotyper, so the set of
fields below is **not** identical to a full-genotyper call:

- **Omitted fields:** `QD`, `LeftFlankNormalizedDepth`, and
  `RightFlankNormalizedDepth` are not computed and are absent from each allele
  object. A consumer that uses these should treat their absence on a
  `QuickGenotype` row as "not available" rather than as a failing (e.g. zero)
  value.
- **Approximated fields:** `Depth`, `HighQualityUnambiguousReads`,
  `StrandBiasBinomialPhred`, `MeanInsertedBasesWithinRepeats`, and
  `MeanDeletedBasesWithinRepeats` are derived from the high-quality spanning-read
  counts rather than from full graph realignment. In particular,
  `HighQualityUnambiguousReads` is just the per-allele spanning-read depth on the
  fast path.

Full-genotyper calls (those where the `QuickGenotype: true` field is absent)
always carry the complete set of fields described below.

## Metric Definitions

### AlleleNumber

1-based allele index. For heterozygous calls:
- `1` = shorter allele (or only allele in homozygous/hemizygous cases)
- `2` = longer allele

### AlleleSize

Size of the allele in repeat units, matching the value in `Genotype`.

### Depth

Allele-specific read depth: the count of reads supporting this allele. This is
derived from the same read assignment used for genotyping.

### QD (Quality by Depth)

The sum of read quality scores divided by allele depth.

**Formula:** `QD = sum(readQuality) / Depth`

Where `readQuality` for each read is computed as:
`readQuality = matchedBases / totalBases`

A value close to 1.0 indicates high-quality alignments with few mismatches,
insertions, or deletions. Lower values may indicate mapping issues or sequence
errors.

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

## Example Filtering Strategies

### Basic quality filter

```python
# Require reasonable depth and quality
if metrics["Depth"] >= 5 and metrics["QD"] >= 0.8:
    # Pass filter
```

### Strand bias filter

```python
# Filter out calls with extreme strand bias
if metrics["StrandBiasBinomialPhred"] < 30:
    # Pass filter
```

### Comprehensive filter

```python
def passes_quality_filter(allele_metrics):
    return (
        allele_metrics["Depth"] >= 5 and
        allele_metrics["QD"] >= 0.75 and
        allele_metrics["StrandBiasBinomialPhred"] < 30 and
        allele_metrics["HighQualityUnambiguousReads"] >= 3 and
        allele_metrics["LeftFlankNormalizedDepth"] >= 0.3 and
        allele_metrics["RightFlankNormalizedDepth"] >= 0.3
    )
```

## Notes

- Metrics are computed using the same read-to-haplotype assignments used for
  genotyping and read visualization
- For homozygous calls, reads from both haplotypes are combined into a single
  allele's metrics
- Values are rounded to 3 decimal places in the JSON output
