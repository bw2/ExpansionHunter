//
// ExpansionHunter
//
// Implementation of the genotype-quality per-allele feature assembler.
//

#include "genotype_quality/GenotypeQualityFeatures.hh"

#include <cassert>
#include <cmath>
#include <limits>

namespace ehunter
{
namespace gq
{

const std::vector<std::string>& featureNamesForGenotypingRegime(GenotypingRegime genotypingRegime)
{
    // features.QUICK_FEATURES order (kept in lockstep with assembleFeatures below).
    static const std::vector<std::string> quick = {
        "motif_size", "num_repeats_in_reference", "ref_size_bp", "eh", "eh_minus_ref",
        "allele_rank", "n_alleles", "n_distinct_alleles", "ci_width", "ci_asymmetry", "ci_over_eh",
        "spanning_total", "hq_unamb_total", "flanking_total",
        "spanning_at_called", "spanning_above_called", "flanking_above_called",
        "support_frac", "flanking_frac", "coverage",
        "depth", "hq_unambiguous_reads", "strand_bias_phred",
        "mean_inserted_bases", "mean_deleted_bases",
        "reference_repeat_purity", "read_repeat_purity"};
    // features.FULL_FEATURES = QUICK_FEATURES + the two flank-normalized depths.
    static const std::vector<std::string> full = [] {
        std::vector<std::string> names = quick;
        names.push_back("left_flank_norm_depth");
        names.push_back("right_flank_norm_depth");
        return names;
    }();
    return genotypingRegime == GenotypingRegime::Quick ? quick : full;
}

namespace
{

// Mirrors JsonWriter's round3: the training parquets are parsed from the emitted JSON, and every
// AlleleQualityMetrics-derived field is written there through that same 3-decimal rounding. Feeding
// the model the unrounded double would differ by up to 5e-4 from what training saw. NaN rounds to NaN.
double roundLikeJson3(double value)
{
    return std::round(value * 1000.0) / 1000.0;
}

int sumCounts(const CountTable& table)
{
    int total = 0;
    for (const auto& entry : table)
    {
        total += entry.second;
    }
    return total;
}

int countsAbove(const CountTable& table, int threshold)
{
    int total = 0;
    for (const auto& entry : table)
    {
        if (entry.first > threshold)
        {
            total += entry.second;
        }
    }
    return total;
}

} // namespace

int numDistinctAllelesOf(bool isHomozygous, int numAlleles)
{
    return isHomozygous ? 1 : numAlleles;
}

GenotypingRegime genotypingRegimeOf(bool quickGenotype, int spanningAtCalled)
{
    if (quickGenotype)
    {
        return GenotypingRegime::Quick;
    }
    return spanningAtCalled >= 1 ? GenotypingRegime::FullSpanning : GenotypingRegime::FullNonspanning;
}

std::vector<double> assembleFeatures(
    const LocusFeatureContext& ctx, int alleleRank, int eh, int ciStart, int ciEnd,
    const AlleleMetrics& aqm, GenotypingRegime genotypingRegime)
{
    const double motifSize = ctx.motifSize;
    const double refSizeBp = ctx.refSizeBp;
    const double numRepeatsInReference = (motifSize > 0) ? (refSizeBp / motifSize) : std::numeric_limits<double>::quiet_NaN();
    const double ehMinusRef = eh - numRepeatsInReference;

    const double ciWidth = ciEnd - ciStart;
    // Engineered CI columns (features.add_engineered); inputs are always present
    // here, so the computed value is always kept (never the missing -> 0 branch).
    const double ciAsymmetry = ((ciEnd - eh) - (eh - ciStart)) / (ciWidth + 1.0);
    const double ciOverEh = ciWidth / (eh + 1.0);

    const int spanningTotal = sumCounts(ctx.spanningReads);
    const int hqUnambTotal = sumCounts(ctx.hqUnambiguousReads);
    const int flankingTotal = sumCounts(ctx.flankingReads);
    const int spanningAtCalled = ctx.spanningReads.countOf(eh);
    const int spanningAboveCalled = countsAbove(ctx.spanningReads, eh);
    const int flankingAboveCalled = countsAbove(ctx.flankingReads, eh);
    const double supportFrac = (spanningTotal > 0)
        ? (static_cast<double>(spanningAtCalled) / spanningTotal)
        : std::numeric_limits<double>::quiet_NaN();
    // Flanking share of the locus's informative reads. 0 (not NaN) when there are no flanking
    // reads at all, which also covers the 0/0 case -- matches eh_json.py's flanking_frac.
    const double flankingFrac = (flankingTotal > 0)
        ? (static_cast<double>(flankingTotal) / (flankingTotal + spanningTotal))
        : 0.0;

    // 0 = not supplied (see LocusFeatureContext). Production always has a genotype in hand here, so
    // this branch is unit-test-only; NaN (not 0) is used so a tree takes its learned missing-value
    // branch instead of reading the sentinel as a real allele count.
    const double numAlleles = (ctx.numAlleles > 0)
        ? static_cast<double>(ctx.numAlleles)
        : std::numeric_limits<double>::quiet_NaN();
    const double numDistinctAlleles = (ctx.numDistinctAlleles > 0)
        ? static_cast<double>(ctx.numDistinctAlleles)
        : std::numeric_limits<double>::quiet_NaN();

    // -1.0 = not computed (JsonWriter's emit-omission sentinel); the model was trained on NaN
    // for missing values (eh_json.py's variant.get(...) -> None -> NaN), so map here.
    const double referenceRepeatPurity = (ctx.referenceRepeatPurity >= 0.0)
        ? roundLikeJson3(ctx.referenceRepeatPurity)
        : std::numeric_limits<double>::quiet_NaN();
    const double readRepeatPurity = (aqm.readRepeatPurity >= 0.0)
        ? roundLikeJson3(aqm.readRepeatPurity)
        : std::numeric_limits<double>::quiet_NaN();

    // features.QUICK_FEATURES order.
    std::vector<double> features = {
        motifSize,
        numRepeatsInReference,
        refSizeBp,
        static_cast<double>(eh),
        ehMinusRef,
        static_cast<double>(alleleRank),
        numAlleles,
        numDistinctAlleles,
        ciWidth,
        ciAsymmetry,
        ciOverEh,
        static_cast<double>(spanningTotal),
        static_cast<double>(hqUnambTotal),
        static_cast<double>(flankingTotal),
        static_cast<double>(spanningAtCalled),
        static_cast<double>(spanningAboveCalled),
        static_cast<double>(flankingAboveCalled),
        supportFrac,
        flankingFrac,
        ctx.coverage,  // already rounded by the caller, to match the emitted LocusResults.Coverage
        roundLikeJson3(aqm.depth),
        static_cast<double>(aqm.highQualityUnambiguousReads),
        roundLikeJson3(aqm.strandBiasBinomialPhred),
        roundLikeJson3(aqm.meanInsertedBasesWithinRepeats),
        roundLikeJson3(aqm.meanDeletedBasesWithinRepeats),
        referenceRepeatPurity,
        readRepeatPurity,
    };

    // The full genotyping_regimes append the two flank-normalized depths (FULL_FEATURES).
    if (genotypingRegime != GenotypingRegime::Quick)
    {
        features.push_back(roundLikeJson3(aqm.leftFlankNormalizedDepth));
        features.push_back(roundLikeJson3(aqm.rightFlankNormalizedDepth));
    }

    // NOTE: values are left in double precision on purpose, and the model must be trained the same way.
    // An earlier pipeline downcast every float64 parquet column to float32 before fitting, so thresholds
    // were learned on float32 values and this function had to reproduce that downcast; otherwise a value
    // not exactly representable in float32 (coverage 47.47, flanking_frac 21/53) reached the model as a
    // different number than it did in training, and any threshold falling in the ~1 ULP gap between the
    // two routed the allele to the wrong child. The training side now keeps float64 end to end, so the
    // downcast is gone from both. Reinstating it on either side alone silently re-opens that gap.
    //
    // The `roundLikeJson3` / pre-rounded-coverage calls above are a separate matter and still required:
    // the training parquets are parsed from the JSON this binary emits, so the model only ever saw the
    // 3-decimal (2 for coverage) emitted values.

    // The value order above must stay in lockstep with the canonical name list.
    assert(features.size() == featureNamesForGenotypingRegime(genotypingRegime).size());
    return features;
}

} // namespace gq
} // namespace ehunter
