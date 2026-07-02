//
// ExpansionHunter
//
// Implementation of the genotype-quality per-allele feature assembler.
//

#include "genotype_quality/GenotypeQualityFeatures.hh"

#include <cassert>
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
        "allele_rank", "ci_width", "ci_asymmetry", "ci_over_eh", "spanning_total",
        "hq_unamb_total", "spanning_at_called", "spanning_above_called", "flanking_above_called",
        "support_frac", "depth", "hq_unambiguous_reads", "strand_bias_phred",
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
    const int spanningAtCalled = ctx.spanningReads.countOf(eh);
    const int spanningAboveCalled = countsAbove(ctx.spanningReads, eh);
    const int flankingAboveCalled = countsAbove(ctx.flankingReads, eh);
    const double supportFrac = (spanningTotal > 0)
        ? (static_cast<double>(spanningAtCalled) / spanningTotal)
        : std::numeric_limits<double>::quiet_NaN();

    // -1.0 = not computed (JsonWriter's emit-omission sentinel); the model was trained on NaN
    // for missing values (eh_json.py's variant.get(...) -> None -> NaN), so map here.
    const double referenceRepeatPurity = (ctx.referenceRepeatPurity >= 0.0)
        ? ctx.referenceRepeatPurity
        : std::numeric_limits<double>::quiet_NaN();
    const double readRepeatPurity = (aqm.readRepeatPurity >= 0.0)
        ? aqm.readRepeatPurity
        : std::numeric_limits<double>::quiet_NaN();

    // features.FAST_FEATURES order.
    std::vector<double> features = {
        motifSize,
        numRepeatsInReference,
        refSizeBp,
        static_cast<double>(eh),
        ehMinusRef,
        static_cast<double>(alleleRank),
        ciWidth,
        ciAsymmetry,
        ciOverEh,
        static_cast<double>(spanningTotal),
        static_cast<double>(hqUnambTotal),
        static_cast<double>(spanningAtCalled),
        static_cast<double>(spanningAboveCalled),
        static_cast<double>(flankingAboveCalled),
        supportFrac,
        aqm.depth,
        static_cast<double>(aqm.highQualityUnambiguousReads),
        aqm.strandBiasBinomialPhred,
        aqm.meanInsertedBasesWithinRepeats,
        aqm.meanDeletedBasesWithinRepeats,
        referenceRepeatPurity,
        readRepeatPurity,
    };

    // The full genotyping_regimes append the two flank-normalized depths (FULL_FEATURES).
    if (genotypingRegime != GenotypingRegime::Quick)
    {
        features.push_back(aqm.leftFlankNormalizedDepth);
        features.push_back(aqm.rightFlankNormalizedDepth);
    }

    // The value order above must stay in lockstep with the canonical name list.
    assert(features.size() == featureNamesForGenotypingRegime(genotypingRegime).size());
    return features;
}

} // namespace gq
} // namespace ehunter
