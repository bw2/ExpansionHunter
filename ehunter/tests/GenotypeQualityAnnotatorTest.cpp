//
// ExpansionHunter
//
// Unit tests for the per-allele genotype-quality prediction glue (genotypingRegime routing
// + feature assembly + model evaluation composed end to end).
//

#include "genotype_quality/GenotypeQualityAnnotator.hh"

#include <cmath>
#include <map>
#include <stdexcept>

#include "gmock/gmock.h"

using namespace ehunter;
using namespace ehunter::gq;

namespace
{

// q_median routes on feature 0 (motif_size) threshold 5 for the fast genotypingRegime; the
// full genotyping_regimes use an empty q head (LCF == 1) and differ only in their direction
// head, so the test can verify which genotypingRegime an allele was routed to.
// feature_names use real assembler feature names (consumed by name); the fast q tree
// routes on "motif_size" (canonical index 0), which the ctx below sets to 3.
const char* kModel = R"JSON(
{
  "format_version": 2,
  "feature_names": { "quick": ["motif_size"], "full": ["motif_size"] },
  "genotyping_regimes": {
    "quick": {
      "q_median": { "baseline": 0.1, "trees": [ { "nodes": [
        { "feature": 0, "threshold": 5.0, "missing_left": true, "left": 1, "right": 2 },
        { "leaf": true, "value": 0.2 },
        { "leaf": true, "value": -0.3 } ] } ] },
      "direction": { "classes": ["OK", "TOO_LONG", "TOO_SHORT"], "baseline": [0.0, 0.0, 0.0], "trees": [],
        "calibrators": [ {"x":[],"y":[]}, {"x":[],"y":[]}, {"x":[],"y":[]} ] }
    },
    "full_spanning": {
      "q_median": { "baseline": 0.0, "trees": [] },
      "direction": { "classes": ["OK", "TOO_LONG", "TOO_SHORT"], "baseline": [1.0, 0.0, 0.0], "trees": [],
        "calibrators": [ {"x":[],"y":[]}, {"x":[],"y":[]}, {"x":[],"y":[]} ] }
    },
    "full_nonspanning": {
      "q_median": { "baseline": 0.0, "trees": [] },
      "direction": { "classes": ["OK", "TOO_LONG", "TOO_SHORT"], "baseline": [0.0, 0.0, 0.0], "trees": [],
        "calibrators": [ {"x":[0.0,1.0],"y":[0.0,0.0]}, {"x":[0.0,1.0],"y":[0.0,0.0]}, {"x":[0.0,1.0],"y":[0.0,0.0]} ] }
    }
  }
}
)JSON";

AlleleMetrics makeAqm()
{
    AlleleMetrics aqm;
    aqm.depth = 30.0;
    return aqm;
}

GenotypeQualityModel model()
{
    return GenotypeQualityModel::fromJson(nlohmann::json::parse(kModel));
}

} // namespace

TEST(GenotypeQualityAnnotator, QuickGenotypingRegimeUsesQTree)
{
    const CountTable spanning;
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{/*motif=*/3, /*ref_bp=*/30, spanning, flanking, hq};

    // quickGenotype == true -> Quick genotypingRegime; motif_size 3 <= 5 -> leaf 0.2, t = 0.3.
    AllelePrediction p = predictAllele(model(), /*quick=*/true, ctx, 0, /*eh=*/20, 18, 24, makeAqm());
    EXPECT_NEAR(p.lengthCorrectionFactor, std::exp(0.3), 1e-12);
    EXPECT_NEAR(p.pOk, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pTooLong, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pTooShort, 1.0 / 3.0, 1e-12);
}

TEST(GenotypeQualityAnnotator, RoutesToFullSpanningWhenSpanningSupportsCall)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 3}}); // >=1 read at eh=20
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    AllelePrediction p = predictAllele(model(), /*quick=*/false, ctx, 0, /*eh=*/20, 18, 24, makeAqm());
    EXPECT_NEAR(p.lengthCorrectionFactor, 1.0, 1e-12); // empty q head -> exp(0)
    const double denom = std::exp(1.0) + 2.0;
    EXPECT_NEAR(p.pTooLong, 1.0 / denom, 1e-12);  // full_spanning baseline [1,0,0]
    EXPECT_NEAR(p.pTooShort, 1.0 / denom, 1e-12);
}

TEST(GenotypeQualityAnnotator, RoutesToFullNonspanningWithoutSpanningSupport)
{
    const CountTable spanning; // no read at eh -> spanning_at_called 0
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    AllelePrediction p = predictAllele(model(), /*quick=*/false, ctx, 0, /*eh=*/20, 18, 24, makeAqm());
    EXPECT_NEAR(p.lengthCorrectionFactor, 1.0, 1e-12);
    EXPECT_NEAR(p.pTooLong, 1.0 / 3.0, 1e-12);  // zero calibrators -> uniform fallback
    EXPECT_NEAR(p.pTooShort, 1.0 / 3.0, 1e-12);
}

TEST(GenotypeQualityAnnotator, RejectsModelFeatureTheAssemblerDoesNotProduce)
{
    // A model that lists a feature name EH does not assemble is rejected at predict time
    // rather than silently mis-evaluated (the by-name feature contract).
    const char* kBadModel = R"JSON(
    {
      "format_version": 2,
      "feature_names": { "quick": ["not_a_real_feature"], "full": ["not_a_real_feature"] },
      "genotyping_regimes": {
        "quick": { "q_median": { "baseline": 0.0, "trees": [] },
          "direction": { "classes": ["OK","TOO_LONG","TOO_SHORT"], "baseline": [0.0,0.0,0.0], "trees": [],
            "calibrators": [ {"x":[],"y":[]}, {"x":[],"y":[]}, {"x":[],"y":[]} ] } },
        "full_spanning": { "q_median": { "baseline": 0.0, "trees": [] },
          "direction": { "classes": ["OK","TOO_LONG","TOO_SHORT"], "baseline": [0.0,0.0,0.0], "trees": [],
            "calibrators": [ {"x":[],"y":[]}, {"x":[],"y":[]}, {"x":[],"y":[]} ] } },
        "full_nonspanning": { "q_median": { "baseline": 0.0, "trees": [] },
          "direction": { "classes": ["OK","TOO_LONG","TOO_SHORT"], "baseline": [0.0,0.0,0.0], "trees": [],
            "calibrators": [ {"x":[],"y":[]}, {"x":[],"y":[]}, {"x":[],"y":[]} ] } }
      }
    }
    )JSON";
    const GenotypeQualityModel bad = GenotypeQualityModel::fromJson(nlohmann::json::parse(kBadModel));

    const CountTable spanning;
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};
    EXPECT_THROW(predictAllele(bad, /*quick=*/true, ctx, 0, 20, 18, 24, makeAqm()), std::runtime_error);
}
