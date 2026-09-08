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
#include <string>
#include <vector>

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

// -- Readiness for a model retrained on the full assembled feature set -------------------------------
//
// The shipped model declares a SUBSET of what the assembler produces (it names 22 of the 27 quick
// features), so nothing above exercises the five that are assembled but currently unconsumed:
// n_alleles, n_distinct_alleles, flanking_total, flanking_frac and coverage. A retrain that adds them
// must not be the first thing to discover whether the binary can serve them. These tests declare the
// FULL canonical list and route a tree on each of the five, so a name-resolution or value-plumbing
// break shows up here rather than after a multi-hour retrain.

namespace
{

// A model declaring every canonical feature name, with one quick q_median tree that routes on
// `featureIndex` at `threshold`: <= goes to leaf +0.2, > goes to leaf -0.3. Names are taken from the
// assembler itself, so this cannot drift out of sync with it.
GenotypeQualityModel modelSplittingOn(int featureIndex, double threshold)
{
    nlohmann::json j;
    j["format_version"] = 2;
    j["feature_names"]["quick"] = featureNamesForGenotypingRegime(GenotypingRegime::Quick);
    j["feature_names"]["full"] = featureNamesForGenotypingRegime(GenotypingRegime::FullSpanning);

    // Built element by element rather than with brace initializer lists: nlohmann cannot tell an
    // object literal from an array of pairs, and silently picks the wrong one for these shapes.
    nlohmann::json calibrators = nlohmann::json::array();
    for (int i = 0; i < 3; ++i)
    {
        nlohmann::json calibrator = nlohmann::json::object();
        calibrator["x"] = nlohmann::json::array();
        calibrator["y"] = nlohmann::json::array();
        calibrators.push_back(calibrator);
    }
    nlohmann::json emptyDirection = nlohmann::json::object();
    emptyDirection["classes"] = std::vector<std::string>{ "OK", "TOO_LONG", "TOO_SHORT" };
    emptyDirection["baseline"] = std::vector<double>{ 0.0, 0.0, 0.0 };
    emptyDirection["trees"] = nlohmann::json::array();
    emptyDirection["calibrators"] = calibrators;

    nlohmann::json split = nlohmann::json::object();
    split["feature"] = featureIndex;
    split["threshold"] = threshold;
    split["missing_left"] = true;
    split["left"] = 1;
    split["right"] = 2;
    nlohmann::json leftLeaf = nlohmann::json::object();
    leftLeaf["leaf"] = true;
    leftLeaf["value"] = 0.2;
    nlohmann::json rightLeaf = nlohmann::json::object();
    rightLeaf["leaf"] = true;
    rightLeaf["value"] = -0.3;
    nlohmann::json nodes = nlohmann::json::array();
    nodes.push_back(split);
    nodes.push_back(leftLeaf);
    nodes.push_back(rightLeaf);
    nlohmann::json tree = nlohmann::json::object();
    tree["nodes"] = nodes;

    j["genotyping_regimes"]["quick"]["q_median"]["baseline"] = 0.0;
    j["genotyping_regimes"]["quick"]["q_median"]["trees"] = nlohmann::json::array({ tree });
    j["genotyping_regimes"]["quick"]["direction"] = emptyDirection;
    for (const char* regime : { "full_spanning", "full_nonspanning" })
    {
        j["genotyping_regimes"][regime]["q_median"]["baseline"] = 0.0;
        j["genotyping_regimes"][regime]["q_median"]["trees"] = nlohmann::json::array();
        j["genotyping_regimes"][regime]["direction"] = emptyDirection;
    }
    return GenotypeQualityModel::fromJson(j);
}

int canonicalIndexOf(const std::string& name)
{
    const std::vector<std::string>& names = featureNamesForGenotypingRegime(GenotypingRegime::Quick);
    for (std::size_t i = 0; i < names.size(); ++i)
    {
        if (names[i] == name)
        {
            return static_cast<int>(i);
        }
    }
    return -1;
}

} // namespace

TEST(GenotypeQualityAnnotator, ModelMayDeclareEveryAssembledFeature)
{
    // Declaring the full list must load and evaluate; nothing may be assembled-but-unresolvable.
    const GenotypeQualityModel full = modelSplittingOn(0, 5.0);
    EXPECT_EQ(full.featureNamesQuick.size(), 27u);
    EXPECT_EQ(full.featureNamesFull.size(), 29u);

    CountTable spanning;
    spanning.setCountOf(20, 4);
    CountTable flanking;
    flanking.setCountOf(22, 6);
    const CountTable hq;
    const LocusFeatureContext ctx{ 3, 30, spanning, flanking, hq, -1.0, 47.47, 2, 2 };
    EXPECT_NO_THROW(predictAllele(full, /*quick=*/true, ctx, 0, /*eh=*/20, 18, 24, makeAqm()));
}

TEST(GenotypeQualityAnnotator, CurrentlyUnconsumedFeaturesReachTheModelWithTheRightValue)
{
    CountTable spanning;
    spanning.setCountOf(20, 4); // spanning_total = 4
    CountTable flanking;
    flanking.setCountOf(22, 6); // flanking_total = 6, so flanking_frac = 6/10 = 0.6
    const CountTable hq;
    const LocusFeatureContext ctx{ 3, 30, spanning, flanking, hq, -1.0, /*coverage=*/47.47,
        /*numAlleles=*/2, /*numDistinctAlleles=*/2 };

    // For each of the five, a threshold just below the expected value must route right (leaf -0.3) and
    // one just above must route left (leaf +0.2). That pins the value, not merely the name lookup.
    const std::map<std::string, double> expectedValue{ { "n_alleles", 2.0 }, { "n_distinct_alleles", 2.0 },
        { "flanking_total", 6.0 }, { "flanking_frac", 0.6 }, { "coverage", 47.47 } };

    for (const auto& nameAndValue : expectedValue)
    {
        const int index = canonicalIndexOf(nameAndValue.first);
        ASSERT_NE(index, -1) << nameAndValue.first << " is not an assembled feature";
        const double value = nameAndValue.second;

        const AllelePrediction below
            = predictAllele(modelSplittingOn(index, value - 0.5), true, ctx, 0, 20, 18, 24, makeAqm());
        EXPECT_NEAR(below.lengthCorrectionFactor, std::exp(-0.3), 1e-12)
            << nameAndValue.first << " did not route right of " << (value - 0.5);

        const AllelePrediction above
            = predictAllele(modelSplittingOn(index, value + 0.5), true, ctx, 0, 20, 18, 24, makeAqm());
        EXPECT_NEAR(above.lengthCorrectionFactor, std::exp(0.2), 1e-12)
            << nameAndValue.first << " did not route left of " << (value + 0.5);
    }
}
