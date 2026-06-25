//
// ExpansionHunter
//
// Unit tests for the genotype-quality model evaluator and JSON(.gz) loader.
//

#include "genotype_quality/GenotypeQualityModel.hh"

#include <cmath>
#include <fstream>
#include <limits>

#include <zlib.h>

#include "gmock/gmock.h"

using namespace ehunter;
using namespace ehunter::gq;

namespace
{

// A fixture model with one genotypingRegime model set per genotypingRegime. The q head adds a tree
// that routes on feature 0 (threshold 5, NaN goes left); the direction head is
// three single-leaf trees (raw logits [0,0,0] -> uniform) with passthrough
// calibrators, so the genotyping_regimes differ only where noted.
const char* kFixtureModel = R"JSON(
{
  "format_version": 2,
  "feature_names": {
    "quick": ["f0", "f1"],
    "full": ["f0", "f1", "f2", "f3"]
  },
  "genotyping_regimes": {
    "quick": {
      "q_median": {
        "baseline": 0.1,
        "trees": [
          { "nodes": [
            { "feature": 0, "threshold": 5.0, "missing_left": true, "left": 1, "right": 2 },
            { "leaf": true, "value": 0.2 },
            { "leaf": true, "value": -0.3 }
          ] }
        ]
      },
      "direction": {
        "classes": ["OK", "TOO_LONG", "TOO_SHORT"],
        "baseline": [0.0, 0.0, 0.0],
        "trees": [
          [
            { "nodes": [ { "leaf": true, "value": 0.0 } ] },
            { "nodes": [ { "leaf": true, "value": 0.0 } ] },
            { "nodes": [ { "leaf": true, "value": 0.0 } ] }
          ]
        ],
        "calibrators": [
          { "x": [], "y": [] },
          { "x": [], "y": [] },
          { "x": [], "y": [] }
        ]
      }
    },
    "full_spanning": {
      "q_median": { "baseline": 0.0, "trees": [] },
      "direction": {
        "classes": ["OK", "TOO_LONG", "TOO_SHORT"],
        "baseline": [1.0, 0.0, 0.0],
        "trees": [],
        "calibrators": [
          { "x": [], "y": [] },
          { "x": [], "y": [] },
          { "x": [], "y": [] }
        ]
      }
    },
    "full_nonspanning": {
      "q_median": { "baseline": 0.0, "trees": [] },
      "direction": {
        "classes": ["OK", "TOO_LONG", "TOO_SHORT"],
        "baseline": [0.0, 0.0, 0.0],
        "trees": [],
        "calibrators": [
          { "x": [0.0, 1.0], "y": [0.0, 0.0] },
          { "x": [0.0, 1.0], "y": [0.0, 0.0] },
          { "x": [0.0, 1.0], "y": [0.0, 0.0] }
        ]
      }
    }
  }
}
)JSON";

GenotypeQualityModel fixtureModel()
{
    return GenotypeQualityModel::fromJson(nlohmann::json::parse(kFixtureModel));
}

} // namespace

TEST(GenotypeQualityModel, TreeRoutesOnThreshold)
{
    Tree tree;
    tree.nodes = {
        TreeNode{}, // root, filled below
        TreeNode{}, TreeNode{},
    };
    tree.nodes[0] = TreeNode{false, 0.0, 0, 5.0, true, 1, 2};
    tree.nodes[1] = TreeNode{true, 0.2, -1, 0.0, false, -1, -1};
    tree.nodes[2] = TreeNode{true, -0.3, -1, 0.0, false, -1, -1};

    EXPECT_DOUBLE_EQ(tree.eval({3.0, 0.0}), 0.2);  // 3 <= 5 -> left
    EXPECT_DOUBLE_EQ(tree.eval({9.0, 0.0}), -0.3); // 9 > 5 -> right
    EXPECT_DOUBLE_EQ(tree.eval({5.0, 0.0}), 0.2);  // boundary is <= -> left
}

TEST(GenotypeQualityModel, TreeRoutesNaNToMissingChild)
{
    Tree tree;
    tree.nodes = {
        TreeNode{false, 0.0, 0, 5.0, true, 1, 2}, // NaN -> left
        TreeNode{true, 0.2, -1, 0.0, false, -1, -1},
        TreeNode{true, -0.3, -1, 0.0, false, -1, -1},
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_DOUBLE_EQ(tree.eval({nan, 0.0}), 0.2);

    tree.nodes[0].missingLeft = false; // NaN -> right
    EXPECT_DOUBLE_EQ(tree.eval({nan, 0.0}), -0.3);
}

TEST(GenotypeQualityModel, IsotonicClipsAndInterpolates)
{
    Isotonic iso;
    iso.x = {0.0, 1.0, 2.0};
    iso.y = {10.0, 20.0, 50.0};

    EXPECT_DOUBLE_EQ(iso.apply(-1.0), 10.0); // below first knot -> clip to y.front
    EXPECT_DOUBLE_EQ(iso.apply(3.0), 50.0);  // above last knot  -> clip to y.back
    EXPECT_DOUBLE_EQ(iso.apply(0.0), 10.0);
    EXPECT_DOUBLE_EQ(iso.apply(2.0), 50.0);
    EXPECT_DOUBLE_EQ(iso.apply(0.5), 15.0);  // halfway on [0,1] of [10,20]
    EXPECT_DOUBLE_EQ(iso.apply(1.5), 35.0);  // halfway on [1,2] of [20,50]
}

TEST(GenotypeQualityModel, IsotonicEmptyIsPassthrough)
{
    Isotonic iso; // no knots
    EXPECT_DOUBLE_EQ(iso.apply(0.42), 0.42);
    EXPECT_DOUBLE_EQ(iso.apply(-7.0), -7.0);
}

TEST(GenotypeQualityModel, QHeadLengthCorrectionFactor)
{
    const GenotypeQualityModel model = fixtureModel();

    // x[0] = 3 (<=5) -> leaf 0.2, t = 0.1 + 0.2 = 0.3
    AllelePrediction lo = model.predict(model.quick, {3.0, 0.0});
    EXPECT_NEAR(lo.lengthCorrectionFactor, std::exp(0.3), 1e-12);

    // x[0] = 9 (>5) -> leaf -0.3, t = 0.1 - 0.3 = -0.2
    AllelePrediction hi = model.predict(model.quick, {9.0, 0.0});
    EXPECT_NEAR(hi.lengthCorrectionFactor, std::exp(-0.2), 1e-12);
}

TEST(GenotypeQualityModel, DirectionUniformWhenLogitsEqual)
{
    const GenotypeQualityModel model = fixtureModel();
    AllelePrediction p = model.predict(model.quick, {3.0, 0.0});
    EXPECT_NEAR(p.pOk, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pTooLong, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pTooShort, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pOk + p.pTooLong + p.pTooShort, 1.0, 1e-12);
}

TEST(GenotypeQualityModel, DirectionSoftmaxOfBaselineLogits)
{
    const GenotypeQualityModel model = fixtureModel();
    // full_spanning baseline logits [1,0,0]; passthrough calibrators -> softmax.
    AllelePrediction p = model.predict(model.fullSpanning, {0.0, 0.0});
    const double denom = std::exp(1.0) + std::exp(0.0) + std::exp(0.0);
    EXPECT_NEAR(p.pOk, std::exp(1.0) / denom, 1e-12);      // class OK
    EXPECT_NEAR(p.pTooLong, std::exp(0.0) / denom, 1e-12);  // class TOO_LONG
    EXPECT_NEAR(p.pTooShort, std::exp(0.0) / denom, 1e-12); // class TOO_SHORT
}

TEST(GenotypeQualityModel, DirectionFallsBackToUniformWhenCalibratorsZero)
{
    const GenotypeQualityModel model = fixtureModel();
    // full_nonspanning calibrators map every probability to 0 -> renorm fallback.
    AllelePrediction p = model.predict(model.fullNonspanning, {0.0, 0.0});
    EXPECT_NEAR(p.pTooLong, 1.0 / 3.0, 1e-12);
    EXPECT_NEAR(p.pTooShort, 1.0 / 3.0, 1e-12);
}

TEST(GenotypeQualityModel, RejectsUnsupportedFormatVersion)
{
    nlohmann::json j = nlohmann::json::parse(kFixtureModel);
    j["format_version"] = 3;
    EXPECT_THROW(GenotypeQualityModel::fromJson(j), std::runtime_error);
}

TEST(GenotypeQualityModel, RejectsOutOfRangeFeatureIndex)
{
    // feature_names.quick has 2 entries, so a quick-genotypingRegime tree routing on feature 5 is invalid.
    nlohmann::json j = nlohmann::json::parse(kFixtureModel);
    j["genotyping_regimes"]["quick"]["q_median"]["trees"][0]["nodes"][0]["feature"] = 5;
    EXPECT_THROW(GenotypeQualityModel::fromJson(j), std::runtime_error);
}

TEST(GenotypeQualityModel, RejectsDirectionMissingClassLabel)
{
    nlohmann::json j = nlohmann::json::parse(kFixtureModel);
    j["genotyping_regimes"]["quick"]["direction"]["classes"] = {"OK", "TOO_LONG", "BOGUS"};
    EXPECT_THROW(GenotypeQualityModel::fromJson(j), std::runtime_error);
}

TEST(GenotypeQualityModel, DirectionRelabelsByClassOrderNotPosition)
{
    // Permute the class labels and the parallel arrays the same way: TOO_LONG first.
    // full_spanning's logits are [1,0,0] for [OK,TOO_LONG,TOO_SHORT]; moving the OK
    // column to position 1 with classes ["TOO_LONG","OK","TOO_SHORT"] must keep P_OK
    // mapped to the high-logit column.
    nlohmann::json j = nlohmann::json::parse(kFixtureModel);
    j["genotyping_regimes"]["full_spanning"]["direction"]["classes"] = {"TOO_LONG", "OK", "TOO_SHORT"};
    j["genotyping_regimes"]["full_spanning"]["direction"]["baseline"] = {0.0, 1.0, 0.0};
    GenotypeQualityModel model = GenotypeQualityModel::fromJson(j);
    AllelePrediction p = model.predict(model.fullSpanning, {0.0, 0.0});
    const double denom = std::exp(1.0) + 2.0;
    EXPECT_NEAR(p.pOk, std::exp(1.0) / denom, 1e-12);   // OK still maps to the high-logit column
    EXPECT_NEAR(p.pTooLong, 1.0 / denom, 1e-12);
    EXPECT_NEAR(p.pTooShort, 1.0 / denom, 1e-12);
}

TEST(GenotypeQualityModel, RejectsOutOfRangeTreeChild)
{
    nlohmann::json j = nlohmann::json::parse(kFixtureModel);
    j["genotyping_regimes"]["quick"]["q_median"]["trees"][0]["nodes"][0]["left"] = 99;
    EXPECT_THROW(GenotypeQualityModel::fromJson(j), std::runtime_error);
}

TEST(GenotypeQualityModel, LoadsFromPlainJsonFile)
{
    const std::string path = std::string(testing::TempDir()) + "/gq_model_plain.json";
    {
        std::ofstream out(path);
        out << kFixtureModel;
    }
    GenotypeQualityModel model = loadGenotypeQualityModel(path);
    EXPECT_EQ(model.formatVersion, 2);
    EXPECT_EQ(model.featureNamesQuick.size(), 2u);
    EXPECT_EQ(model.featureNamesFull.size(), 4u);
}

TEST(GenotypeQualityModel, LoadsFromGzippedJsonFile)
{
    const std::string path = std::string(testing::TempDir()) + "/gq_model.json.gz";
    {
        gzFile out = gzopen(path.c_str(), "wb");
        ASSERT_NE(out, nullptr);
        gzwrite(out, kFixtureModel, static_cast<unsigned>(std::char_traits<char>::length(kFixtureModel)));
        gzclose(out);
    }
    GenotypeQualityModel model = loadGenotypeQualityModel(path);
    EXPECT_EQ(model.formatVersion, 2);
    AllelePrediction p = model.predict(model.quick, {3.0, 0.0});
    EXPECT_NEAR(p.lengthCorrectionFactor, std::exp(0.3), 1e-12);
}

TEST(GenotypeQualityModel, ThrowsOnMissingFile)
{
    EXPECT_THROW(loadGenotypeQualityModel("/nonexistent/path/model.json.gz"), std::runtime_error);
}
