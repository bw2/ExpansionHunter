//
// ExpansionHunter
//
// Genotype-quality model: a per-allele quality predictor exported from the
// python `genotype_quality/` training harness. Given a feature vector it returns
// a PredictedLengthCorrectionFactor (= eh/true, the q-head median) and the direction
// probabilities pTooShort / pTooLong. The model is a small ensemble of
// gradient-boosted trees plus per-class isotonic calibrators, evaluated here in
// raw-threshold space so it reproduces sklearn's predictions without any heavy
// runtime dependency. See GENOTYPE_QUALITY_CPP_INTEGRATION_PLAN.md.
//

#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "thirdparty/json/json.hpp"

namespace ehunter
{
namespace gq
{

// One node of a gradient-boosted tree, evaluated in raw (un-binned) feature
// space. A leaf carries a value; an internal node routes on a single feature.
struct TreeNode
{
    bool isLeaf = false;
    double value = 0.0; // leaf output

    int feature = -1;       // internal: feature index into the genotypingRegime's feature vector
    double threshold = 0.0; // internal: x[feature] <= threshold goes left
    bool missingLeft = false; // internal: NaN x[feature] routes left when true
    int left = -1;
    int right = -1;
};

// A single regression/classification tree (root == nodes[0]).
struct Tree
{
    std::vector<TreeNode> nodes;

    // Returns the leaf value reached for feature vector x. NaN feature values
    // follow the node's missing-left flag (matching sklearn HistGradientBoosting).
    double eval(const std::vector<double>& x) const;
};

// The q head: a quantile (median) regressor predicting t = log(eh/true).
// PredictedLengthCorrectionFactor = exp(predictT).
struct QHead
{
    double baseline = 0.0;
    std::vector<Tree> trees;

    double predictT(const std::vector<double>& x) const;
};

// A fitted sklearn IsotonicRegression (out_of_bounds='clip', increasing). Empty
// knots == passthrough (used where a calib class was absent and python kept the
// raw probability unchanged).
struct Isotonic
{
    std::vector<double> x; // strictly ascending knots; empty => passthrough
    std::vector<double> y;

    double apply(double v) const;
};

// The direction head: a 3-class boosted classifier (one tree triple per boosting
// iteration) plus per-class isotonic calibrators. The model carries an explicit
// `classes` label array (e.g. ["OK","TOO_LONG","TOO_SHORT"]) that defines what each
// position of baseline/trees/calibrators means; posOk/posTooLong/posTooShort are
// resolved from those labels at load so the consumer never assumes a positional
// order. predictProba returns the calibrated, simplex-renormalized probabilities in
// the model's stored (label) order; predict() relabels them by the resolved indices.
struct DirectionHead
{
    std::array<double, 3> baseline = {{0.0, 0.0, 0.0}};
    std::vector<std::array<Tree, 3>> trees;
    std::array<Isotonic, 3> calibrators;

    // Positions of the OK / TOO_LONG / TOO_SHORT classes within the parallel arrays
    // above, resolved from the model's explicit `classes` labels.
    int posOk = 0;
    int posTooLong = 1;
    int posTooShort = 2;

    std::array<double, 3> predictProba(const std::vector<double>& x) const;
};

// The model set for one genotyping genotypingRegime (fast / full_spanning / full_nonspanning).
struct GenotypingRegimeModel
{
    QHead qMedian;
    DirectionHead direction;
};

// Per-allele prediction emitted into the JSON output.
struct AllelePrediction
{
    double lengthCorrectionFactor = 1.0; // eh/true
    double pOk = 1.0;
    double pTooShort = 0.0;
    double pTooLong = 0.0;
};

// The full genotype-quality model: per-genotypingRegime model sets plus the per-branch
// feature-name lists. The feature names are an explicit part of the contract: the
// annotator builds each allele's feature vector in the order named here, looking up
// each name in the C++ feature assembler (an unknown name is rejected), so the model
// and binary never rely on a shared positional convention. Tree `feature` indices
// reference positions in the matching list (quick genotypingRegime -> featureNamesQuick, full
// genotyping_regimes -> featureNamesFull); fromJson validates they are in range.
struct GenotypeQualityModel
{
    int formatVersion = 0;
    // Model identity (the model file's basename without extension), reported in the output JSON as
    // GenotypeQualityModelVersion. Set by the loader from the file name; not part of the model JSON.
    std::string version;
    std::vector<std::string> featureNamesQuick;
    std::vector<std::string> featureNamesFull;
    GenotypingRegimeModel quick;
    GenotypingRegimeModel fullSpanning;
    GenotypingRegimeModel fullNonspanning;

    // Runs both heads for the given genotypingRegime and feature vector.
    AllelePrediction predict(const GenotypingRegimeModel& genotypingRegime, const std::vector<double>& x) const;

    // Parses a model from its JSON representation (see the plan's schema). Throws
    // std::runtime_error on a malformed/unsupported document.
    static GenotypeQualityModel fromJson(const nlohmann::json& j);
};

// Loads a model from a .json or .json.gz file (gzip auto-detected). Throws
// std::runtime_error if the file cannot be read or parsed.
GenotypeQualityModel loadGenotypeQualityModel(const std::string& path);

// Loads the model compiled into the binary at build time from the gzip named by
// GQ_MODEL_FILE in ehunter/CMakeLists.txt. Returns nullptr when no model was
// embedded (the checked-in placeholder is empty until the trained model is dropped
// in), so the caller emits no quality fields by default. Throws on a corrupt blob.
std::shared_ptr<const GenotypeQualityModel> loadEmbeddedGenotypeQualityModel();

} // namespace gq
} // namespace ehunter
