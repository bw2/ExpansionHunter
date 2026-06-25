//
// ExpansionHunter
//
// Implementation of the genotype-quality model evaluator and JSON(.gz) loader.
//

#include "genotype_quality/GenotypeQualityModel.hh"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <zlib.h>

namespace ehunter
{
namespace gq
{

double Tree::eval(const std::vector<double>& x) const
{
    int n = 0;
    // Bounded by node count; a leaf terminates. Malformed (cyclic) trees are
    // rejected at load time, so this cannot loop forever on a valid model.
    while (!nodes[n].isLeaf)
    {
        const TreeNode& node = nodes[n];
        const double xv = x[node.feature];
        if (std::isnan(xv))
        {
            n = node.missingLeft ? node.left : node.right;
        }
        else
        {
            n = (xv <= node.threshold) ? node.left : node.right;
        }
    }
    return nodes[n].value;
}

double QHead::predictT(const std::vector<double>& x) const
{
    double t = baseline;
    for (const Tree& tree : trees)
    {
        t += tree.eval(x);
    }
    return t;
}

double Isotonic::apply(double v) const
{
    if (x.empty())
    {
        return v; // passthrough calibrator
    }
    if (v <= x.front())
    {
        return y.front();
    }
    if (v >= x.back())
    {
        return y.back();
    }
    // x is strictly ascending; find the bracketing interval [lo, hi].
    const std::size_t hi = static_cast<std::size_t>(
        std::upper_bound(x.begin(), x.end(), v) - x.begin());
    const std::size_t lo = hi - 1;
    const double span = x[hi] - x[lo];
    if (span <= 0.0)
    {
        return y[lo];
    }
    const double frac = (v - x[lo]) / span;
    return y[lo] + frac * (y[hi] - y[lo]);
}

std::array<double, 3> DirectionHead::predictProba(const std::vector<double>& x) const
{
    std::array<double, 3> raw = baseline;
    for (const std::array<Tree, 3>& triple : trees)
    {
        for (int c = 0; c < 3; ++c)
        {
            raw[c] += triple[c].eval(x);
        }
    }

    // Softmax (shift by max for numerical stability).
    const double m = std::max({raw[0], raw[1], raw[2]});
    std::array<double, 3> p;
    double sum = 0.0;
    for (int c = 0; c < 3; ++c)
    {
        p[c] = std::exp(raw[c] - m);
        sum += p[c];
    }
    for (int c = 0; c < 3; ++c)
    {
        p[c] /= sum;
    }

    // Per-class isotonic calibration, then renormalize onto the simplex (matching
    // model_direction.predict_proba). All-zero calibrated rows fall back to uniform.
    std::array<double, 3> cal;
    double calSum = 0.0;
    for (int c = 0; c < 3; ++c)
    {
        cal[c] = std::min(1.0, std::max(0.0, calibrators[c].apply(p[c])));
        calSum += cal[c];
    }
    std::array<double, 3> out;
    if (calSum > 0.0)
    {
        for (int c = 0; c < 3; ++c)
        {
            out[c] = cal[c] / calSum;
        }
    }
    else
    {
        out = {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}};
    }
    return out;
}

AllelePrediction GenotypeQualityModel::predict(const GenotypingRegimeModel& genotypingRegime, const std::vector<double>& x) const
{
    AllelePrediction prediction;
    prediction.lengthCorrectionFactor = std::exp(genotypingRegime.qMedian.predictT(x));
    const std::array<double, 3> proba = genotypingRegime.direction.predictProba(x);
    // proba is in the model's stored class order; relabel by the resolved positions.
    prediction.pOk = proba[genotypingRegime.direction.posOk];
    prediction.pTooLong = proba[genotypingRegime.direction.posTooLong];
    prediction.pTooShort = proba[genotypingRegime.direction.posTooShort];
    return prediction;
}

namespace
{

Tree parseTree(const nlohmann::json& jt)
{
    Tree tree;
    const nlohmann::json& jnodes = jt.at("nodes");
    tree.nodes.reserve(jnodes.size());
    for (const nlohmann::json& jn : jnodes)
    {
        TreeNode node;
        if (jn.value("leaf", false))
        {
            node.isLeaf = true;
            node.value = jn.at("value").get<double>();
        }
        else
        {
            node.isLeaf = false;
            node.feature = jn.at("feature").get<int>();
            node.threshold = jn.at("threshold").get<double>();
            node.missingLeft = jn.at("missing_left").get<bool>();
            node.left = jn.at("left").get<int>();
            node.right = jn.at("right").get<int>();
        }
        tree.nodes.push_back(node);
    }
    if (tree.nodes.empty())
    {
        throw std::runtime_error("genotype-quality model: tree with no nodes");
    }
    // Validate child indices and acyclicity (children strictly greater than parent
    // is not required by sklearn, so just bound-check and reject self/empty cycles
    // by ensuring every internal node points to in-range distinct children).
    const int count = static_cast<int>(tree.nodes.size());
    for (int i = 0; i < count; ++i)
    {
        const TreeNode& node = tree.nodes[i];
        if (node.isLeaf)
        {
            continue;
        }
        if (node.left < 0 || node.left >= count || node.right < 0 || node.right >= count
            || node.left == i || node.right == i)
        {
            throw std::runtime_error("genotype-quality model: tree child index out of range");
        }
    }
    return tree;
}

QHead parseQHead(const nlohmann::json& jq)
{
    QHead head;
    head.baseline = jq.at("baseline").get<double>();
    for (const nlohmann::json& jt : jq.at("trees"))
    {
        head.trees.push_back(parseTree(jt));
    }
    return head;
}

Isotonic parseIsotonic(const nlohmann::json& ji)
{
    Isotonic iso;
    iso.x = ji.at("x").get<std::vector<double>>();
    iso.y = ji.at("y").get<std::vector<double>>();
    if (iso.x.size() != iso.y.size())
    {
        throw std::runtime_error("genotype-quality model: isotonic x/y length mismatch");
    }
    return iso;
}

// Resolves the position of a required class label within the direction head's
// explicit `classes` array. Throws if it is absent (which also rejects duplicate
// labels, since a duplicate displaces one of the three required names).
int classPosition(const std::vector<std::string>& classes, const std::string& name)
{
    const auto it = std::find(classes.begin(), classes.end(), name);
    if (it == classes.end())
    {
        throw std::runtime_error("genotype-quality model: direction classes missing '" + name + "'");
    }
    return static_cast<int>(it - classes.begin());
}

DirectionHead parseDirectionHead(const nlohmann::json& jd)
{
    DirectionHead head;
    const std::vector<double> baseline = jd.at("baseline").get<std::vector<double>>();
    if (baseline.size() != 3)
    {
        throw std::runtime_error("genotype-quality model: direction baseline must have 3 entries");
    }
    for (int c = 0; c < 3; ++c)
    {
        head.baseline[c] = baseline[c];
    }

    // Explicit class labels for the parallel baseline/trees/calibrators arrays, so the
    // consumer maps probabilities to OK / TOO_LONG / TOO_SHORT by name, not by position.
    const std::vector<std::string> classes = jd.at("classes").get<std::vector<std::string>>();
    if (classes.size() != 3)
    {
        throw std::runtime_error("genotype-quality model: direction classes must have 3 entries");
    }
    head.posOk = classPosition(classes, "OK");
    head.posTooLong = classPosition(classes, "TOO_LONG");
    head.posTooShort = classPosition(classes, "TOO_SHORT");
    for (const nlohmann::json& jtriple : jd.at("trees"))
    {
        if (jtriple.size() != 3)
        {
            throw std::runtime_error("genotype-quality model: direction tree group must have 3 trees");
        }
        std::array<Tree, 3> triple = {{parseTree(jtriple[0]), parseTree(jtriple[1]), parseTree(jtriple[2])}};
        head.trees.push_back(std::move(triple));
    }
    const nlohmann::json& jcal = jd.at("calibrators");
    if (jcal.size() != 3)
    {
        throw std::runtime_error("genotype-quality model: direction must have 3 calibrators");
    }
    for (int c = 0; c < 3; ++c)
    {
        head.calibrators[c] = parseIsotonic(jcal[c]);
    }
    return head;
}

GenotypingRegimeModel parseGenotypingRegime(const nlohmann::json& jr)
{
    GenotypingRegimeModel genotypingRegime;
    genotypingRegime.qMedian = parseQHead(jr.at("q_median"));
    genotypingRegime.direction = parseDirectionHead(jr.at("direction"));
    return genotypingRegime;
}

void validateTreeFeatureIndices(const Tree& tree, std::size_t featureCount)
{
    for (const TreeNode& node : tree.nodes)
    {
        if (!node.isLeaf && (node.feature < 0 || static_cast<std::size_t>(node.feature) >= featureCount))
        {
            throw std::runtime_error(
                "genotype-quality model: tree feature index " + std::to_string(node.feature)
                + " out of range for a " + std::to_string(featureCount) + "-feature genotypingRegime");
        }
    }
}

// Rejects a model whose tree feature indices fall outside its declared feature_names
// list for that genotypingRegime, so a layout-skewed model fails at load instead of indexing
// past the assembled feature vector at eval time.
void validateGenotypingRegimeFeatureIndices(const GenotypingRegimeModel& genotypingRegime, std::size_t featureCount)
{
    for (const Tree& tree : genotypingRegime.qMedian.trees)
    {
        validateTreeFeatureIndices(tree, featureCount);
    }
    for (const std::array<Tree, 3>& triple : genotypingRegime.direction.trees)
    {
        for (const Tree& tree : triple)
        {
            validateTreeFeatureIndices(tree, featureCount);
        }
    }
}

std::string readMaybeGzip(const std::string& path)
{
    // gzread reads both gzip-compressed and plain files transparently.
    gzFile file = gzopen(path.c_str(), "rb");
    if (file == nullptr)
    {
        throw std::runtime_error("genotype-quality model: cannot open " + path);
    }
    std::string contents;
    char buffer[1 << 16];
    int n;
    while ((n = gzread(file, buffer, sizeof(buffer))) > 0)
    {
        contents.append(buffer, static_cast<std::size_t>(n));
    }
    const bool readError = (n < 0);
    gzclose(file);
    if (readError)
    {
        throw std::runtime_error("genotype-quality model: error reading " + path);
    }
    return contents;
}

} // namespace

GenotypeQualityModel GenotypeQualityModel::fromJson(const nlohmann::json& j)
{
    GenotypeQualityModel model;
    model.formatVersion = j.at("format_version").get<int>();
    if (model.formatVersion != 2)
    {
        throw std::runtime_error(
            "genotype-quality model: unsupported format_version " + std::to_string(model.formatVersion));
    }
    const nlohmann::json& jfeat = j.at("feature_names");
    model.featureNamesQuick = jfeat.at("quick").get<std::vector<std::string>>();
    model.featureNamesFull = jfeat.at("full").get<std::vector<std::string>>();

    const nlohmann::json& jgenotyping_regimes = j.at("genotyping_regimes");
    model.quick = parseGenotypingRegime(jgenotyping_regimes.at("quick"));
    model.fullSpanning = parseGenotypingRegime(jgenotyping_regimes.at("full_spanning"));
    model.fullNonspanning = parseGenotypingRegime(jgenotyping_regimes.at("full_nonspanning"));

    // The quick genotypingRegime indexes the quick feature_names; both full genotyping_regimes index the full
    // list. Reject any tree whose feature index is out of range for its genotypingRegime.
    validateGenotypingRegimeFeatureIndices(model.quick, model.featureNamesQuick.size());
    validateGenotypingRegimeFeatureIndices(model.fullSpanning, model.featureNamesFull.size());
    validateGenotypingRegimeFeatureIndices(model.fullNonspanning, model.featureNamesFull.size());
    return model;
}

GenotypeQualityModel loadGenotypeQualityModel(const std::string& path)
{
    const std::string contents = readMaybeGzip(path);
    nlohmann::json j;
    try
    {
        j = nlohmann::json::parse(contents);
    }
    catch (const nlohmann::json::parse_error& e)
    {
        throw std::runtime_error("genotype-quality model: JSON parse error in " + path + ": " + e.what());
    }
    return GenotypeQualityModel::fromJson(j);
}

namespace
{

// Inflates an in-memory gzip buffer (used for the embedded model blob).
std::string gunzip(const unsigned char* data, std::size_t length)
{
    z_stream stream{};
    // 15 + 16 selects the gzip wrapper (vs raw/zlib) for the window size.
    if (inflateInit2(&stream, 15 + 16) != Z_OK)
    {
        throw std::runtime_error("genotype-quality model: failed to initialize gzip inflate");
    }
    stream.next_in = const_cast<Bytef*>(data);
    stream.avail_in = static_cast<uInt>(length);

    std::string out;
    char buffer[1 << 16];
    int ret;
    do
    {
        stream.next_out = reinterpret_cast<Bytef*>(buffer);
        stream.avail_out = sizeof(buffer);
        ret = inflate(&stream, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END)
        {
            inflateEnd(&stream);
            throw std::runtime_error("genotype-quality model: corrupt embedded gzip blob");
        }
        out.append(buffer, sizeof(buffer) - stream.avail_out);
    } while (ret != Z_STREAM_END);
    inflateEnd(&stream);
    return out;
}

} // namespace

// Defined in the generated GenotypeQualityModelData.cpp (byte array of the
// checked-in ehunter/data/genotype_quality_model.json.gz; length 0 when the
// placeholder is empty).
extern const unsigned char kEmbeddedModelGz[];
extern const std::size_t kEmbeddedModelGzLen;

std::shared_ptr<const GenotypeQualityModel> loadEmbeddedGenotypeQualityModel()
{
    if (kEmbeddedModelGzLen == 0)
    {
        return nullptr;
    }
    const std::string contents = gunzip(kEmbeddedModelGz, kEmbeddedModelGzLen);
    nlohmann::json j;
    try
    {
        j = nlohmann::json::parse(contents);
    }
    catch (const nlohmann::json::parse_error& e)
    {
        throw std::runtime_error(
            std::string("genotype-quality model: JSON parse error in embedded model: ") + e.what());
    }
    return std::make_shared<const GenotypeQualityModel>(GenotypeQualityModel::fromJson(j));
}

} // namespace gq
} // namespace ehunter
