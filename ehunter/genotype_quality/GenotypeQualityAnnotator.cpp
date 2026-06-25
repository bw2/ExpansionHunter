//
// ExpansionHunter
//
// Implementation of the per-allele genotype-quality prediction glue.
//

#include "genotype_quality/GenotypeQualityAnnotator.hh"

#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ehunter
{
namespace gq
{

namespace
{

const GenotypingRegimeModel& genotypingRegimeModelFor(const GenotypeQualityModel& model, GenotypingRegime genotypingRegime)
{
    switch (genotypingRegime)
    {
    case GenotypingRegime::Quick:
        return model.quick;
    case GenotypingRegime::FullSpanning:
        return model.fullSpanning;
    case GenotypingRegime::FullNonspanning:
        return model.fullNonspanning;
    }
    return model.quick; // unreachable; silences -Werror
}

} // namespace

AllelePrediction predictAllele(
    const GenotypeQualityModel& model, bool quickGenotype, const LocusFeatureContext& ctx,
    int alleleRank, int eh, int ciStart, int ciEnd, const AlleleMetrics& aqm)
{
    const int spanningAtCalled = ctx.spanningReads.countOf(eh);
    const GenotypingRegime genotypingRegime = genotypingRegimeOf(quickGenotype, spanningAtCalled);

    // Assemble this allele's features in the C++ canonical order, then materialize the
    // model's feature vector in the order the model itself declares (looked up by name),
    // so the model and binary never rely on a shared positional layout. A model that asks
    // for a feature this build does not produce is rejected rather than mis-evaluated.
    const std::vector<double> canonValues = assembleFeatures(ctx, alleleRank, eh, ciStart, ciEnd, aqm, genotypingRegime);
    const std::vector<std::string>& canonNames = featureNamesForGenotypingRegime(genotypingRegime);
    std::unordered_map<std::string, double> valueByName;
    valueByName.reserve(canonNames.size());
    for (std::size_t i = 0; i < canonNames.size(); ++i)
    {
        valueByName.emplace(canonNames[i], canonValues[i]);
    }

    const std::vector<std::string>& modelNames
        = (genotypingRegime == GenotypingRegime::Quick) ? model.featureNamesQuick : model.featureNamesFull;
    std::vector<double> x;
    x.reserve(modelNames.size());
    for (const std::string& name : modelNames)
    {
        const auto it = valueByName.find(name);
        if (it == valueByName.end())
        {
            throw std::runtime_error(
                "genotype-quality model requires feature '" + name
                + "' that ExpansionHunter does not produce");
        }
        x.push_back(it->second);
    }

    return model.predict(genotypingRegimeModelFor(model, genotypingRegime), x);
}

} // namespace gq
} // namespace ehunter
