//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// Adapted from REViewer (Copyright 2020 Illumina, Inc., GPL-3.0 license)
//

#include "reviewer/Phasing.hh"

#include <algorithm>
#include <cassert>
#include <map>

#include "reviewer/Projection.hh"

namespace ehunter
{
namespace reviewer
{

using std::pair;
using std::string;
using std::vector;

static int scorePath(int pathIndex, const PairPathAlignById& pairPathAlignById)
{
    int pathScore = 0;
    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& pairAlign = idAndPairPathAlign.second;
        for (const auto& readAlign : pairAlign.readAligns)
        {
            if (readAlign.pathIndex == pathIndex)
            {
                pathScore += score(*readAlign.align);
                break;
            }
        }

        for (const auto& mateAlign : pairAlign.mateAligns)
        {
            if (mateAlign.pathIndex == pathIndex)
            {
                pathScore += score(*mateAlign.align);
                break;
            }
        }
    }

    return pathScore;
}

ScoredDiplotypes scoreDiplotypes(const FragById& fragById, const vector<Diplotype>& diplotypes)
{
    vector<ScoredDiplotype> scoredDiplotypes;

    for (const auto& diplotype : diplotypes)
    {
        assert(diplotype.size() == 1 || diplotype.size() == 2);

        auto pairPathAlignById = project(diplotype, fragById);

        int genotypeScore = scorePath(0, pairPathAlignById);
        if (diplotype.size() == 2)
        {
            genotypeScore += scorePath(1, pairPathAlignById);
        }

        scoredDiplotypes.emplace_back(diplotype, genotypeScore);
    }

    // stable_sort so tied diplotypes keep the (sorted) candidate order, making the top pick the same on every
    // platform; std::sort's order for ties depends on the standard library implementation.
    std::stable_sort(
        scoredDiplotypes.begin(), scoredDiplotypes.end(),
        [](const ScoredDiplotype& gt1, const ScoredDiplotype& gt2) { return gt1.second > gt2.second; });

    assert(!scoredDiplotypes.empty());
    return scoredDiplotypes;
}

/// Best pair alignment score of each fragment across the haplotypes of the given diplotype. project() keeps only
/// each fragment's best-scoring haplotype alignments, so the first read and mate alignments carry that score.
static std::map<string, int> getBestPairScoreByFragId(const Diplotype& diplotype, const FragById& fragById)
{
    std::map<string, int> bestPairScoreByFragId;
    for (const auto& fragIdAndPairPathAlign : project(diplotype, fragById))
    {
        const PairPathAlign& pairPathAlign = fragIdAndPairPathAlign.second;
        bestPairScoreByFragId.emplace(
            fragIdAndPairPathAlign.first,
            score(*pairPathAlign.readAligns.front().align) + score(*pairPathAlign.mateAligns.front().align));
    }

    return bestPairScoreByFragId;
}

DiplotypeChoiceSupport compareTopTwoDiplotypes(const FragById& fragById, const ScoredDiplotypes& scoredDiplotypes)
{
    DiplotypeChoiceSupport support;
    support.numberOfCandidateDiplotypes = static_cast<int>(scoredDiplotypes.size());
    if (scoredDiplotypes.size() < 2)
    {
        return support;
    }

    support.topDiplotypeIsTiedWithNextBest = scoredDiplotypes[0].second == scoredDiplotypes[1].second;

    const auto topScoreByFragId = getBestPairScoreByFragId(scoredDiplotypes[0].first, fragById);
    const auto nextBestScoreByFragId = getBestPairScoreByFragId(scoredDiplotypes[1].first, fragById);
    for (const auto& fragIdAndFrag : fragById)
    {
        // A fragment that projects onto no haplotype of a diplotype counts as fitting that diplotype worst
        const auto topIt = topScoreByFragId.find(fragIdAndFrag.first);
        const auto nextBestIt = nextBestScoreByFragId.find(fragIdAndFrag.first);
        const bool fitsTop = topIt != topScoreByFragId.end();
        const bool fitsNextBest = nextBestIt != nextBestScoreByFragId.end();

        if (fitsTop && (!fitsNextBest || topIt->second > nextBestIt->second))
        {
            ++support.fragmentsFavoringTopDiplotype;
        }
        else if (fitsNextBest && (!fitsTop || nextBestIt->second > topIt->second))
        {
            ++support.fragmentsFavoringNextBestDiplotype;
        }
    }

    return support;
}

std::optional<RepeatAllelePhasing> summarizeRepeatAllelePhasing(
    const LocusSpecification& locusSpec, const LocusFindings& findings, const Diplotype& chosenDiplotype,
    const DiplotypeChoiceSupport& support)
{
    // Repeat sizes on each haplotype, in catalog variant order
    vector<vector<int>> repeatSizesOnEachHaplotype(chosenDiplotype.size());
    vector<string> repeatVariantIds;
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type != VariantType::kRepeat)
        {
            continue;
        }

        const auto findingsIt = findings.findingsForEachVariant.find(variantSpec.id());
        const RepeatFindings* repeatFindings = findingsIt == findings.findingsForEachVariant.end()
            ? nullptr
            : dynamic_cast<const RepeatFindings*>(findingsIt->second.get());
        if (!repeatFindings || !repeatFindings->optionalGenotype())
        {
            return std::nullopt;
        }
        const RepeatGenotype& genotype = *repeatFindings->optionalGenotype();
        repeatVariantIds.push_back(variantSpec.id());

        // The haplotype paths hold the genotype's allele sizes, possibly capped (see capLengths in
        // GenotypePaths.cpp), so the shorter path count marks the haplotype carrying the short allele.
        const graphtools::NodeId repeatNode = variantSpec.nodes().front();
        vector<int> repeatNodeCountOnEachHaplotype;
        for (const auto& haplotypePath : chosenDiplotype)
        {
            const auto& nodeIds = haplotypePath.nodeIds();
            repeatNodeCountOnEachHaplotype.push_back(
                static_cast<int>(std::count(nodeIds.begin(), nodeIds.end(), repeatNode)));
        }

        const bool isHeterozygous = genotype.shortAlleleSizeInUnits() != genotype.longAlleleSizeInUnits();
        if (isHeterozygous
            && (chosenDiplotype.size() != 2 || repeatNodeCountOnEachHaplotype[0] == repeatNodeCountOnEachHaplotype[1]))
        {
            return std::nullopt;
        }

        for (size_t haplotypeIndex = 0; haplotypeIndex != chosenDiplotype.size(); ++haplotypeIndex)
        {
            const bool carriesShortAllele = !isHeterozygous
                || repeatNodeCountOnEachHaplotype[haplotypeIndex]
                    < repeatNodeCountOnEachHaplotype[1 - haplotypeIndex];
            repeatSizesOnEachHaplotype[haplotypeIndex].push_back(
                carriesShortAllele ? genotype.shortAlleleSizeInUnits() : genotype.longAlleleSizeInUnits());
        }
    }

    std::sort(repeatSizesOnEachHaplotype.begin(), repeatSizesOnEachHaplotype.end());

    RepeatAllelePhasing phasing;
    for (const auto& repeatSizes : repeatSizesOnEachHaplotype)
    {
        std::map<string, int> repeatSizeByVariantId;
        for (size_t variantIndex = 0; variantIndex != repeatVariantIds.size(); ++variantIndex)
        {
            repeatSizeByVariantId.emplace(repeatVariantIds[variantIndex], repeatSizes[variantIndex]);
        }
        phasing.repeatSizeByVariantIdOnEachHaplotype.push_back(std::move(repeatSizeByVariantId));
    }
    phasing.numberOfPossiblePairings = support.numberOfCandidateDiplotypes;
    phasing.fragmentsSupportingChosenPairingOverNextBest = support.fragmentsFavoringTopDiplotype;
    phasing.fragmentsSupportingNextBestPairingOverChosen = support.fragmentsFavoringNextBestDiplotype;
    phasing.chosenPairingIsTiedWithNextBest = support.topDiplotypeIsTiedWithNextBest;

    return phasing;
}

}  // namespace reviewer
}  // namespace ehunter
