//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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
//

#include "io/LocusSpecDecoding.hh"

#include <algorithm>
#include <memory>
#include <stdexcept>

#include <boost/optional.hpp>

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using boost::optional;
using graphtools::Graph;
using graphtools::NodeId;
using std::string;
using std::to_string;
using std::vector;

namespace ehunter
{

static vector<GenomicRegion> addFlankingRegions(int kExtensionLength, const vector<GenomicRegion>& referenceRegions)
{
    const GenomicRegion& firstRegion = referenceRegions.front();
    const int64_t leftFlankStart = firstRegion.start() - kExtensionLength;
    GenomicRegion leftFlank(firstRegion.contigIndex(), leftFlankStart, firstRegion.start());

    const GenomicRegion& lastRegion = referenceRegions.back();
    const int64_t rightFlankEnd = lastRegion.end() + kExtensionLength;
    GenomicRegion rightFlank(lastRegion.contigIndex(), lastRegion.end(), rightFlankEnd);

    auto regions = referenceRegions;
    regions.insert(regions.begin(), std::move(leftFlank));
    regions.push_back(std::move(rightFlank));

    return regions;
}

static string extendLocusStructure(
    Reference& reference, const vector<GenomicRegion>& referenceRegions, const string& flanklessLocusStructure)
{

    const auto& leftFlankRegion = referenceRegions.front();
    string leftFlank = reference.getSequence(leftFlankRegion);

    const auto& rightFlankRegion = referenceRegions.back();
    string rightFlank = reference.getSequence(rightFlankRegion);

    const size_t maxNsAllowedInFlanks = 5;
    const size_t numNsInLeftFlank = std::count(leftFlank.begin(), leftFlank.end(), 'N');
    const size_t numNsInRightFlank = std::count(rightFlank.begin(), rightFlank.end(), 'N');

    if (numNsInLeftFlank + numNsInRightFlank > maxNsAllowedInFlanks)
    {
        const string message = "Flanks can contain at most " + to_string(maxNsAllowedInFlanks)
            + " characters N but found " + to_string(numNsInLeftFlank + numNsInRightFlank) + " Ns";
        throw std::runtime_error(message);
    }

    return leftFlank + flanklessLocusStructure + rightFlank;
}

static vector<GenomicRegion> addReferenceRegionsForInterruptions(const GraphBlueprint& blueprint, const vector<GenomicRegion>& referenceRegions)
{
    int regionIndex = 0;
    vector<GenomicRegion> completedReferenceRegions;

    for (const auto& feature : blueprint)
    {
        if (feature.type == GraphBlueprintFeatureType::kInterruption)
        {
            assert(regionIndex != 0 && regionIndex < static_cast<int>(referenceRegions.size()));
            const auto& leftRegion = referenceRegions[regionIndex];
            const auto& rightRegion = referenceRegions[regionIndex + 1];
            completedReferenceRegions.emplace_back(leftRegion.contigIndex(), leftRegion.end(), rightRegion.start());
        }
        else
        {
            completedReferenceRegions.push_back(referenceRegions[regionIndex]);
            ++regionIndex;
        }
    }

    assert(blueprint.size() == completedReferenceRegions.size());
    return completedReferenceRegions;
}

static GenomicRegion mergeRegions(const vector<GenomicRegion>& regions)
{
    const int kMaxMergeDistance = 500;
    vector<GenomicRegion> mergedReferenceRegions = merge(regions, kMaxMergeDistance);
    if (mergedReferenceRegions.size() != 1)
    {
        std::stringstream out;
        for (const GenomicRegion& region : regions)
        {
            out << region << " ";
        }
        throw std::runtime_error(
            "Expected reference regions to be closer than " + to_string(kMaxMergeDistance)
            + " from one another: " + out.str());
    }

    return mergedReferenceRegions.front();
}

static NodeToRegionAssociation associateNodesWithReferenceRegions(
    const GraphBlueprint& blueprint, const Graph& graph, const vector<GenomicRegion>& referenceRegions)
{
    assert(blueprint.size() == referenceRegions.size());

    NodeToRegionAssociation referenceRegionsOfGraphNodes;

    for (int featureIndex = 0; featureIndex != static_cast<int>(blueprint.size()); ++featureIndex)
    {
        const auto& feature = blueprint[featureIndex];
        const auto& referenceRegion = referenceRegions[featureIndex];

        for (const auto& nodeId : feature.nodeIds)
        {
            const int nodeLength = graph.nodeSeq(nodeId).length();
            GenomicRegion referenceRegionForNode(
                referenceRegion.contigIndex(), referenceRegion.start(), referenceRegion.start() + nodeLength);
            referenceRegionsOfGraphNodes.emplace(std::make_pair(nodeId, std::move(referenceRegionForNode)));
        }
    }

    return referenceRegionsOfGraphNodes;
}

static VariantType determineVariantType(GraphBlueprintFeatureType featureType)
{
    switch (featureType)
    {
    case GraphBlueprintFeatureType::kInsertionOrDeletion:
    case GraphBlueprintFeatureType::kSwap:
        return VariantType::kSmallVariant;
    case GraphBlueprintFeatureType::kSkippableRepeat:
    case GraphBlueprintFeatureType::kUnskippableRepeat:
        return VariantType::kRepeat;
    default:
        std::ostringstream encoding;
        encoding << featureType;
        throw std::logic_error("Feature of type " + encoding.str() + " does not define a variant");
    }
}

static VariantSubtype determineVariantSubtype(
    GraphBlueprintFeatureType featureType, VariantTypeFromUser variantTypeFromUser, const GenomicRegion referenceRegion)
{
    if (featureType == GraphBlueprintFeatureType::kInsertionOrDeletion)
    {
        if (referenceRegion.length() == 0)
        {
            return VariantSubtype::kInsertion;
        }
        else
        {
            return VariantSubtype::kDeletion;
        }
    }
    else if (featureType == GraphBlueprintFeatureType::kSwap)
    {
        return VariantSubtype::kSwap;
    }
    else if (variantTypeFromUser == VariantTypeFromUser::kCommonRepeat)
    {
        return VariantSubtype::kCommonRepeat;
    }
    else if (variantTypeFromUser == VariantTypeFromUser::kRareRepeat)
    {
        return VariantSubtype::kRareRepeat;
    }
    else
    {
        std::ostringstream message;
        message << featureType;
        throw std::logic_error("Feature " + message.str() + " does not correspond to variant");
    }
}

static optional<NodeId> determineReferenceNode(
    const GraphBlueprintFeature& feature, Reference& reference, const GenomicRegion& referenceRegion)
{

    if (feature.type == GraphBlueprintFeatureType::kSkippableRepeat
        || feature.type == GraphBlueprintFeatureType::kUnskippableRepeat)
    {
        return feature.nodeIds.front();
    }

    const auto& contigName = reference.contigInfo().getContigName(referenceRegion.contigIndex());
    const string refSequence = reference.getSequence(contigName, referenceRegion.start(), referenceRegion.end());

    optional<NodeId> optionalReferenceNode;
    for (int index = 0; index != static_cast<int>(feature.nodeIds.size()); ++index)
    {
        if (refSequence == feature.sequences[index])
        {
            optionalReferenceNode = feature.nodeIds[index];
            break;
        }
    }

    return optionalReferenceNode;
}

LocusSpecification decodeLocusSpecification(
    const LocusDescription& locusDescription, Reference& reference,
    const HeuristicParameters& heuristicParams,
    const bool extendFlanks)
{
    /** Params:
    * 		extendFlanks (bool): if true, extend the locus structure to add flanking sequences from the reference.
    *            if false, skip this time- and memory-consuming step and create a LocusSpecification with just the
    *            repeat regions as specified in the LocusDescription. This is useful for creating a stub LocusSpec for
    *            testing or fast genotyping without using the full EH genotyping algorithm.
    */
    try
    {
        assertValidity(locusDescription);

        const int kExtensionLength = heuristicParams.regionExtensionLength();

        std::string locusStructure;
        std::vector<GenomicRegion> referenceRegions;
        if(extendFlanks) {
            referenceRegions = addFlankingRegions(kExtensionLength, locusDescription.referenceRegions());
            locusStructure = extendLocusStructure(reference, referenceRegions, locusDescription.locusStructure());
        } else {
            locusStructure = locusDescription.locusStructure();
            referenceRegions = locusDescription.referenceRegions();
        }

        GraphBlueprint blueprint = decodeFeaturesFromRegex(locusStructure);
        graphtools::Graph locusGraph = makeRegionGraph(blueprint, locusDescription.locusId());
        auto completeReferenceRegions = addReferenceRegionsForInterruptions(blueprint, referenceRegions);

        GenomicRegion referenceRegionForEntireLocus = mergeRegions(locusDescription.referenceRegions());

        vector<GenomicRegion> targetReadExtractionRegions;
        for (const GenomicRegion& region : locusDescription.targetRegions())
        {
            targetReadExtractionRegions.push_back(region.extend(kExtensionLength));
        }

        if (targetReadExtractionRegions.empty())
        {
            targetReadExtractionRegions.push_back(referenceRegionForEntireLocus.extend(kExtensionLength));
        }

        NodeToRegionAssociation referenceRegionsOfGraphNodes
            = associateNodesWithReferenceRegions(blueprint, locusGraph, completeReferenceRegions);

        GenotyperParameters parameters(heuristicParams.minLocusCoverage());
        if (locusDescription.errorRate())
        {
            parameters.errorRate = *locusDescription.errorRate();
        }
        if (locusDescription.likelihoodRatioThreshold())
        {
            parameters.likelihoodRatioThreshold = *locusDescription.likelihoodRatioThreshold();
        }
        if (locusDescription.minLocusCoverage())
        {
            parameters.minLocusCoverage = *locusDescription.minLocusCoverage();
        }

        LocusSpecification locusSpec(
            locusDescription.locusId(), locusDescription.chromType(), std::move(targetReadExtractionRegions),
            std::move(locusGraph), std::move(referenceRegionsOfGraphNodes), std::move(parameters),
            locusDescription.useRFC1MotifAnalysis(),
            locusDescription.plotConditions(),
            locusDescription.extraFields());
        locusSpec.setOfftargetReadExtractionRegions(locusDescription.offtargetRegions());

        int variantIndex = 0;
        for (const auto& feature : blueprint)
        {
            if (doesFeatureDefineVariant(feature.type))
            {
                const GenomicRegion& referenceRegion = locusDescription.referenceRegions().at(variantIndex);

                VariantTypeFromUser variantDescription = locusDescription.variantTypesFromUser().at(variantIndex);
                const string& variantId = locusDescription.variantIds()[variantIndex];
                VariantType variantType = determineVariantType(feature.type);
                VariantSubtype variantSubtype
                    = determineVariantSubtype(feature.type, variantDescription, referenceRegion);

                optional<NodeId> optionalReferenceNode = determineReferenceNode(feature, reference, referenceRegion);

                VariantClassification classification(variantType, variantSubtype);

                locusSpec.addVariantSpecification(
                    variantId, classification, referenceRegion, feature.nodeIds, optionalReferenceNode);

                ++variantIndex;
            }
        }
        return locusSpec;
    }
    catch (const std::exception& e)
    {
        throw std::runtime_error("Error loading locus " + locusDescription.locusId() + ": " + e.what());
    }
}

void assertValidity(const LocusDescription& locusDescription)
{
    const GraphBlueprint blueprint = decodeFeaturesFromRegex(locusDescription.locusStructure());
    int numVariants = 0;
    for (const GraphBlueprintFeature& feature : blueprint)
    {
        if (doesFeatureDefineVariant(feature.type))
        {
            ++numVariants;
        }
    }

    if (numVariants == 0)
    {
        throw std::runtime_error(
            "Locus " + locusDescription.locusId() + " must encode at least one variant " + locusDescription.locusStructure());
    }

    if (numVariants != static_cast<int>(locusDescription.referenceRegions().size()))
    {
        throw std::runtime_error(
            "Locus " + locusDescription.locusId() + " must specify reference regions for " + to_string(numVariants)
            + " variants");
    }

    if (numVariants != static_cast<int>(locusDescription.variantTypesFromUser().size()))
    {
        throw std::runtime_error(
            "Locus " + locusDescription.locusId() + " must specify variant types for " + to_string(numVariants)
            + " variants");
    }

    if (locusDescription.useRFC1MotifAnalysis())
    {
        if ((numVariants != 1) or (locusDescription.variantTypesFromUser()[0] != VariantTypeFromUser::kCommonRepeat))
        {
            throw std::runtime_error(
                "Locus " + locusDescription.locusId()
                + " has option 'useRFC1MotifAnalysis' enabled, which requires that"
                  " exactly one variant of type 'Repeat' is defined.");
        }
    }
}

}
