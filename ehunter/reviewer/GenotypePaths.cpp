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

#include "reviewer/GenotypePaths.hh"

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

namespace ehunter
{
namespace reviewer
{

using graphtools::Graph;
using graphtools::NodeId;
using graphtools::Path;
using std::map;
using std::pair;
using std::runtime_error;
using std::string;
using std::vector;

using NodeRange = pair<NodeId, NodeId>;
using Nodes = std::vector<graphtools::NodeId>;
using NodeVector = std::vector<graphtools::NodeId>;
using NodeVectors = std::vector<NodeVector>;

/// Extract repeat lengths from LocusFindings (instead of VCF file)
static vector<int> extractRepeatLengths(const string& variantId, const LocusFindings& findings)
{
    auto it = findings.findingsForEachVariant.find(variantId);
    if (it == findings.findingsForEachVariant.end())
    {
        throw runtime_error("No findings for variant " + variantId);
    }

    const VariantFindings* findingsPtr = it->second.get();
    const RepeatFindings* repeatFindings = dynamic_cast<const RepeatFindings*>(findingsPtr);
    if (!repeatFindings)
    {
        throw runtime_error("Variant " + variantId + " is not a repeat variant");
    }

    const auto& optionalGenotype = repeatFindings->optionalGenotype();
    if (!optionalGenotype)
    {
        throw runtime_error("Cannot create a plot because the genotype of " + variantId + " is missing");
    }

    vector<int> sizes;
    sizes.push_back(optionalGenotype->shortAlleleSizeInUnits());
    if (optionalGenotype->numAlleles() == 2)
    {
        sizes.push_back(optionalGenotype->longAlleleSizeInUnits());
    }
    return sizes;
}

static vector<int> capLengths(int upperBound, const vector<int>& lengths)
{
    vector<int> cappedLength;
    cappedLength.reserve(lengths.size());
    for (int length : lengths)
    {
        cappedLength.push_back(length <= upperBound ? length : upperBound);
    }

    return cappedLength;
}

/// Determine sequences of nodes corresponding to each allele of the given variant
/// @param meanFragLen Mean fragment length
/// @param locusSpec Description of the target locus
/// @param findings Locus findings with genotypes
/// @return Sequences of nodes for each allele indexed by the range of nodes corresponding to the entire variant
///
/// Assumption: Locus contains only STRs
/// Detail: STR lengths are capped by fragment length (see implementation)
/// Example:
///  An STR corresponding to RE (CAG)* with genotype 3/4 corresponds to the
///  output {{1, 1}: {{1, 1, 1}, {1, 1, 1, 1}}
static map<NodeRange, NodeVectors>
getGenotypeNodesByNodeRange(int meanFragLen, const LocusSpecification& locusSpec, const LocusFindings& findings)
{
    map<NodeRange, NodeVectors> genotypeNodesByNodeRange;
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            throw std::logic_error("REViewer does not accept locus definitions containing small variants (e.g. '(A|T)').");
        }
        assert(variantSpec.classification().type == VariantType::kRepeat);
        assert(variantSpec.nodes().size() == 1);
        NodeId repeatNode = variantSpec.nodes().front();
        NodeVectors genotypeNodes;

        auto repeatLens = extractRepeatLengths(variantSpec.id(), findings);
        repeatLens = capLengths(meanFragLen, repeatLens);

        genotypeNodes.reserve(repeatLens.size());
        for (int repeatLen : repeatLens)
        {
            genotypeNodes.emplace_back(repeatLen, repeatNode);
        }

        NodeId nodeRangeFrom = std::numeric_limits<NodeId>::max();
        NodeId nodeRangeTo = std::numeric_limits<NodeId>::lowest();

        for (NodeId node : variantSpec.nodes())
        {
            if (node != static_cast<NodeId>(-1))
            {
                nodeRangeFrom = std::min(nodeRangeFrom, node);
                nodeRangeTo = std::max(nodeRangeTo, node);
            }
        }

        assert(genotypeNodes.size() <= 2);
        genotypeNodesByNodeRange.emplace(std::make_pair(nodeRangeFrom, nodeRangeTo), genotypeNodes);
    }

    return genotypeNodesByNodeRange;
}

static boost::optional<pair<NodeVectors, NodeId>>
getVariantGenotypeNodes(const map<NodeRange, NodeVectors>& nodeRangeToPaths, NodeId node)
{
    for (const auto& nodeRangeAndPaths : nodeRangeToPaths)
    {
        const auto& nodeRange = nodeRangeAndPaths.first;
        const auto& paths = nodeRangeAndPaths.second;
        if (nodeRange.first <= node && node <= nodeRange.second)
        {
            assert(paths.size() <= 2);
            return std::make_pair(paths, nodeRange.second);
        }
    }

    return boost::none;
}

static vector<NodeVectors> extendDiplotype(const vector<NodeVectors>& genotypes, const NodeVectors& genotypeExtension)
{
    vector<NodeVectors> extendedGenotype;
    for (const auto& genotype : genotypes)
    {
        assert(genotype.size() == genotypeExtension.size());
        if (genotype.size() == 1)
        {
            const Nodes& haplotypeExtension = genotypeExtension.front();
            Nodes extendedHaplotype = genotype.front();
            extendedHaplotype.insert(extendedHaplotype.end(), haplotypeExtension.begin(), haplotypeExtension.end());
            extendedGenotype.push_back({ extendedHaplotype });
        }
        else
        {
            assert(genotype.size() == 2);
            Nodes hap1Ext1 = genotype.front();
            hap1Ext1.insert(hap1Ext1.end(), genotypeExtension.front().begin(), genotypeExtension.front().end());

            Nodes hap2Ext2 = genotype.back();
            hap2Ext2.insert(hap2Ext2.end(), genotypeExtension.back().begin(), genotypeExtension.back().end());

            extendedGenotype.push_back({ hap1Ext1, hap2Ext2 });

            Nodes hap1Ext2 = genotype.front();
            hap1Ext2.insert(hap1Ext2.end(), genotypeExtension.back().begin(), genotypeExtension.back().end());

            Nodes hap2Ext1 = genotype.back();
            hap2Ext1.insert(hap2Ext1.end(), genotypeExtension.front().begin(), genotypeExtension.front().end());

            extendedGenotype.push_back({ hap1Ext2, hap2Ext1 });
        }
    }

    return extendedGenotype;
}

std::vector<Diplotype> getCandidateDiplotypes(
    int meanFragLen,
    const LocusSpecification& locusSpec,
    const LocusFindings& findings)
{
    auto genotypeNodesByNodeRange = getGenotypeNodesByNodeRange(meanFragLen, locusSpec, findings);

    // Assume that all variants have the same number of alleles
    const auto numAlleles = genotypeNodesByNodeRange.empty() ? 2 : genotypeNodesByNodeRange.begin()->second.size();

    vector<NodeVectors> nodesByDiplotype = { NodeVectors(numAlleles, { 0 }) };

    NodeId node = 1;
    while (node != locusSpec.regionGraph().numNodes())
    {
        auto variantPathNodesAndLastNode = getVariantGenotypeNodes(genotypeNodesByNodeRange, node);
        if (variantPathNodesAndLastNode)
        {
            nodesByDiplotype = extendDiplotype(nodesByDiplotype, variantPathNodesAndLastNode->first);
            node = variantPathNodesAndLastNode->second;
        }
        else
        {
            for (auto& genotypeNodes : nodesByDiplotype)
            {
                for (auto& haplotypeNodes : genotypeNodes)
                {
                    haplotypeNodes.push_back(node);
                }
            }
        }

        ++node;
    }

    vector<Diplotype> diplotypes;
    const NodeId rightFlankNode = locusSpec.regionGraph().numNodes() - 1;
    const int rightFlankLength = static_cast<int>(locusSpec.regionGraph().nodeSeq(rightFlankNode).length());
    for (const auto& diplotypeNodes : nodesByDiplotype)
    {
        Diplotype diplotype;
        for (const auto& haplotypeNodes : diplotypeNodes)
        {
            diplotype.emplace_back(&locusSpec.regionGraph(), 0, haplotypeNodes, rightFlankLength);
        }

        // The code so far considers diplotypes that differ by the order of constituent haplotypes to be distinct.
        // To overcome this issue, we enforce a consistent haplotype order.
        if (diplotype.front() < diplotype.back())
        {
            std::iter_swap(diplotype.begin(), diplotype.end() - 1);
        }

        diplotypes.push_back(diplotype);
    }

    std::sort(diplotypes.begin(), diplotypes.end());
    diplotypes.erase(std::unique(diplotypes.begin(), diplotypes.end()), diplotypes.end());

    return diplotypes;
}

static string summarizePath(const graphtools::Path& path)
{
    const auto& graph = *path.graphRawPtr();
    string summary;
    std::set<graphtools::NodeId> observedNodes;
    for (const auto nodeId : path.nodeIds())
    {
        if (observedNodes.find(nodeId) != observedNodes.end())
        {
            continue;
        }

        const bool isLoopNode = graph.hasEdge(nodeId, nodeId);

        if (nodeId == 0)
        {
            summary += "(LF)";
        }
        else if (nodeId + 1 == graph.numNodes())
        {
            summary += "(RF)";
        }
        else
        {
            assert(graph.numNodes() != 0);
            const string& nodeSeq = graph.nodeSeq(nodeId);
            summary += "(" + nodeSeq + ")";

            if (isLoopNode)
            {
                int numMotifs = static_cast<int>(std::count(path.nodeIds().begin(), path.nodeIds().end(), nodeId));
                summary += "{" + std::to_string(numMotifs) + "}";
            }
        }

        observedNodes.emplace(nodeId);
    }

    return summary;
}

std::ostream& operator<<(std::ostream& out, const Diplotype& diplotype)
{
    out << summarizePath(diplotype.front());
    if (diplotype.size() == 2)
    {
        out << "/" << summarizePath(diplotype.back());
    }
    return out;
}

}  // namespace reviewer
}  // namespace ehunter
