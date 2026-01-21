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

#include "reviewer/Projection.hh"

#include <algorithm>
#include <set>
#include <string>

#include "graphalign/Operation.hh"

namespace ehunter
{
namespace reviewer
{

using boost::optional;
using graphtools::Alignment;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::string;
using std::vector;

int score(const GraphAlignment& alignment, int matchScore, int mismatchScore, int gapScore)
{
    int totalScore = 0;
    for (const auto& nodeAlign : alignment.alignments())
    {
        for (const Operation& operation : nodeAlign)
        {
            switch (operation.type())
            {
            case OperationType::kMatch:
                totalScore += matchScore * static_cast<int>(operation.referenceLength());
                break;
            case OperationType::kMismatch:
                totalScore += mismatchScore * static_cast<int>(operation.referenceLength());
                break;
            case OperationType::kInsertionToRef:
                totalScore += gapScore * static_cast<int>(operation.queryLength());
                break;
            case OperationType::kDeletionFromRef:
                totalScore += gapScore * static_cast<int>(operation.referenceLength());
                break;
            default:
                break;
            }
        }
    }

    return totalScore;
}

static bool checkCommonElement(const vector<NodeId>& nodesA, size_t startIndexA, const vector<NodeId>& nodesB, size_t startIndexB)
{
    assert(startIndexA < nodesA.size() && startIndexB < nodesB.size());
    while (startIndexA != nodesA.size() && startIndexB != nodesB.size())
    {
        if (nodesA[startIndexA] < nodesB[startIndexB])
        {
            ++startIndexA;
        }
        else if (nodesB[startIndexB] < nodesA[startIndexA])
        {
            ++startIndexB;
        }
        else
        {
            assert(nodesA[startIndexA] == nodesB[startIndexB]);
            return true;
        }
    }

    return false;
}

static vector<int> getStartIndexes(const Path& projPath, int initialStartIndex, const GraphAlignment& align)
{
    // Check if IRR
    std::set<NodeId> alignNodes(align.path().nodeIds().begin(), align.path().nodeIds().end());
    if (alignNodes.size() != 1)
    {
        return { initialStartIndex };
    }

    NodeId repeatNode = *alignNodes.begin();
    if (!projPath.graphRawPtr()->hasEdge(repeatNode, repeatNode))
    {
        return { initialStartIndex };
    }

    int numMotifsInPath = static_cast<int>(std::count(projPath.begin(), projPath.end(), repeatNode));
    int numMotifsInAlign = static_cast<int>(std::count(align.path().nodeIds().begin(), align.path().nodeIds().end(), repeatNode));

    if (numMotifsInPath <= numMotifsInAlign)
    {
        return { initialStartIndex };
    }

    auto repeatStartIndexIt = std::find(projPath.begin(), projPath.end(), repeatNode);
    assert(repeatStartIndexIt != projPath.end());
    int repeatStartIndex = static_cast<int>(repeatStartIndexIt - projPath.begin());

    vector<int> startIndexes;
    for (int index = 0; index != numMotifsInPath - numMotifsInAlign + 1; ++index)
    {
        startIndexes.push_back(repeatStartIndex + index);
    }

    return startIndexes;
}

struct AlignProj
{
    AlignProj(vector<int> startIndexes, GraphAlignPtr align)
        : startIndexes(std::move(startIndexes))
        , align(std::move(align))
    {
    }

    vector<int> startIndexes;
    GraphAlignPtr align;
};

static optional<AlignProj> projectAlign(const GraphAlign& align, const Path& projPath)
{
    // Find the first common node
    size_t alignNodeIndex = 0;
    size_t projNodeIndex = 0;

    const auto& alignNodes = align.path().nodeIds();
    const auto& projPathNodes = projPath.nodeIds();
    while (alignNodes[alignNodeIndex] != projPathNodes[projNodeIndex])
    {
        if (alignNodes[alignNodeIndex] < projPathNodes[projNodeIndex])
        {
            ++alignNodeIndex;
            if (alignNodeIndex == alignNodes.size())
            {
                return boost::none;
            }
        }
        else if (projPathNodes[projNodeIndex] < alignNodes[alignNodeIndex])
        {
            ++projNodeIndex;
            if (projNodeIndex == projPathNodes.size())
            {
                return boost::none;
            }
        }
        else
        {
            assert(false);
        }
    }

    // Line up the nodes
    graphtools::NodeId commonNode = alignNodes[alignNodeIndex];
    const size_t alignCommonNodeCount = static_cast<size_t>(std::count(alignNodes.begin(), alignNodes.end(), commonNode));
    const size_t projCommonNodeCount = static_cast<size_t>(std::count(projPathNodes.begin(), projPathNodes.end(), commonNode));

    if (alignCommonNodeCount < projCommonNodeCount)
    {
        projNodeIndex += projCommonNodeCount - alignCommonNodeCount;
    }
    else
    {
        alignNodeIndex += alignCommonNodeCount - projCommonNodeCount;
    }

    // Project starting from the common node
    vector<graphtools::NodeId> projNodes;
    vector<Alignment> projAligns;

    size_t projStartNodeIndex = projNodeIndex;
    size_t startingAlignIndex = alignNodeIndex;
    size_t endingAlignIndex = alignNodeIndex;

    while (alignNodeIndex != alignNodes.size())
    {
        if (!checkCommonElement(alignNodes, alignNodeIndex, projPathNodes, projNodeIndex))
        {
            break;
        }

        NodeId alignNode = alignNodes[alignNodeIndex];
        NodeId targetNode = projPathNodes[projNodeIndex];

        if (alignNode == targetNode)
        {
            projAligns.push_back(align.alignments()[alignNodeIndex]);
            projNodes.push_back(targetNode);
            endingAlignIndex = alignNodeIndex;

            ++alignNodeIndex;
            ++projNodeIndex;
            if (alignNodeIndex == alignNodes.size() || projNodeIndex == projPathNodes.size())
            {
                break;
            }
        }
        else if (alignNode < targetNode) // Add an insertion
        {
            auto queryLen = align.alignments()[alignNodeIndex].queryLength();
            auto refStart = projAligns.back().referenceStart();
            auto operations = projAligns.back().operations();
            operations.emplace_back(OperationType::kInsertionToRef, static_cast<uint32_t>(queryLen));
            projAligns.back() = Alignment(refStart, operations);

            ++alignNodeIndex;
            if (alignNodeIndex == alignNodes.size())
            {
                break;
            }
        }
        else if (targetNode < alignNode) // Add a deletion
        {
            auto refLen = align.path().graphRawPtr()->nodeSeq(targetNode).length();
            projAligns.emplace_back(0, std::to_string(refLen) + "D");
            projNodes.push_back(targetNode);
            ++projNodeIndex;

            if (projNodeIndex == projPathNodes.size())
            {
                break;
            }
        }
        else
        {
            assert(false);
        }
    }

    endingAlignIndex = alignNodeIndex - 1;

    // Add left soft clip if needed
    size_t leftSoftclipLen = 0;
    for (size_t i = 0; i != startingAlignIndex; ++i)
    {
        const auto& alignNode = align.alignments()[i];
        leftSoftclipLen += alignNode.queryLength();
    }

    if (leftSoftclipLen)
    {
        size_t refStart = projAligns.front().referenceStart();
        auto operations = projAligns.front().operations();
        operations.emplace_front(OperationType::kSoftclip, static_cast<uint32_t>(leftSoftclipLen));
        projAligns.front() = Alignment(refStart, operations);
    }

    // Add right soft clip if needed
    size_t rightSoftclipLen = 0;
    for (size_t i = endingAlignIndex + 1; i != alignNodes.size(); ++i)
    {
        const auto& nodeAlign = align.alignments()[i];
        rightSoftclipLen += nodeAlign.queryLength();
    }

    if (rightSoftclipLen)
    {
        auto refStart = projAligns.back().referenceStart();
        auto operations = projAligns.back().operations();
        operations.emplace_back(OperationType::kSoftclip, static_cast<uint32_t>(rightSoftclipLen));
        projAligns.back() = Alignment(refStart, operations);
    }

    // Initialize projected alignment
    auto projPathStart = projAligns.front().referenceStart();
    auto projPathEnd = projAligns.back().referenceStart() + projAligns.back().referenceLength();
    Path projAlignPath(projPath.graphRawPtr(), static_cast<int>(projPathStart), projNodes, static_cast<int>(projPathEnd));
    GraphAlignPtr projAlign(new GraphAlignment(projAlignPath, projAligns));

    auto startIndexes = getStartIndexes(projPath, static_cast<int>(projStartNodeIndex), *projAlign);
    return AlignProj(startIndexes, std::move(projAlign));
}

static vector<ReadPathAlign> projectRead(const GraphAlign& align, int pathIndex, const Path& path)
{
    vector<ReadPathAlign> pathAligns;
    optional<AlignProj> alignProj = projectAlign(align, path);
    if (alignProj)
    {
        for (auto startIndex : alignProj->startIndexes)
        {
            pathAligns.emplace_back(path, pathIndex, startIndex, alignProj->align);
        }
    }

    return pathAligns;
}

PairPathAlignById project(const vector<Path>& genotypePaths, const FragById& fragById)
{
    PairPathAlignById pairPathAlignById;
    for (const auto& idAndFrag : fragById)
    {
        const auto& fragId = idAndFrag.first;
        const auto& frag = idAndFrag.second;
        PairPathAlign pathAlign;
        int bestPairScore = std::numeric_limits<int>::lowest();
        for (size_t pathIndex = 0; pathIndex != genotypePaths.size(); ++pathIndex)
        {
            const auto& path = genotypePaths[pathIndex];
            vector<ReadPathAlign> readPathAligns = projectRead(frag.read.align, static_cast<int>(pathIndex), path);
            vector<ReadPathAlign> matePathAligns = projectRead(frag.mate.align, static_cast<int>(pathIndex), path);
            if (readPathAligns.empty() || matePathAligns.empty())
            {
                continue;
            }

            const int pairScore = score(*readPathAligns.front().align) + score(*matePathAligns.front().align);
            if (pairScore > bestPairScore)
            {
                bestPairScore = pairScore;
                pathAlign.readAligns.clear();
                pathAlign.mateAligns.clear();
            }

            if (pairScore == bestPairScore)
            {
                pathAlign.readAligns.insert(pathAlign.readAligns.end(), readPathAligns.begin(), readPathAligns.end());
                pathAlign.mateAligns.insert(pathAlign.mateAligns.end(), matePathAligns.begin(), matePathAligns.end());
            }
        }

        if (!pathAlign.readAligns.empty() && !pathAlign.mateAligns.empty())
        {
            pairPathAlignById.emplace(fragId, pathAlign);
        }
    }

    return pairPathAlignById;
}

}  // namespace reviewer
}  // namespace ehunter
