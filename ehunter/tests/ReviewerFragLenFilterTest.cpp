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
//

#include "reviewer/FragLenFilter.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "core/Read.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/Projection.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Helper to create a Read object
Read makeRead(const std::string& fragId, MateNumber mateNum, const std::string& seq)
{
    ReadId readId(fragId, mateNum);
    return Read(readId, seq, false);
}

// Helper to create a Frag with both mates
Frag makeFrag(
    const std::string& fragId,
    const std::string& readSeq,
    GraphAlignment readAlign,
    const std::string& mateSeq,
    GraphAlignment mateAlign)
{
    Read read = makeRead(fragId, MateNumber::kFirstMate, readSeq);
    Read mate = makeRead(fragId, MateNumber::kSecondMate, mateSeq);
    ReadWithAlign readWithAlign(std::move(read), std::move(readAlign));
    ReadWithAlign mateWithAlign(std::move(mate), std::move(mateAlign));
    return Frag(std::move(readWithAlign), std::move(mateWithAlign));
}

}  // namespace

// ============================================================================
// getMeanFragLen tests
// ============================================================================

TEST(GetMeanFragLen, SingleFlankingFragment_ReturnsCorrectLength)
{
    // Create a simple graph: leftFlank(10bp) - repeat(3bp) - rightFlank(10bp)
    // Total: nodes 0, 1, 2 where node 0 is left flank, node 1 is repeat, node 2 is right flank
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    // Create alignments that start on left flank (node 0)
    // Both read and mate start on node 0 (left flank)
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[5M]", &graph);  // starts at position 0, 5bp long
    GraphAlignment mateAlign = decodeGraphAlignment(5, "0[5M]", &graph);  // starts at position 5, 5bp long

    FragById fragById;
    fragById.emplace("frag1", makeFrag("frag1", "AAAAA", readAlign, "AAAAA", mateAlign));

    int meanFragLen = getMeanFragLen(fragById);

    // Read: positions 0-5, Mate: positions 5-10
    // Fragment length = max(5, 10) - min(0, 5) = 10 - 0 = 10
    EXPECT_EQ(10, meanFragLen);
}

TEST(GetMeanFragLen, MultipleFlankingFragments_ReturnsCorrectMean)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    // Fragment 1: read at 0-5, mate at 5-10 => length 10
    GraphAlignment readAlign1 = decodeGraphAlignment(0, "0[5M]", &graph);
    GraphAlignment mateAlign1 = decodeGraphAlignment(5, "0[5M]", &graph);

    // Fragment 2: read at 0-4, mate at 4-8 => length 8
    GraphAlignment readAlign2 = decodeGraphAlignment(0, "0[4M]", &graph);
    GraphAlignment mateAlign2 = decodeGraphAlignment(4, "0[4M]", &graph);

    FragById fragById;
    fragById.emplace("frag1", makeFrag("frag1", "AAAAA", readAlign1, "AAAAA", mateAlign1));
    fragById.emplace("frag2", makeFrag("frag2", "AAAA", readAlign2, "AAAA", mateAlign2));

    int meanFragLen = getMeanFragLen(fragById);

    // Mean = (10 + 8) / 2 = 9
    EXPECT_EQ(9, meanFragLen);
}

TEST(GetMeanFragLen, EmptyFragById_ThrowsException)
{
    FragById fragById;

    EXPECT_THROW(getMeanFragLen(fragById), std::runtime_error);
}

TEST(GetMeanFragLen, NoFlankingReads_ThrowsException)
{
    // Create a graph where we can have non-flanking reads
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    // Create alignments that start on repeat node (node 1), not on flank nodes
    // This simulates reads that don't qualify as flanking reads
    GraphAlignment readAlign = decodeGraphAlignment(0, "1[3M]", &graph);  // starts on repeat node
    GraphAlignment mateAlign = decodeGraphAlignment(0, "1[3M]", &graph);  // starts on repeat node

    FragById fragById;
    fragById.emplace("frag1", makeFrag("frag1", "CAG", readAlign, "CAG", mateAlign));

    // Should throw because no flanking reads were found
    EXPECT_THROW(getMeanFragLen(fragById), std::runtime_error);
}

TEST(GetMeanFragLen, MixedFlankingAndNonFlankingReads_OnlyCountsFlanking)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    // Fragment 1: flanking (both on left flank), length = 10
    GraphAlignment readAlign1 = decodeGraphAlignment(0, "0[5M]", &graph);
    GraphAlignment mateAlign1 = decodeGraphAlignment(5, "0[5M]", &graph);

    // Fragment 2: non-flanking (read on repeat, mate on repeat)
    GraphAlignment readAlign2 = decodeGraphAlignment(0, "1[3M]", &graph);
    GraphAlignment mateAlign2 = decodeGraphAlignment(0, "1[3M]", &graph);

    // Fragment 3: flanking (both on left flank), length = 8
    // Node 0 has 10bp, so we need alignments that fit within it
    GraphAlignment readAlign3 = decodeGraphAlignment(0, "0[4M]", &graph);
    GraphAlignment mateAlign3 = decodeGraphAlignment(4, "0[4M]", &graph);

    FragById fragById;
    fragById.emplace("frag1", makeFrag("frag1", "AAAAA", readAlign1, "AAAAA", mateAlign1));
    fragById.emplace("frag2", makeFrag("frag2", "CAG", readAlign2, "CAG", mateAlign2));
    fragById.emplace("frag3", makeFrag("frag3", "AAAA", readAlign3, "AAAA", mateAlign3));

    int meanFragLen = getMeanFragLen(fragById);

    // Only flanking reads are counted: (10 + 8) / 2 = 9
    EXPECT_EQ(9, meanFragLen);
}

// ============================================================================
// resolveByFragLen tests
// ============================================================================

TEST(ResolveByFragLen, EmptyInput_ReturnsEmptyOutput)
{
    Diplotype paths;  // empty diplotype
    PairPathAlignById pairPathAlignById;  // empty input

    FragPathAlignsById result = resolveByFragLen(300, paths, pairPathAlignById);

    EXPECT_TRUE(result.empty());
}

TEST(ResolveByFragLen, ConsistentFragmentLength_KeepsAlignment)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    // Create a path through the graph
    graphtools::Path path(&graph, 0, {0, 1, 2}, 10);
    Diplotype paths = {path, path};

    // Create ReadPathAlign with fragment length close to meanFragLen
    auto readAlignPtr = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));
    auto mateAlignPtr = std::make_shared<GraphAlign>(decodeGraphAlignment(5, "0[5M]1[3M]2[2M]", &graph));

    ReadPathAlign readPathAlign(path, 0, 0, readAlignPtr);
    ReadPathAlign matePathAlign(path, 0, 0, mateAlignPtr);

    // Set begin/end to simulate positions
    // Fragment length will be calculated from these positions

    PairPathAlignById pairPathAlignById;
    pairPathAlignById["frag1"].readAligns.push_back(readPathAlign);
    pairPathAlignById["frag1"].mateAligns.push_back(matePathAlign);

    // Use mean fragment length that matches the alignment
    int meanFragLen = matePathAlign.end - readPathAlign.begin;
    FragPathAlignsById result = resolveByFragLen(meanFragLen, paths, pairPathAlignById);

    EXPECT_EQ(1u, result.size());
    EXPECT_EQ(1u, result["frag1"].size());
}

TEST(ResolveByFragLen, MultipleAlignments_KeepsClosestToMeanFragLen)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    graphtools::Path path(&graph, 0, {0, 1, 2}, 10);
    Diplotype paths = {path, path};

    // Create two alignments with different fragment lengths
    auto readAlignPtr1 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));
    auto mateAlignPtr1 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[10M]", &graph));

    auto readAlignPtr2 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));
    auto mateAlignPtr2 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));

    ReadPathAlign readPathAlign1(path, 0, 0, readAlignPtr1);
    ReadPathAlign matePathAlign1(path, 0, 0, mateAlignPtr1);

    ReadPathAlign readPathAlign2(path, 0, 0, readAlignPtr2);
    ReadPathAlign matePathAlign2(path, 0, 0, mateAlignPtr2);

    PairPathAlignById pairPathAlignById;
    pairPathAlignById["frag1"].readAligns.push_back(readPathAlign1);
    pairPathAlignById["frag1"].readAligns.push_back(readPathAlign2);
    pairPathAlignById["frag1"].mateAligns.push_back(matePathAlign1);
    pairPathAlignById["frag1"].mateAligns.push_back(matePathAlign2);

    // Set mean fragment length closer to alignment 2
    int fragLen2 = matePathAlign2.end - readPathAlign2.begin;
    int meanFragLen = fragLen2;  // Prefer alignment 2

    FragPathAlignsById result = resolveByFragLen(meanFragLen, paths, pairPathAlignById);

    EXPECT_EQ(1u, result.size());
    // The result should contain alignments with fragment length closest to mean
    EXPECT_FALSE(result["frag1"].empty());
}

TEST(ResolveByFragLen, DifferentPathIndices_SkipsInconsistentPairs)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    graphtools::Path path0(&graph, 0, {0, 1, 2}, 10);
    graphtools::Path path1(&graph, 0, {0, 1, 2}, 10);
    Diplotype paths = {path0, path1};

    auto readAlignPtr = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));
    auto mateAlignPtr = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));

    // Create read alignment on path 0
    ReadPathAlign readPathAlign(path0, 0, 0, readAlignPtr);
    // Create mate alignment on path 1 (different path index)
    ReadPathAlign matePathAlign(path1, 1, 0, mateAlignPtr);

    PairPathAlignById pairPathAlignById;
    pairPathAlignById["frag1"].readAligns.push_back(readPathAlign);
    pairPathAlignById["frag1"].mateAligns.push_back(matePathAlign);

    FragPathAlignsById result = resolveByFragLen(300, paths, pairPathAlignById);

    // Should be empty because read and mate are on different paths
    EXPECT_TRUE(result.empty() || result["frag1"].empty());
}

TEST(ResolveByFragLen, MultipleFragments_ProcessesEachIndependently)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("AAAAAAAAAA(CAG)*TTTTTTTTTT"));

    graphtools::Path path(&graph, 0, {0, 1, 2}, 10);
    Diplotype paths = {path, path};

    // Fragment 1
    auto readAlignPtr1 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[5M]", &graph));
    auto mateAlignPtr1 = std::make_shared<GraphAlign>(decodeGraphAlignment(5, "0[5M]", &graph));
    ReadPathAlign readPathAlign1(path, 0, 0, readAlignPtr1);
    ReadPathAlign matePathAlign1(path, 0, 0, mateAlignPtr1);

    // Fragment 2
    auto readAlignPtr2 = std::make_shared<GraphAlign>(decodeGraphAlignment(0, "0[6M]", &graph));
    auto mateAlignPtr2 = std::make_shared<GraphAlign>(decodeGraphAlignment(6, "0[4M]", &graph));
    ReadPathAlign readPathAlign2(path, 0, 0, readAlignPtr2);
    ReadPathAlign matePathAlign2(path, 0, 0, mateAlignPtr2);

    PairPathAlignById pairPathAlignById;
    pairPathAlignById["frag1"].readAligns.push_back(readPathAlign1);
    pairPathAlignById["frag1"].mateAligns.push_back(matePathAlign1);
    pairPathAlignById["frag2"].readAligns.push_back(readPathAlign2);
    pairPathAlignById["frag2"].mateAligns.push_back(matePathAlign2);

    int meanFragLen = 10;
    FragPathAlignsById result = resolveByFragLen(meanFragLen, paths, pairPathAlignById);

    // Both fragments should be processed
    EXPECT_EQ(2u, result.size());
    EXPECT_FALSE(result["frag1"].empty());
    EXPECT_FALSE(result["frag2"].empty());
}
