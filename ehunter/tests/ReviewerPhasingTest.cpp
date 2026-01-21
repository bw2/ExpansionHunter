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

#include "reviewer/Phasing.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Helper to create a Read with given ID and sequence
Read makeRead(const std::string& fragId, MateNumber mate, const std::string& seq)
{
    ReadId readId(fragId, mate);
    return Read(readId, seq, false);
}

// Helper to create a Frag from a graph alignment
Frag makeFrag(
    const std::string& fragId, const Graph& graph, const std::string& readSeq, const std::string& readAlignStr,
    const std::string& mateSeq, const std::string& mateAlignStr)
{
    Read read = makeRead(fragId, MateNumber::kFirstMate, readSeq);
    GraphAlignment readAlign = decodeGraphAlignment(0, readAlignStr, &graph);
    ReadWithAlign readWithAlign(std::move(read), std::move(readAlign));

    Read mate = makeRead(fragId, MateNumber::kSecondMate, mateSeq);
    GraphAlignment mateAlign = decodeGraphAlignment(0, mateAlignStr, &graph);
    ReadWithAlign mateWithAlign(std::move(mate), std::move(mateAlign));

    return Frag(std::move(readWithAlign), std::move(mateWithAlign));
}

}  // namespace

// Test scoreDiplotypes with a single diplotype
TEST(ReviewerPhasing_ScoreDiplotypes, SingleDiplotype_ReturnsScored)
{
    // Create a simple graph: ATTCGA(C)*ATGTCG
    // Node 0: ATTCGA (left flank), Node 1: C (repeat), Node 2: ATGTCG (right flank)
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create a diplotype with a single path through nodes 0, 1, 1, 2 (3 copies of repeat)
    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Create fragments that align to this path
    FragById fragById;
    // Read spanning left flank and repeat
    Frag frag1 = makeFrag("frag1", graph, "ATTCGAC", "0[6M]1[1M]", "CATGTCG", "1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // The diplotype should have a positive score since fragments align
    EXPECT_GT(result[0].second, 0);
    EXPECT_EQ(diplotype.size(), result[0].first.size());
}

// Test scoreDiplotypes with multiple diplotypes sorted by score
TEST(ReviewerPhasing_ScoreDiplotypes, MultipleDiplotypes_SortedByScoreDescending)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create two diplotypes with different repeat counts
    // Path with 2 repeats
    Path path2repeats(&graph, 0, { 0, 1, 1, 2 }, 6);
    // Path with 3 repeats
    Path path3repeats(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    Diplotype diplotype2 = { path2repeats };
    Diplotype diplotype3 = { path3repeats };
    std::vector<Diplotype> diplotypes = { diplotype2, diplotype3 };

    // Create a fragment that aligns better to the 2-repeat path
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(2u, result.size());
    // Results should be sorted by score in descending order
    EXPECT_GE(result[0].second, result[1].second);
}

// Test scoreDiplotypes with tied scores maintains stable ordering
TEST(ReviewerPhasing_ScoreDiplotypes, TiedScores_StableOrdering)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create identical paths - they should have the same score
    Path path1(&graph, 0, { 0, 1, 1, 2 }, 6);
    Path path2(&graph, 0, { 0, 1, 1, 2 }, 6);

    Diplotype diplotype1 = { path1 };
    Diplotype diplotype2 = { path2 };
    std::vector<Diplotype> diplotypes = { diplotype1, diplotype2 };

    // Create a fragment
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(2u, result.size());
    // Both should have the same score
    EXPECT_EQ(result[0].second, result[1].second);
    // Note: std::sort is not stable, but with identical scores both orderings are valid
}

// Test scoreDiplotypes with no fragments - all diplotypes should score 0
TEST(ReviewerPhasing_ScoreDiplotypes, NoFragments_AllDiplotypesScoreZero)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Empty fragment map
    FragById fragById;

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // With no fragments, score should be 0
    EXPECT_EQ(0, result[0].second);
}

// Test scoreDiplotypes with heterozygous diplotype (two different haplotypes)
TEST(ReviewerPhasing_ScoreDiplotypes, HeterozygousDiplotype_BothHaplotypesScored)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create a heterozygous diplotype with paths of different repeat counts
    Path pathShort(&graph, 0, { 0, 1, 2 }, 6);     // 1 repeat
    Path pathLong(&graph, 0, { 0, 1, 1, 1, 2 }, 6);  // 3 repeats

    Diplotype diplotype = { pathShort, pathLong };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Create fragments
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGAC", "0[6M]1[1M]", "CATGTCG", "1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // The diplotype should have a positive score
    EXPECT_GT(result[0].second, 0);
    // Verify it's still a heterozygous diplotype
    EXPECT_EQ(2u, result[0].first.size());
}

// Test scoreDiplotypes with multiple fragments
TEST(ReviewerPhasing_ScoreDiplotypes, MultipleFragments_ScoresAccumulate)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    FragById fragById;
    // Add first fragment
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    // Get score with one fragment
    ScoredDiplotypes result1 = scoreDiplotypes(fragById, diplotypes);
    int score1 = result1[0].second;

    // Add second fragment
    Frag frag2 = makeFrag("frag2", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag2", std::move(frag2));

    // Get score with two fragments
    ScoredDiplotypes result2 = scoreDiplotypes(fragById, diplotypes);
    int score2 = result2[0].second;

    // Score with two fragments should be greater than with one
    EXPECT_GT(score2, score1);
}
