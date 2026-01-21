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

#include "reviewer/Projection.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"
#include "core/Read.hh"
#include "reviewer/Aligns.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using std::string;
using std::vector;

using namespace ehunter;
using namespace ehunter::reviewer;

// =============================================================================
// Tests for score() function
// =============================================================================

TEST(ReviewerProjectionScore, PerfectMatchAlignment_ReturnsPositiveScore)
{
    // Graph: TAAT(CCG)*CCTT
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // A perfect match alignment: 3 bases match in node 0, 3 bases match in node 1
    GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);

    // Default scoring: match=5, mismatch=-4, gap=-8
    // 6 matches * 5 = 30
    int alignScore = score(alignment);
    EXPECT_EQ(30, alignScore);
}

TEST(ReviewerProjectionScore, AlignmentWithMismatches_ReturnsLowerScore)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Alignment with 2 matches and 1 mismatch on node 0
    GraphAlignment alignment = decodeGraphAlignment(1, "0[2M1X]1[3M]", &graph);

    // 5 matches * 5 + 1 mismatch * (-4) = 25 - 4 = 21
    int alignScore = score(alignment);
    EXPECT_EQ(21, alignScore);
}

TEST(ReviewerProjectionScore, AlignmentWithInsertion_ReturnsLowerScore)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Alignment with insertion: 3M on node 0, then 2M1I on node 1
    // This means 3 matches on node 0, 2 matches + 1 insertion on node 1
    GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[2M1I]", &graph);

    // 5 matches * 5 + 1 insertion * (-8) = 25 - 8 = 17
    int alignScore = score(alignment);
    EXPECT_EQ(17, alignScore);
}

TEST(ReviewerProjectionScore, AlignmentWithDeletion_ReturnsLowerScore)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Alignment with deletion: 3M on node 0, then 2M1D on node 1
    GraphAlignment alignment = decodeGraphAlignment(1, "0[3M]1[2M1D]", &graph);

    // 5 matches * 5 + 1 deletion * (-8) = 25 - 8 = 17
    int alignScore = score(alignment);
    EXPECT_EQ(17, alignScore);
}

TEST(ReviewerProjectionScore, AlignmentWithCustomScoring_UsesProvidedScores)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    GraphAlignment alignment = decodeGraphAlignment(1, "0[2M1X]1[3M]", &graph);

    // Custom scoring: match=10, mismatch=-5, gap=-10
    // 5 matches * 10 + 1 mismatch * (-5) = 50 - 5 = 45
    int alignScore = score(alignment, 10, -5, -10);
    EXPECT_EQ(45, alignScore);
}

TEST(ReviewerProjectionScore, LongerAlignment_ScoreAccumulatesCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Longer alignment spanning multiple nodes
    GraphAlignment alignment = decodeGraphAlignment(0, "0[4M]1[3M]1[3M]2[4M]", &graph);

    // 14 matches * 5 = 70
    int alignScore = score(alignment);
    EXPECT_EQ(70, alignScore);
}

TEST(ReviewerProjectionScore, AlignmentWithMultipleOperations_ScoresCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Alignment with mix of operations: match, mismatch, insertion
    GraphAlignment alignment = decodeGraphAlignment(1, "0[2M1X]1[1M1I1M]", &graph);

    // 4 matches * 5 + 1 mismatch * (-4) + 1 insertion * (-8) = 20 - 4 - 8 = 8
    int alignScore = score(alignment);
    EXPECT_EQ(8, alignScore);
}

// =============================================================================
// Tests for project() function
// =============================================================================

// Helper function to create a fragment with reads and alignments
static FragById createFragById(
    const string& fragId,
    const GraphAlignment& readAlign,
    const GraphAlignment& mateAlign)
{
    // Create Read objects
    ReadId readId(fragId, MateNumber::kFirstMate);
    ReadId mateId(fragId, MateNumber::kSecondMate);

    // Generate dummy sequences of appropriate length
    string readSeq(readAlign.queryLength(), 'A');
    string mateSeq(mateAlign.queryLength(), 'A');

    Read read(readId, readSeq, false);
    Read mate(mateId, mateSeq, false);

    ReadWithAlign readWithAlign(std::move(read), readAlign);
    ReadWithAlign mateWithAlign(std::move(mate), mateAlign);

    Frag frag(std::move(readWithAlign), std::move(mateWithAlign));

    FragById fragById;
    fragById.emplace(fragId, std::move(frag));
    return fragById;
}

TEST(ReviewerProjectionProject, SingleFragmentOntoSinglePath_ProjectsCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Create alignments for read and mate
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[4M]1[3M]", &graph);
    GraphAlignment mateAlign = decodeGraphAlignment(0, "1[3M]2[4M]", &graph);

    FragById fragById = createFragById("frag1", readAlign, mateAlign);

    // Create a genotype path spanning the locus
    Path genotypePath(&graph, 0, {0, 1, 1, 2}, 4);
    vector<Path> genotypePaths = { genotypePath };

    PairPathAlignById result = project(genotypePaths, fragById);

    // Should have one fragment projected
    EXPECT_EQ(1u, result.size());
    EXPECT_TRUE(result.find("frag1") != result.end());

    // Both read and mate should have alignments
    const auto& pairAlign = result.at("frag1");
    EXPECT_FALSE(pairAlign.readAligns.empty());
    EXPECT_FALSE(pairAlign.mateAligns.empty());
}

TEST(ReviewerProjectionProject, MultipleFragments_AllProjected)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Create first fragment
    GraphAlignment readAlign1 = decodeGraphAlignment(0, "0[4M]1[3M]", &graph);
    GraphAlignment mateAlign1 = decodeGraphAlignment(0, "1[3M]2[4M]", &graph);
    FragById fragById1 = createFragById("frag1", readAlign1, mateAlign1);

    // Create second fragment
    GraphAlignment readAlign2 = decodeGraphAlignment(1, "0[3M]1[3M]", &graph);
    GraphAlignment mateAlign2 = decodeGraphAlignment(0, "1[3M]1[3M]2[2M]", &graph);
    FragById fragById2 = createFragById("frag2", readAlign2, mateAlign2);

    // Merge fragments
    FragById fragById;
    fragById.insert(fragById1.begin(), fragById1.end());
    fragById.insert(fragById2.begin(), fragById2.end());

    // Create genotype path
    Path genotypePath(&graph, 0, {0, 1, 1, 2}, 4);
    vector<Path> genotypePaths = { genotypePath };

    PairPathAlignById result = project(genotypePaths, fragById);

    // Should have both fragments projected
    EXPECT_EQ(2u, result.size());
    EXPECT_TRUE(result.find("frag1") != result.end());
    EXPECT_TRUE(result.find("frag2") != result.end());
}

TEST(ReviewerProjectionProject, TwoHaplotypePaths_FragmentsAlignToBest)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Create fragment that aligns to both paths
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[4M]1[3M]", &graph);
    GraphAlignment mateAlign = decodeGraphAlignment(0, "1[3M]2[4M]", &graph);

    FragById fragById = createFragById("frag1", readAlign, mateAlign);

    // Create two genotype paths (diplotype) with different repeat counts
    Path path1(&graph, 0, {0, 1, 2}, 4);        // 1 repeat
    Path path2(&graph, 0, {0, 1, 1, 1, 2}, 4);  // 3 repeats
    vector<Path> genotypePaths = { path1, path2 };

    PairPathAlignById result = project(genotypePaths, fragById);

    // Fragment should be projected
    EXPECT_EQ(1u, result.size());
    EXPECT_TRUE(result.find("frag1") != result.end());

    // Should have alignments
    const auto& pairAlign = result.at("frag1");
    EXPECT_FALSE(pairAlign.readAligns.empty());
    EXPECT_FALSE(pairAlign.mateAligns.empty());
}

TEST(ReviewerProjectionProject, EmptyFragmentMap_ReturnsEmptyResult)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    FragById emptyFragById;

    Path genotypePath(&graph, 0, {0, 1, 1, 2}, 4);
    vector<Path> genotypePaths = { genotypePath };

    PairPathAlignById result = project(genotypePaths, emptyFragById);

    EXPECT_TRUE(result.empty());
}

TEST(ReviewerProjectionProject, FragmentNotAlignableToPath_NotIncludedInResult)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    // Create a fragment where read aligns but mate does not overlap with path nodes
    // The read aligns to node 0, but the mate aligns only to node 2
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[4M]", &graph);
    GraphAlignment mateAlign = decodeGraphAlignment(0, "2[4M]", &graph);

    FragById fragById = createFragById("frag1", readAlign, mateAlign);

    // Create a path that only includes node 1 (repeat node)
    // This path won't connect the read (node 0) and mate (node 2)
    Path genotypePath(&graph, 0, {1, 1}, 3);
    vector<Path> genotypePaths = { genotypePath };

    PairPathAlignById result = project(genotypePaths, fragById);

    // Fragment cannot be properly projected because it doesn't align to the path
    // The result should be empty since the fragment can't align
    EXPECT_TRUE(result.empty());
}

TEST(ReviewerProjectionProject, PathIndexRecordedCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    GraphAlignment readAlign = decodeGraphAlignment(0, "0[4M]1[3M]", &graph);
    GraphAlignment mateAlign = decodeGraphAlignment(0, "1[3M]2[4M]", &graph);

    FragById fragById = createFragById("frag1", readAlign, mateAlign);

    // Create two paths
    Path path0(&graph, 0, {0, 1, 2}, 4);
    Path path1(&graph, 0, {0, 1, 1, 2}, 4);
    vector<Path> genotypePaths = { path0, path1 };

    PairPathAlignById result = project(genotypePaths, fragById);

    EXPECT_EQ(1u, result.size());

    // Check that path indices are recorded
    const auto& pairAlign = result.at("frag1");
    for (const auto& readPathAlign : pairAlign.readAligns)
    {
        EXPECT_GE(readPathAlign.pathIndex, 0);
        EXPECT_LT(readPathAlign.pathIndex, 2);
    }
}
