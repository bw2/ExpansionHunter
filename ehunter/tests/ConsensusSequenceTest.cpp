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

#include "reviewer/ConsensusSequence.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using ehunter::MateNumber;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using namespace ehunter;
using namespace ehunter::reviewer;

// Helper function to create a FragAssignment from fragPathAlignsById
// Uses index 0 for each fragment (first alignment)
static FragAssignment createFragAssignment(const FragPathAlignsById& fragPathAlignsById)
{
    std::vector<std::string> fragIds;
    std::vector<int> alignIndexByFrag;
    for (const auto& pair : fragPathAlignsById)
    {
        fragIds.push_back(pair.first);
        alignIndexByFrag.push_back(0); // Use first alignment
    }
    return FragAssignment(std::move(fragIds), std::move(alignIndexByFrag));
}

// =============================================================================
// PositionVotes Tests
// =============================================================================

// Test that PositionVotes starts with no votes
TEST(PositionVotesTest, InitiallyEmpty)
{
    PositionVotes votes;
    EXPECT_FALSE(votes.hasVotes());
    EXPECT_EQ(0, votes.coverage);
    EXPECT_DOUBLE_EQ(0.0, votes.totalWeight());
}

// Test that addVote correctly accumulates weights for different bases
TEST(PositionVotesTest, AddVoteAccumulatesWeights)
{
    PositionVotes votes;

    // Add votes for A
    votes.addVote('A', 1.0);
    votes.addVote('A', 0.5);
    EXPECT_DOUBLE_EQ(1.5, votes.weights[0]); // A is at index 0

    // Add votes for C
    votes.addVote('C', 2.0);
    EXPECT_DOUBLE_EQ(2.0, votes.weights[1]); // C is at index 1

    // Check coverage count
    EXPECT_EQ(3, votes.coverage);
    EXPECT_TRUE(votes.hasVotes());
}

// Test that consensusBase returns the highest-voted base
TEST(PositionVotesTest, ConsensusBaseReturnsHighestVoted)
{
    PositionVotes votes;
    votes.addVote('A', 1.0);
    votes.addVote('C', 3.0);
    votes.addVote('G', 2.0);
    votes.addVote('T', 0.5);

    EXPECT_EQ('C', votes.consensusBase());
}

// Test consensusBase with ties (any of the tied bases is acceptable)
TEST(PositionVotesTest, ConsensusBaseHandlesTies)
{
    PositionVotes votes;
    votes.addVote('A', 2.0);
    votes.addVote('G', 2.0);

    char consensus = votes.consensusBase();
    // Either A or G is acceptable
    EXPECT_TRUE(consensus == 'A' || consensus == 'G');
}

// Test that totalWeight returns the sum of all weights
TEST(PositionVotesTest, TotalWeightSumsAllBases)
{
    PositionVotes votes;
    votes.addVote('A', 1.0);
    votes.addVote('C', 2.0);
    votes.addVote('G', 3.0);
    votes.addVote('T', 4.0);

    EXPECT_DOUBLE_EQ(10.0, votes.totalWeight());
}

// Test handling of lowercase bases
TEST(PositionVotesTest, HandlesLowercaseBases)
{
    PositionVotes votes;
    votes.addVote('a', 1.0);
    votes.addVote('c', 2.0);

    EXPECT_DOUBLE_EQ(1.0, votes.weights[0]); // 'a' -> A at index 0
    EXPECT_DOUBLE_EQ(2.0, votes.weights[1]); // 'c' -> C at index 1
    EXPECT_EQ(2, votes.coverage);
}

// Test that invalid bases are ignored
TEST(PositionVotesTest, IgnoresInvalidBases)
{
    PositionVotes votes;
    votes.addVote('N', 1.0); // N is not A, C, G, or T
    votes.addVote('X', 1.0); // X is not valid

    EXPECT_EQ(0, votes.coverage);
    EXPECT_DOUBLE_EQ(0.0, votes.totalWeight());
}

// Test that consensusBaseCount only counts reads matching the consensus base
TEST(PositionVotesTest, ConsensusBaseCountOnlyCountsMatchingReads)
{
    PositionVotes votes;

    // Add 5 votes for C (the consensus)
    votes.addVote('C', 1.0);
    votes.addVote('C', 1.0);
    votes.addVote('C', 1.0);
    votes.addVote('C', 1.0);
    votes.addVote('C', 1.0);

    // Add 2 votes for A (mismatches)
    votes.addVote('A', 1.0);
    votes.addVote('A', 1.0);

    // Add 1 vote for T (mismatch)
    votes.addVote('T', 1.0);

    // Total coverage is 8
    EXPECT_EQ(8, votes.coverage);

    // Consensus base should be C (highest count)
    EXPECT_EQ('C', votes.consensusBase());

    // consensusBaseCount should be 5 (only the C votes), NOT 8
    EXPECT_EQ(5, votes.consensusBaseCount());

    // Verify counts array
    EXPECT_EQ(2, votes.counts[0]); // A
    EXPECT_EQ(5, votes.counts[1]); // C
    EXPECT_EQ(0, votes.counts[2]); // G
    EXPECT_EQ(1, votes.counts[3]); // T
}

// Test that toReadSupportString uses consensusBaseCount, not coverage
TEST(AlleleConsensusTest, ReadSupportStringOnlyCountsMatchingReads)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 3;
    consensus.positions.resize(3);

    // Position 0: 5 votes for C, 2 votes for A
    consensus.positions[0] = PositionVotes();
    consensus.positions[0]->addVote('C', 1.0);
    consensus.positions[0]->addVote('C', 1.0);
    consensus.positions[0]->addVote('C', 1.0);
    consensus.positions[0]->addVote('C', 1.0);
    consensus.positions[0]->addVote('C', 1.0);
    consensus.positions[0]->addVote('A', 1.0); // mismatch
    consensus.positions[0]->addVote('A', 1.0); // mismatch

    // Position 1: 3 votes for A, 2 votes for G
    consensus.positions[1] = PositionVotes();
    consensus.positions[1]->addVote('A', 1.0);
    consensus.positions[1]->addVote('A', 1.0);
    consensus.positions[1]->addVote('A', 1.0);
    consensus.positions[1]->addVote('G', 1.0); // mismatch
    consensus.positions[1]->addVote('G', 1.0); // mismatch

    // Position 2: 4 votes for G, 1 vote for T
    consensus.positions[2] = PositionVotes();
    consensus.positions[2]->addVote('G', 1.0);
    consensus.positions[2]->addVote('G', 1.0);
    consensus.positions[2]->addVote('G', 1.0);
    consensus.positions[2]->addVote('G', 1.0);
    consensus.positions[2]->addVote('T', 1.0); // mismatch

    // Consensus sequence should be "CAG"
    std::string seqStr = consensus.toString();
    EXPECT_EQ("CAG", seqStr);

    // Read support string should be "534" (not "755" which would be total coverage)
    // Position 0: 5 reads matching C
    // Position 1: 3 reads matching A
    // Position 2: 4 reads matching G
    std::string supportStr = consensus.toReadSupportString();
    EXPECT_EQ("534", supportStr);

    // Verify that the coverage values are different from support values
    EXPECT_EQ(7, consensus.positions[0]->coverage); // 5 C + 2 A = 7 total
    EXPECT_EQ(5, consensus.positions[1]->coverage); // 3 A + 2 G = 5 total
    EXPECT_EQ(5, consensus.positions[2]->coverage); // 4 G + 1 T = 5 total
}

// =============================================================================
// AlleleConsensus Tests
// =============================================================================

// Test toString returns all N's for empty consensus
TEST(AlleleConsensusTest, EmptyConsensusReturnsAllNs)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 5;
    consensus.positions.resize(5); // All nullopt

    std::string result = consensus.toString();
    EXPECT_EQ("NNNNN", result);
    EXPECT_EQ(5u, result.length());
}

// Test toString returns correct bases for fully covered consensus
TEST(AlleleConsensusTest, FullyCoveredConsensusReturnsBases)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 4;
    consensus.positions.resize(4);

    // Add votes for each position
    consensus.positions[0] = PositionVotes();
    consensus.positions[0]->addVote('C', 1.0);

    consensus.positions[1] = PositionVotes();
    consensus.positions[1]->addVote('A', 1.0);

    consensus.positions[2] = PositionVotes();
    consensus.positions[2]->addVote('G', 1.0);

    consensus.positions[3] = PositionVotes();
    consensus.positions[3]->addVote('T', 1.0);

    std::string result = consensus.toString();
    EXPECT_EQ("CAGT", result);
}

// Test toString returns mixed bases and N's for partial coverage
TEST(AlleleConsensusTest, PartialCoverageReturnsMixedBasesAndNs)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 5;
    consensus.positions.resize(5);

    // Only cover positions 0 and 2
    consensus.positions[0] = PositionVotes();
    consensus.positions[0]->addVote('A', 1.0);

    consensus.positions[2] = PositionVotes();
    consensus.positions[2]->addVote('G', 1.0);

    std::string result = consensus.toString();
    EXPECT_EQ("ANGNN", result); // A at pos 0, N at pos 1, G at pos 2, N at pos 3 and 4
    EXPECT_EQ(5u, result.length());
}

// Test knownPositions counts correctly
TEST(AlleleConsensusTest, KnownPositionsCountsNonNPositions)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 5;
    consensus.positions.resize(5);

    // Add votes at positions 0, 2, and 3
    consensus.positions[0] = PositionVotes();
    consensus.positions[0]->addVote('A', 1.0);

    consensus.positions[2] = PositionVotes();
    consensus.positions[2]->addVote('G', 1.0);

    consensus.positions[3] = PositionVotes();
    consensus.positions[3]->addVote('T', 1.0);

    EXPECT_EQ(3, consensus.knownPositions());
}

// Test knownPositions returns 0 for empty consensus
TEST(AlleleConsensusTest, KnownPositionsReturnsZeroForEmpty)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 5;
    consensus.positions.resize(5);

    EXPECT_EQ(0, consensus.knownPositions());
}

// =============================================================================
// AnchorClassification Tests
// =============================================================================

// Test AnchorClassification default values
TEST(AnchorClassificationTest, DefaultValues)
{
    AnchorClassification classification;

    EXPECT_FALSE(classification.isAnchor);
    EXPECT_EQ(-1, classification.assignedAllele);
}

// =============================================================================
// RepeatAlignedBases Tests
// =============================================================================

// Test RepeatAlignedBases default values
TEST(RepeatAlignedBasesTest, DefaultValues)
{
    RepeatAlignedBases rab;

    EXPECT_TRUE(rab.sequence.empty());
    EXPECT_TRUE(rab.weights.empty());
    EXPECT_EQ(0, rab.repeatStartPos);
    EXPECT_FALSE(rab.valid);
}

// =============================================================================
// ConsensusResult Tests
// =============================================================================

// Test ConsensusResult default values
TEST(ConsensusResultTest, DefaultValues)
{
    ConsensusResult result;

    EXPECT_TRUE(result.alleleConsensuses.empty());
    EXPECT_EQ(0, result.totalAnchors);
    EXPECT_EQ(0, result.totalFragments);
}

// =============================================================================
// computeRepeatRegionLength Tests
// =============================================================================

// Test computing repeat region length for a path with repeat nodes
TEST(ComputeRepeatRegionLengthTest, PathWithRepeatNodes)
{
    // Create graph: "ATTCGA(CAG)*ATGTCG"
    // Node 0: ATTCGA (6 bases) - left flank
    // Node 1: CAG (3 bases, repeat)
    // Node 2: ATGTCG (6 bases) - right flank
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(CAG)*ATGTCG"));

    // Path with 3 repeat units: 0 -> 1 -> 1 -> 1 -> 2
    Path path(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    int repeatLength = computeRepeatRegionLength(path);
    EXPECT_EQ(9, repeatLength); // 3 units * 3 bases per unit = 9
}

// Test computing repeat region length for a path with no repeat nodes
TEST(ComputeRepeatRegionLengthTest, PathWithNoRepeatNodes)
{
    // Create a simple linear graph with no repeats
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGAATGTCG"));

    // Path through linear region
    Path path(&graph, 0, { 0 }, 12);

    int repeatLength = computeRepeatRegionLength(path);
    EXPECT_EQ(0, repeatLength);
}

// =============================================================================
// AnchorClassification Tests
// =============================================================================

// Helper function to create a ReadWithAlign object
static ReadWithAlign makeReadWithAlign(
    const std::string& fragId, MateNumber mateNum, const std::string& sequence, const Graph& graph,
    const std::string& alignEncoding)
{
    ehunter::ReadId readId(fragId, mateNum);
    ehunter::Read read(readId, sequence, false);
    GraphAlignment align = decodeGraphAlignment(0, alignEncoding, &graph);
    return ReadWithAlign(std::move(read), std::move(align));
}

// Helper function to create a Frag object
static Frag makeFragment(
    const std::string& fragId, const std::string& readSeq, const std::string& mateSeq, const Graph& graph,
    const std::string& readAlignEncoding, const std::string& mateAlignEncoding)
{
    return Frag(
        makeReadWithAlign(fragId, MateNumber::kFirstMate, readSeq, graph, readAlignEncoding),
        makeReadWithAlign(fragId, MateNumber::kSecondMate, mateSeq, graph, mateAlignEncoding));
}

// Test that a fragment with alignments to a single allele (single path) is classified as anchor
TEST(ClassifyAsAnchorTest, SinglePathIsAnchor)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with two paths
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12);

    // Create a fragment with alignment to only one allele (path 0)
    auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]", &graph));

    ReadPathAlign readPathAlign(path0, 0, 0, readAlign);
    ReadPathAlign matePathAlign(path0, 0, 1, mateAlign);
    FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"].push_back(std::move(fragPathAlign));
    // Note: Only one alignment - fragment aligns to only one allele

    // Classify - should be anchor since all alignments have same pathIndex
    AnchorClassification result = classifyAsAnchor("frag1", fragPathAlignsById);

    EXPECT_TRUE(result.isAnchor);
    EXPECT_EQ(0, result.assignedAllele);
}

// Test that a fragment with alignments to multiple alleles (multiple paths) is NOT an anchor
TEST(ClassifyAsAnchorTest, MultiplePathsIsNotAnchor)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with two paths
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12);

    // Create alignments to BOTH alleles
    auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

    auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

    ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
    ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
    FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

    ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
    ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
    FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"].push_back(std::move(fragPathAlign0));
    fragPathAlignsById["frag1"].push_back(std::move(fragPathAlign1));

    // Classify - should NOT be anchor since alignments have different pathIndex values
    AnchorClassification result = classifyAsAnchor("frag1", fragPathAlignsById);

    EXPECT_FALSE(result.isAnchor);
    EXPECT_EQ(-1, result.assignedAllele);
}

// Test that a fragment with no alignments is not an anchor
TEST(ClassifyAsAnchorTest, NoAlignmentsIsNotAnchor)
{
    FragPathAlignsById fragPathAlignsById;
    // Empty alignments

    AnchorClassification result = classifyAsAnchor("frag1", fragPathAlignsById);

    EXPECT_FALSE(result.isAnchor);
    EXPECT_EQ(-1, result.assignedAllele);
}

// Test that a fragment with multiple alignments to the same allele IS an anchor
TEST(ClassifyAsAnchorTest, MultipleAlignmentsSamePathIsAnchor)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);

    // Create two different alignments, both to path 0
    auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

    auto readAlign2 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[10M]", &graph));
    auto mateAlign2 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]", &graph));

    ReadPathAlign readPathAlign1(path0, 0, 0, readAlign1);
    ReadPathAlign matePathAlign1(path0, 0, 1, mateAlign1);
    FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

    ReadPathAlign readPathAlign2(path0, 0, 0, readAlign2);
    ReadPathAlign matePathAlign2(path0, 0, 1, mateAlign2);
    FragPathAlign fragPathAlign2(std::move(readPathAlign2), std::move(matePathAlign2));

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"].push_back(std::move(fragPathAlign1));
    fragPathAlignsById["frag1"].push_back(std::move(fragPathAlign2));

    // Classify - should be anchor since all alignments have pathIndex=0
    AnchorClassification result = classifyAsAnchor("frag1", fragPathAlignsById);

    EXPECT_TRUE(result.isAnchor);
    EXPECT_EQ(0, result.assignedAllele);
}

// =============================================================================
// Integration Tests: buildConsensusFromAnchors
// =============================================================================

// Helper function to create complete fragment and alignment data for integration tests
static void createFragmentWithAlignment(
    const std::string& fragId, const std::string& readSeq, const std::string& mateSeq, const Graph& graph,
    const std::string& readAlignEncoding, const std::string& mateAlignEncoding, const Path& haplotypePath,
    int pathIndex, int readStartIdx, int mateStartIdx, FragById& fragById, FragPathAlignsById& fragPathAlignsById)
{
    // Create the fragment
    Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, readAlignEncoding, mateAlignEncoding);
    fragById.emplace(fragId, std::move(frag));

    // Create projected alignments
    auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, readAlignEncoding, &graph));
    auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, mateAlignEncoding, &graph));

    ReadPathAlign readPathAlign(haplotypePath, pathIndex, readStartIdx, readAlign);
    ReadPathAlign matePathAlign(haplotypePath, pathIndex, mateStartIdx, mateAlign);
    FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

    fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
}

// Test buildConsensusFromAnchors with a simple homozygous case (single allele)
TEST(BuildConsensusFromAnchorsTest, HomozygousSingleAllele)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    // Node 0: ATTCGATTCGA (11 bases) - left flank
    // Node 1: CAG (3 bases, repeat)
    // Node 2: ATGTCGATGTCG (12 bases) - right flank
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with single path (homozygous, 3 repeat units = 9 bp repeat region)
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create anchor fragment 1: spans left flank and all 3 repeats
    // Read: 11 bp left flank + 9 bp repeat (CAGCAGCAG)
    // Mate: 12 bp right flank
    createFragmentWithAlignment(
        "frag1",
        "ATTCGATTCGACAGCAGCAG", // Full left flank + 3 repeat units
        "ATGTCGATGTCG", // Full right flank
        graph,
        "0[11M]1[3M]1[3M]1[3M]", // Read alignment
        "2[12M]", // Mate alignment
        path0, 0, 0, 4, // pathIndex=0, read starts at node 0, mate starts at node 4
        fragById, fragPathAlignsById);

    // Create anchor fragment 2: spans left flank and first 2 repeats
    // Read: 11 bp left flank + 6 bp repeat (CAGCAG)
    // Mate: overlaps last repeat + right flank
    createFragmentWithAlignment(
        "frag2",
        "ATTCGATTCGACAGCAG", // Full left flank + 2 repeat units
        "CAGATGTCGATGTCG", // 1 repeat unit + full right flank
        graph,
        "0[11M]1[3M]1[3M]", // Read alignment
        "1[3M]2[12M]", // Mate alignment
        path0, 0, 0, 3, // pathIndex=0, read starts at node 0, mate starts at node 3
        fragById, fragPathAlignsById);

    // Build consensus with minimum flank overlap of 10
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Verify we have 1 allele consensus
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    // Verify the consensus for allele 0
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength); // 3 units * 3 bp = 9 bp

    // The consensus should be "CAGCAGCAG" (or with some N's if coverage is incomplete)
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ(9u, consensusStr.length());

    // First 6 bases should be covered by both fragments (positions 0-5)
    // We expect "CAGCAG" for the first 6 positions
    EXPECT_EQ('C', consensusStr[0]);
    EXPECT_EQ('A', consensusStr[1]);
    EXPECT_EQ('G', consensusStr[2]);
    EXPECT_EQ('C', consensusStr[3]);
    EXPECT_EQ('A', consensusStr[4]);
    EXPECT_EQ('G', consensusStr[5]);

    // Last 3 bases (positions 6-8) should be covered by frag1 and frag2's mate
    EXPECT_EQ('C', consensusStr[6]);
    EXPECT_EQ('A', consensusStr[7]);
    EXPECT_EQ('G', consensusStr[8]);

    // Verify anchor counts
    EXPECT_EQ(2, result.totalAnchors);
    EXPECT_EQ(2, result.totalFragments);
    EXPECT_GE(consensus0.anchorReadCount, 1); // At least 1 anchor contributed
}

// Test buildConsensusFromAnchors with heterozygous case (two alleles of different sizes)
TEST(BuildConsensusFromAnchorsTest, HeterozygousTwoAlleles)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create heterozygous diplotype:
    // Allele 0: 2 repeat units (6 bp repeat region)
    // Allele 1: 4 repeat units (12 bp repeat region)
    Path path0(&graph, 0, { 0, 1, 1, 2 }, 12); // 2 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12); // 4 repeat units
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create anchor fragment for shorter allele (allele 0):
    // Read spans left flank and 2 repeats, mate spans right flank
    createFragmentWithAlignment(
        "frag_short",
        "ATTCGATTCGACAGCAG", // Full left flank + 2 repeat units
        "ATGTCGATGTCG", // Full right flank
        graph,
        "0[11M]1[3M]1[3M]", // Read alignment
        "2[12M]", // Mate alignment
        path0, 0, 0, 3, // pathIndex=0 (shorter allele)
        fragById, fragPathAlignsById);

    // Create anchor fragment for longer allele (allele 1):
    // This fragment ONLY aligns to path1 (we don't add it to path0's alignments)
    // Read spans left flank and 4 repeats (extends beyond shorter allele)
    {
        std::string fragId = "frag_long";
        std::string readSeq = "ATTCGATTCGACAGCAGCAGCAG"; // Full left flank + 4 repeat units
        std::string mateSeq = "ATGTCGATGTCG"; // Full right flank
        std::string readAlignEnc = "0[11M]1[3M]1[3M]1[3M]1[3M]";
        std::string mateAlignEnc = "2[12M]";

        // Create the fragment
        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, readAlignEnc, mateAlignEnc);
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignment ONLY for path1 (longer allele)
        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, readAlignEnc, &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, mateAlignEnc, &graph));

        ReadPathAlign readPathAlign(path1, 1, 0, readAlign); // pathIndex=1
        ReadPathAlign matePathAlign(path1, 1, 5, mateAlign); // after 4 repeat nodes
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Build consensus with minimum flank overlap of 10
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Verify we have 2 allele consensuses
    ASSERT_EQ(2u, result.alleleConsensuses.size());

    // Verify consensus for allele 0 (shorter, 2 repeats = 6 bp)
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(6, consensus0.repeatLength); // 2 units * 3 bp = 6 bp

    std::string consensus0Str = consensus0.toString();
    EXPECT_EQ(6u, consensus0Str.length());
    // Should be "CAGCAG"
    EXPECT_EQ('C', consensus0Str[0]);
    EXPECT_EQ('A', consensus0Str[1]);
    EXPECT_EQ('G', consensus0Str[2]);
    EXPECT_EQ('C', consensus0Str[3]);
    EXPECT_EQ('A', consensus0Str[4]);
    EXPECT_EQ('G', consensus0Str[5]);

    // Verify consensus for allele 1 (longer, 4 repeats = 12 bp)
    const AlleleConsensus& consensus1 = result.alleleConsensuses[1];
    EXPECT_EQ(1, consensus1.alleleIndex);
    EXPECT_EQ(12, consensus1.repeatLength); // 4 units * 3 bp = 12 bp

    std::string consensus1Str = consensus1.toString();
    EXPECT_EQ(12u, consensus1Str.length());
    // Should be "CAGCAGCAGCAG"
    EXPECT_EQ('C', consensus1Str[0]);
    EXPECT_EQ('A', consensus1Str[1]);
    EXPECT_EQ('G', consensus1Str[2]);
    EXPECT_EQ('C', consensus1Str[3]);
    EXPECT_EQ('A', consensus1Str[4]);
    EXPECT_EQ('G', consensus1Str[5]);
    EXPECT_EQ('C', consensus1Str[6]);
    EXPECT_EQ('A', consensus1Str[7]);
    EXPECT_EQ('G', consensus1Str[8]);
    EXPECT_EQ('C', consensus1Str[9]);
    EXPECT_EQ('A', consensus1Str[10]);
    EXPECT_EQ('G', consensus1Str[11]);

    // Verify total anchors
    EXPECT_EQ(2, result.totalAnchors);
}

// Test buildConsensusFromAnchors with no anchor reads (all fragments align to multiple paths)
TEST(BuildConsensusFromAnchorsTest, NoAnchorReads)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with two paths (heterozygous)
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12); // 4 repeat units
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create a fragment that aligns to BOTH paths (ambiguous)
    // With the new classification logic, ambiguous fragments are NOT anchors
    {
        std::string fragId = "frag_ambiguous";
        std::string readSeq = "ATTCGATTCGACAGCAG"; // left flank + 2 repeat units (fits on both alleles)
        std::string mateSeq = "ATGTCGATGTCG"; // right flank

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignments to BOTH alleles
        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
        ReadPathAlign matePathAlign0(path0, 0, 3, mateAlign0);
        FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

        ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
        ReadPathAlign matePathAlign1(path1, 1, 4, mateAlign1);
        FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign0));
        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign1));
    }

    // Build consensus
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Verify we have 2 allele consensuses
    ASSERT_EQ(2u, result.alleleConsensuses.size());

    // Both consensuses should be all N's since no anchors qualified (all fragments ambiguous)
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength); // 3 units * 3 bp = 9 bp
    std::string consensusStr0 = consensus0.toString();
    EXPECT_EQ("NNNNNNNNN", consensusStr0);
    EXPECT_EQ(0, consensus0.knownPositions());
    EXPECT_EQ(0, consensus0.anchorReadCount);

    const AlleleConsensus& consensus1 = result.alleleConsensuses[1];
    EXPECT_EQ(1, consensus1.alleleIndex);
    EXPECT_EQ(12, consensus1.repeatLength); // 4 units * 3 bp = 12 bp
    std::string consensusStr1 = consensus1.toString();
    EXPECT_EQ("NNNNNNNNNNNN", consensusStr1);
    EXPECT_EQ(0, consensus1.knownPositions());
    EXPECT_EQ(0, consensus1.anchorReadCount);

    // No anchors should be counted
    EXPECT_EQ(0, result.totalAnchors);
    EXPECT_EQ(1, result.totalFragments);
}

// Test buildConsensusFromAnchors with repeat interruption (variant in repeat)
TEST(BuildConsensusFromAnchorsTest, RepeatInterruption)
{
    // Create graph with CAG repeat
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with single path (3 repeat units = 9 bp)
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create fragment 1: normal CAG repeat
    createFragmentWithAlignment(
        "frag1",
        "ATTCGATTCGACAGCAGCAG", // Full left flank + 3 CAG
        "ATGTCGATGTCG", graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]", path0, 0, 0, 4, fragById, fragPathAlignsById);

    // Create fragment 2: with a CAA interruption at position 4-5 (middle of second unit)
    // The sequence "CAGCAACAG" has the middle CAG changed to CAA
    // Note: The alignment will still work with mismatches
    createFragmentWithAlignment(
        "frag2",
        "ATTCGATTCGACAGCAACAG", // Full left flank + CAG + CAA + CAG (interruption!)
        "ATGTCGATGTCG", graph,
        "0[11M]1[3M]1[3M]1[3M]", // Will have a mismatch at position 4
        "2[12M]", path0, 0, 0, 4, fragById, fragPathAlignsById);

    // Build consensus
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ(9u, consensusStr.length());

    // Position 0: C from both fragments -> C
    EXPECT_EQ('C', consensusStr[0]);
    // Position 1: A from both fragments -> A
    EXPECT_EQ('A', consensusStr[1]);
    // Position 2: G from both fragments -> G
    EXPECT_EQ('G', consensusStr[2]);
    // Position 3: C from both fragments -> C
    EXPECT_EQ('C', consensusStr[3]);
    // Position 4: A from frag1, A from frag2 -> A (both have A here!)
    EXPECT_EQ('A', consensusStr[4]);
    // Position 5: G from frag1, A from frag2 -> G or A (tie possible, G wins due to order or weight)
    // With equal weights, either G or A is acceptable
    EXPECT_TRUE(consensusStr[5] == 'G' || consensusStr[5] == 'A');

    EXPECT_EQ(2, result.totalAnchors);
}

// Test buildConsensusFromAnchors with empty input
TEST(BuildConsensusFromAnchorsTest, EmptyInput)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with single path
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0 };

    // Empty fragment maps
    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Build consensus with no fragments
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Should still produce an allele consensus (with all N's)
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength);
    EXPECT_EQ("NNNNNNNNN", consensus0.toString());
    EXPECT_EQ(0, consensus0.knownPositions());
    EXPECT_EQ(0, consensus0.anchorReadCount);
    EXPECT_EQ(0, result.totalAnchors);
    EXPECT_EQ(0, result.totalFragments);
}

// =============================================================================
// RepeatAlignedBases Deletion Handling Tests
// =============================================================================

// Test that extractRepeatAlignedBases correctly maps positions when there is a deletion inside the repeat
// This verifies that bases after a deletion land on the correct repeat coordinates (no left shift)
TEST(RepeatAlignedBasesTest, DeletionInRepeatPreservesPositions)
{
    // Create graph: "ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"
    // Node 0: ATTCGATTCGA (11 bases) - left flank
    // Node 1: CAGCAG (6 bases, repeat unit)
    // Node 2: ATGTCGATGTCG (12 bases) - right flank
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"));

    // Create a path with 2 repeat units (12 bp repeat region total)
    // Path: node 0 -> node 1 -> node 1 -> node 2
    Path haplotypePath(&graph, 0, { 0, 1, 1, 2 }, 12);

    // Create a read that spans the first repeat unit with a deletion
    // The read has sequence "CACAG" which would align as:
    // Reference: CAGCAG (repeat unit)
    // Read:      CA-CAG (2M1D3M - 2 match, 1 deletion of G, 3 match)
    // Positions: 0,1, ,3,4,5 (position 2 is deleted, no base)
    //
    // With the old bug: bases after deletion would shift left to positions 0,1,2,3,4
    // With the fix: bases should be at positions 0,1,3,4,5 (position 2 is skipped)
    std::string readSeq = "CACAG"; // 5 bases that align to 6bp repeat unit with 1 deletion

    // Create the read with alignment containing a deletion in the repeat
    ehunter::ReadId readId("frag_del", MateNumber::kFirstMate);
    ehunter::Read read(readId, readSeq, false);

    // Graph alignment: starts at position 0 of node 1, has CIGAR 2M1D3M
    // This means: 2 matches (CA), 1 deletion (missing G), 3 matches (CAG)
    GraphAlignment graphAlign = decodeGraphAlignment(0, "1[2M1D3M]", &graph);
    ReadWithAlign readWithAlign(std::move(read), std::move(graphAlign));

    // Create projected alignment onto the haplotype path
    // The read aligns to the first repeat node (index 1 in the path)
    auto projectedAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[2M1D3M]", &graph));

    // pathAlign.startIndexOnPath = 1 means alignment starts at the second node in path (first repeat node)
    // pathAlign.begin should be 11 (start of repeat region on linearized path, after 11bp left flank)
    // pathAlign.end should be 11 + 6 = 17 (end of first repeat unit)
    ReadPathAlign pathAlign(haplotypePath, 0, 1, projectedAlign);

    // Extract the repeat-aligned bases
    RepeatAlignedBases result = extractRepeatAlignedBases(readWithAlign, pathAlign, haplotypePath);

    // Verify extraction was successful
    EXPECT_TRUE(result.valid);

    // Verify the extracted sequence
    EXPECT_EQ("CACAG", result.sequence);
    EXPECT_EQ(5u, result.sequence.size());

    // KEY TEST: Verify the repeat positions are correct (no left shift after deletion)
    // The extracted bases should map to repeat coordinates 0, 1, 3, 4, 5
    // NOT to coordinates 0, 1, 2, 3, 4 (which would be wrong - left shifted)
    ASSERT_EQ(5u, result.repeatPositions.size());
    EXPECT_EQ(0, result.repeatPositions[0]); // 'C' at repeat position 0
    EXPECT_EQ(1, result.repeatPositions[1]); // 'A' at repeat position 1
    // Position 2 is deleted (no base in read)
    EXPECT_EQ(3, result.repeatPositions[2]); // 'C' at repeat position 3
    EXPECT_EQ(4, result.repeatPositions[3]); // 'A' at repeat position 4
    EXPECT_EQ(5, result.repeatPositions[4]); // 'G' at repeat position 5

    // Verify repeatStartPos is set to the first position
    EXPECT_EQ(0, result.repeatStartPos);
}

// Test that extractRepeatAlignedBases handles a deletion at the end of a repeat unit
TEST(RepeatAlignedBasesTest, DeletionAtEndOfRepeatUnit)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Path with 3 repeat units (9 bp repeat region)
    Path haplotypePath(&graph, 0, { 0, 1, 1, 1, 2 }, 12);

    // Create a read that spans 2 repeat units with a deletion at end of first unit
    // Read sequence "CACAGCAG" aligns as:
    // Reference: CAGCAGCAG (3 repeat units)
    // Read:      CA-CAGCAG (2M1D6M - 2 match, 1 deletion of G, 6 match)
    // Positions: 0,1, ,3,4,5,6,7,8
    std::string readSeq = "CACAGCAG"; // 8 bases that align to 9bp repeat with 1 deletion

    ehunter::ReadId readId("frag_del2", MateNumber::kFirstMate);
    ehunter::Read read(readId, readSeq, false);

    // Alignment: 2M1D6M across two repeat nodes
    // Node 1 (first): 2M1D (covers 3bp reference)
    // Node 1 (second): 3M (covers full 3bp)
    // Node 1 (third): 3M (covers full 3bp)
    GraphAlignment graphAlign = decodeGraphAlignment(0, "1[2M1D]1[3M]1[3M]", &graph);
    ReadWithAlign readWithAlign(std::move(read), std::move(graphAlign));

    auto projectedAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[2M1D]1[3M]1[3M]", &graph));
    ReadPathAlign pathAlign(haplotypePath, 0, 1, projectedAlign);

    RepeatAlignedBases result = extractRepeatAlignedBases(readWithAlign, pathAlign, haplotypePath);

    EXPECT_TRUE(result.valid);
    EXPECT_EQ("CACAGCAG", result.sequence);
    EXPECT_EQ(8u, result.sequence.size());

    // Verify positions: deletion at position 2
    ASSERT_EQ(8u, result.repeatPositions.size());
    EXPECT_EQ(0, result.repeatPositions[0]); // 'C' at position 0
    EXPECT_EQ(1, result.repeatPositions[1]); // 'A' at position 1
    // Position 2 is deleted
    EXPECT_EQ(3, result.repeatPositions[2]); // 'C' at position 3
    EXPECT_EQ(4, result.repeatPositions[3]); // 'A' at position 4
    EXPECT_EQ(5, result.repeatPositions[4]); // 'G' at position 5
    EXPECT_EQ(6, result.repeatPositions[5]); // 'C' at position 6
    EXPECT_EQ(7, result.repeatPositions[6]); // 'A' at position 7
    EXPECT_EQ(8, result.repeatPositions[7]); // 'G' at position 8
}

// Test that consensus building correctly uses position-accurate mapping with deletions
TEST(BuildConsensusFromAnchorsTest, DeletionDoesNotShiftConsensusPositions)
{
    // Create graph: "ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"
    // Node 0: ATTCGATTCGA (11 bases) - left flank
    // Node 1: CAGCAG (6 bases, repeat unit)
    // Node 2: ATGTCGATGTCG (12 bases) - right flank
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"));

    // Create diplotype with single path (1 repeat unit = 6 bp repeat region)
    Path path0(&graph, 0, { 0, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create anchor fragment 1: perfect alignment (no deletion)
    // Read: 11 bp left flank + 6 bp repeat (CAGCAG)
    // Mate: 12 bp right flank
    createFragmentWithAlignment(
        "frag_perfect",
        "ATTCGATTCGACAGCAG", // Full left flank + full repeat unit
        "ATGTCGATGTCG", // Full right flank
        graph,
        "0[11M]1[6M]", // Read alignment
        "2[12M]", // Mate alignment
        path0, 0, 0, 2, // pathIndex=0, read starts at node 0, mate starts at node 2
        fragById, fragPathAlignsById);

    // Create anchor fragment 2: alignment with deletion at position 2 (missing G)
    // The read covers left flank + repeat with a deletion
    // Read sequence: ATTCGATTCGA + CACAG (left flank + 5bp repeat with deletion)
    // This should contribute votes to positions 0, 1, 3, 4, 5 (NOT 0, 1, 2, 3, 4)
    {
        std::string fragId = "frag_deletion";
        std::string readSeq = "ATTCGATTCGACACAG"; // 11bp flank + 5bp repeat (deletion of position 2)
        std::string mateSeq = "ATGTCGATGTCG"; // Full right flank

        ehunter::ReadId readIdR(fragId, MateNumber::kFirstMate);
        ehunter::Read readR(readIdR, readSeq, false);
        // CIGAR: 11M on left flank, then 2M1D3M on repeat (2 match, 1 del, 3 match)
        GraphAlignment readGraphAlign = decodeGraphAlignment(0, "0[11M]1[2M1D3M]", &graph);

        ehunter::ReadId mateIdR(fragId, MateNumber::kSecondMate);
        ehunter::Read mateR(mateIdR, mateSeq, false);
        GraphAlignment mateGraphAlign = decodeGraphAlignment(0, "2[12M]", &graph);

        ReadWithAlign readWithAlign(std::move(readR), std::move(readGraphAlign));
        ReadWithAlign mateWithAlign(std::move(mateR), std::move(mateGraphAlign));

        Frag frag(std::move(readWithAlign), std::move(mateWithAlign));
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignments
        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[2M1D3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path0, 0, 0, readAlign);
        ReadPathAlign matePathAlign(path0, 0, 2, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Build consensus with minimum flank overlap of 10
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Verify we have 1 allele consensus
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(6, consensus0.repeatLength); // 1 unit * 6 bp = 6 bp

    // The consensus should be "CAGCAG"
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ(6u, consensusStr.length());

    // Position 0: Both fragments vote 'C' -> C
    EXPECT_EQ('C', consensusStr[0]);
    // Position 1: Both fragments vote 'A' -> A
    EXPECT_EQ('A', consensusStr[1]);
    // Position 2: Only frag_perfect votes 'G' (frag_deletion has deletion here)
    // Should still be 'G' since we have at least one vote
    EXPECT_EQ('G', consensusStr[2]);
    // Position 3: Both fragments vote 'C' -> C
    EXPECT_EQ('C', consensusStr[3]);
    // Position 4: Both fragments vote 'A' -> A
    EXPECT_EQ('A', consensusStr[4]);
    // Position 5: Both fragments vote 'G' -> G
    EXPECT_EQ('G', consensusStr[5]);

    // Verify the full consensus is correct
    EXPECT_EQ("CAGCAG", consensusStr);

    // Both anchors should be counted
    EXPECT_EQ(2, result.totalAnchors);
    EXPECT_EQ(2, result.totalFragments);

    // Verify read support - position 2 should have lower support (only 1 read vs 2 for others)
    std::string supportStr = consensus0.toReadSupportString();
    EXPECT_EQ(6u, supportStr.length());

    // Positions 0, 1, 3, 4, 5 should have support from both reads (2)
    // Position 2 should have support from only 1 read (the one without deletion)
    EXPECT_EQ('2', supportStr[0]); // Position 0: 2 reads
    EXPECT_EQ('2', supportStr[1]); // Position 1: 2 reads
    EXPECT_EQ('1', supportStr[2]); // Position 2: 1 read (other has deletion)
    EXPECT_EQ('2', supportStr[3]); // Position 3: 2 reads
    EXPECT_EQ('2', supportStr[4]); // Position 4: 2 reads
    EXPECT_EQ('2', supportStr[5]); // Position 5: 2 reads
}

// =============================================================================
// Homozygous Genotype Tests
// =============================================================================

// Test buildConsensusFromAnchors with homozygous genotype represented as two identical-length paths
// This is the common case in practice: diplotype.size() == 2 but both paths have the same repeat length
// The fix ensures we only produce ONE consensus sequence with all anchors combined
TEST(BuildConsensusFromAnchorsTest, HomozygousTwoIdenticalPaths)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with TWO paths of the SAME repeat length (both have 3 repeats)
    // This simulates a homozygous 3/3 genotype
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units (identical to path0)
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create anchor fragment 1: aligned to path0
    createFragmentWithAlignment(
        "frag1",
        "ATTCGATTCGACAGCAGCAG", // Full left flank + 3 CAG repeats
        "ATGTCGATGTCG", // Full right flank
        graph,
        "0[11M]1[3M]1[3M]1[3M]", // Read alignment
        "2[12M]", // Mate alignment
        path0, 0, 0, 4, // pathIndex=0
        fragById, fragPathAlignsById);

    // Create anchor fragment 2: aligned to path1
    // (In the old buggy behavior, this anchor would be lost/unused)
    createFragmentWithAlignment(
        "frag2",
        "ATTCGATTCGACAGCAGCAG", // Full left flank + 3 CAG repeats
        "ATGTCGATGTCG", // Full right flank
        graph,
        "0[11M]1[3M]1[3M]1[3M]", // Read alignment
        "2[12M]", // Mate alignment
        path1, 1, 0, 4, // pathIndex=1 (second allele)
        fragById, fragPathAlignsById);

    // Build consensus with minimum flank overlap of 10
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // KEY ASSERTION: Should produce only ONE consensus (not two with one all-N's)
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength); // 3 units * 3 bp = 9 bp

    // Consensus should be "CAGCAGCAG"
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ(9u, consensusStr.length());
    EXPECT_EQ("CAGCAGCAG", consensusStr);

    // BOTH anchors should be counted (all contributing to the single consensus)
    EXPECT_EQ(2, result.totalAnchors);
    EXPECT_EQ(2, result.totalFragments);

    // The consensus should have been built from both fragments
    EXPECT_GE(consensus0.anchorReadCount, 1);
}

// Test that a fragment aligned ONLY to allele 1 (not allele 0) still contributes to the single
// consensus in a homozygous genotype (equal-length alleles).
// This tests the fix for the anchor filtering bug where anchors classified for allele 1
// were being dropped even in homozygous mode.
TEST(BuildConsensusFromAnchorsTest, HomozygousAllele1OnlyAnchorContributes)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with TWO paths of the SAME repeat length (both have 3 repeats)
    // This simulates a homozygous 3/3 genotype
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units (identical to path0)
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create anchor fragment that ONLY aligns to allele 1 (NOT allele 0)
    // This simulates the case where the aligner only produced an alignment to path1
    // In a homozygous genotype, this should still contribute to the consensus
    {
        std::string fragId = "frag_allele1_only";
        std::string readSeq = "ATTCGATTCGACAGCAGCAG"; // Full left flank + 3 CAG repeats
        std::string mateSeq = "ATGTCGATGTCG"; // Full right flank

        // Create the fragment
        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignment ONLY to path1 (allele 1), NOT to path0
        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]1[3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path1, 1, 0, readAlign); // pathIndex=1 (allele 1 only!)
        ReadPathAlign matePathAlign(path1, 1, 4, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        // Note: We're NOT adding any alignment to path0/allele 0
        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Build consensus with minimum flank overlap of 10
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Should produce only ONE consensus (homozygous)
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength); // 3 units * 3 bp = 9 bp

    // KEY ASSERTION: The fragment should have contributed to the consensus
    // even though it only aligned to allele 1.
    // With the fix, classification.assignedAllele is overridden to 0 for homozygous,
    // so the fragment passes the filter in buildAlleleConsensus.
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ(9u, consensusStr.length());
    EXPECT_EQ("CAGCAGCAG", consensusStr);

    // The anchor should be counted
    EXPECT_EQ(1, result.totalAnchors);
    EXPECT_EQ(1, result.totalFragments);

    // The single anchor should have contributed to the consensus
    EXPECT_EQ(1, consensus0.anchorReadCount);

    // Verify read support - all positions should have support of 1
    std::string supportStr = consensus0.toReadSupportString();
    EXPECT_EQ("111111111", supportStr);
}

// =============================================================================
// Ambiguous Anchor Exclusion Tests
// =============================================================================

// Test that ambiguous fragments (aligning equally well to both alleles) are NOT included
// in totalAnchors count when building consensus.
// This verifies the policy that ambiguous fragments are excluded from consensus building entirely.
TEST(BuildConsensusFromAnchorsTest, AmbiguousFragmentsNotCountedInTotalAnchors)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create heterozygous diplotype with similar-length alleles
    // Both alleles have enough repeats that the fragment fits on both
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12); // 4 repeat units
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create one AMBIGUOUS fragment that aligns equally well to both alleles
    // (2 repeat units covered, both alleles have >= 2 repeats)
    {
        std::string fragId = "frag_ambiguous";
        std::string readSeq = "ATTCGATTCGA"; // Full left flank
        std::string mateSeq = "CAGCAG"; // 2 repeat units (fits on both alleles)

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]", "1[3M]1[3M]");
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignments to BOTH alleles
        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
        ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
        FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

        ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
        ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
        FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign0));
        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign1));
    }

    // Create one UNAMBIGUOUS fragment that ONLY aligns to allele 0
    {
        std::string fragId = "frag_unambiguous";
        std::string readSeq = "ATTCGATTCGA"; // Full left flank
        std::string mateSeq = "CAGCAGCAG"; // 3 repeat units

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]", "1[3M]1[3M]1[3M]");
        fragById.emplace(fragId, std::move(frag));

        // Create projected alignment ONLY to allele 0
        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]1[3M]", &graph));

        ReadPathAlign readPathAlign(path0, 0, 0, readAlign);
        ReadPathAlign matePathAlign(path0, 0, 1, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
        // Note: No alignment to allele 1
    }

    // Build consensus
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // KEY ASSERTION: totalAnchors should be 1 (only the unambiguous fragment)
    // The ambiguous fragment should NOT be counted
    EXPECT_EQ(1, result.totalAnchors);
    EXPECT_EQ(2, result.totalFragments);

    // Verify we have 2 allele consensuses (heterozygous)
    ASSERT_EQ(2u, result.alleleConsensuses.size());

    // Allele 0 should have the unambiguous fragment's contribution
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(9, consensus0.repeatLength); // 3 units * 3 bp
    // The unambiguous fragment should have contributed
    std::string consensusStr0 = consensus0.toString();
    EXPECT_EQ(9u, consensusStr0.length());

    // Allele 1 should have no anchors (the ambiguous one was excluded)
    const AlleleConsensus& consensus1 = result.alleleConsensuses[1];
    EXPECT_EQ(12, consensus1.repeatLength); // 4 units * 3 bp
    // All N's because the ambiguous fragment was excluded
    EXPECT_EQ("NNNNNNNNNNNN", consensus1.toString());
}

// Test that ambiguous fragment classification is order-independent.
// When the same fragment has alignments to both alleles, the ConsensusResult should be
// identical regardless of which order the alignments are presented.
TEST(BuildConsensusFromAnchorsTest, AmbiguousFragmentOrderIndependent)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create heterozygous diplotype
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12); // 4 repeat units
    std::vector<Path> diplotype = { path0, path1 };

    // Helper lambda to create fragment data
    auto createFragmentData = [&](const std::string& fragId, bool allele0First) {
        FragById fragById;
        FragPathAlignsById fragPathAlignsById;

        std::string readSeq = "ATTCGATTCGA";
        std::string mateSeq = "CAGCAG";

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]", "1[3M]1[3M]");
        fragById.emplace(fragId, std::move(frag));

        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
        ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
        FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

        ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
        ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
        FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

        if (allele0First)
        {
            fragPathAlignsById[fragId].push_back(std::move(fragPathAlign0));
            fragPathAlignsById[fragId].push_back(std::move(fragPathAlign1));
        }
        else
        {
            fragPathAlignsById[fragId].push_back(std::move(fragPathAlign1));
            fragPathAlignsById[fragId].push_back(std::move(fragPathAlign0));
        }

        return std::make_pair(std::move(fragById), std::move(fragPathAlignsById));
    };

    // Build consensus with allele 0 alignment first
    auto [fragById0First, fragPathAlignsById0First] = createFragmentData("frag_ambiguous", true);
    ConsensusResult result0First = buildConsensusFromAnchors(diplotype, fragById0First, fragPathAlignsById0First, createFragAssignment(fragPathAlignsById0First));

    // Build consensus with allele 1 alignment first
    auto [fragById1First, fragPathAlignsById1First] = createFragmentData("frag_ambiguous", false);
    ConsensusResult result1First = buildConsensusFromAnchors(diplotype, fragById1First, fragPathAlignsById1First, createFragAssignment(fragPathAlignsById1First));

    // KEY ASSERTIONS: Results should be IDENTICAL regardless of order
    EXPECT_EQ(result0First.totalAnchors, result1First.totalAnchors);
    EXPECT_EQ(result0First.totalFragments, result1First.totalFragments);

    // Both should have 0 anchors (ambiguous fragment excluded)
    EXPECT_EQ(0, result0First.totalAnchors);
    EXPECT_EQ(0, result1First.totalAnchors);

    // Both should have 2 allele consensuses
    ASSERT_EQ(result0First.alleleConsensuses.size(), result1First.alleleConsensuses.size());
    ASSERT_EQ(2u, result0First.alleleConsensuses.size());

    // Consensus sequences should be identical (all N's since no anchors)
    for (size_t i = 0; i < result0First.alleleConsensuses.size(); ++i)
    {
        EXPECT_EQ(
            result0First.alleleConsensuses[i].toString(), result1First.alleleConsensuses[i].toString())
            << "Consensus for allele " << i << " differs based on alignment order";
        EXPECT_EQ(
            result0First.alleleConsensuses[i].repeatLength, result1First.alleleConsensuses[i].repeatLength);
        EXPECT_EQ(
            result0First.alleleConsensuses[i].anchorReadCount, result1First.alleleConsensuses[i].anchorReadCount);
    }
}

// Test that classifyAsAnchor returns the correct values for ambiguous fragments
// A fragment with alignments to multiple paths (different pathIndex values) is ambiguous
TEST(ClassifyAsAnchorTest, AmbiguousFragmentHasCorrectClassificationValues)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with two paths
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12);

    std::string fragId = "frag_ambiguous";

    // Create projected alignments to BOTH alleles (different pathIndex values)
    auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

    auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
    auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

    ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
    ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
    FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

    ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
    ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
    FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById[fragId].push_back(std::move(fragPathAlign0));
    fragPathAlignsById[fragId].push_back(std::move(fragPathAlign1));

    // Classify the fragment - should NOT be anchor since it aligns to multiple paths
    AnchorClassification result = classifyAsAnchor(fragId, fragPathAlignsById);

    EXPECT_FALSE(result.isAnchor) << "Ambiguous fragment should NOT be an anchor";
    EXPECT_EQ(-1, result.assignedAllele) << "Ambiguous fragment should have assignedAllele = -1";
}

// Test that classifying the same ambiguous fragment with different alignment orders
// produces identical classification results
TEST(ClassifyAsAnchorTest, AmbiguousClassificationOrderIndependent)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12);

    std::string fragId = "frag_ambiguous";

    // Setup 1: Allele 0 alignment first
    FragPathAlignsById fragPathAlignsById1;
    {
        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));
        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
        ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
        FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

        ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
        ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
        FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

        fragPathAlignsById1[fragId].push_back(std::move(fragPathAlign0)); // Allele 0 first
        fragPathAlignsById1[fragId].push_back(std::move(fragPathAlign1));
    }

    // Setup 2: Allele 1 alignment first
    FragPathAlignsById fragPathAlignsById2;
    {
        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));
        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));

        ReadPathAlign readPathAlign0(path0, 0, 0, readAlign0);
        ReadPathAlign matePathAlign0(path0, 0, 1, mateAlign0);
        FragPathAlign fragPathAlign0(std::move(readPathAlign0), std::move(matePathAlign0));

        ReadPathAlign readPathAlign1(path1, 1, 0, readAlign1);
        ReadPathAlign matePathAlign1(path1, 1, 1, mateAlign1);
        FragPathAlign fragPathAlign1(std::move(readPathAlign1), std::move(matePathAlign1));

        fragPathAlignsById2[fragId].push_back(std::move(fragPathAlign1)); // Allele 1 first
        fragPathAlignsById2[fragId].push_back(std::move(fragPathAlign0));
    }

    // Classify using each setup - should give same result regardless of order
    AnchorClassification result1 = classifyAsAnchor(fragId, fragPathAlignsById1);
    AnchorClassification result2 = classifyAsAnchor(fragId, fragPathAlignsById2);

    // Both should produce the same classification
    EXPECT_EQ(result1.isAnchor, result2.isAnchor);
    EXPECT_EQ(result1.assignedAllele, result2.assignedAllele);

    // Both should be excluded as ambiguous
    EXPECT_FALSE(result1.isAnchor);
    EXPECT_FALSE(result2.isAnchor);
    EXPECT_EQ(-1, result1.assignedAllele);
    EXPECT_EQ(-1, result2.assignedAllele);
}

// Test that multiple fragments, each aligning to different alleles (one to 0, one to 1),
// both contribute to the single consensus in a homozygous genotype.
// This is a more comprehensive test of the homozygous anchor retention fix.
TEST(BuildConsensusFromAnchorsTest, HomozygousMixedAlleleAnchorsAllContribute)
{
    // Create graph: "ATTCGATTCGA(CAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with TWO paths of the SAME repeat length (both have 3 repeats)
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units
    Path path1(&graph, 0, { 0, 1, 1, 1, 2 }, 12); // 3 repeat units (identical)
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Fragment 1: ONLY aligns to allele 0
    {
        std::string fragId = "frag_allele0";
        std::string readSeq = "ATTCGATTCGACAGCAGCAG";
        std::string mateSeq = "ATGTCGATGTCG";

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]1[3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path0, 0, 0, readAlign); // pathIndex=0 only
        ReadPathAlign matePathAlign(path0, 0, 4, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Fragment 2: ONLY aligns to allele 1
    {
        std::string fragId = "frag_allele1";
        std::string readSeq = "ATTCGATTCGACAGCAGCAG";
        std::string mateSeq = "ATGTCGATGTCG";

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]1[3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path1, 1, 0, readAlign); // pathIndex=1 only
        ReadPathAlign matePathAlign(path1, 1, 4, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Fragment 3: ONLY aligns to allele 1 (another one, to verify multiple allele-1 anchors work)
    {
        std::string fragId = "frag_allele1_b";
        std::string readSeq = "ATTCGATTCGACAGCAGCAG";
        std::string mateSeq = "ATGTCGATGTCG";

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]1[3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path1, 1, 0, readAlign); // pathIndex=1 only
        ReadPathAlign matePathAlign(path1, 1, 4, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Build consensus
    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Should produce only ONE consensus (homozygous)
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];

    // ALL THREE fragments should have contributed, regardless of their original allele assignment
    EXPECT_EQ(3, result.totalAnchors);
    EXPECT_EQ(3, result.totalFragments);

    // Consensus should be correct
    std::string consensusStr = consensus0.toString();
    EXPECT_EQ("CAGCAGCAG", consensusStr);

    // All 3 anchors should have contributed to the consensus
    EXPECT_EQ(3, consensus0.anchorReadCount);

    // Verify read support - all positions should have support of 3
    std::string supportStr = consensus0.toReadSupportString();
    EXPECT_EQ("333333333", supportStr);
}

// =============================================================================
// Integration Tests: Combined Behavior Tests
// =============================================================================

// Integration test: Verify JSON output format compatibility
// This test ensures that toString() and toReadSupportString() produce valid
// output that can be directly used in JSON arrays.
TEST(IntegrationTest, JsonOutputFieldsFormat)
{
    // Create graph
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Create diplotype with single path (homozygous, 3 repeat units)
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Create multiple anchors with varying coverage
    createFragmentWithAlignment(
        "frag1",
        "ATTCGATTCGACAGCAGCAG", // Full coverage
        "ATGTCGATGTCG", graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]", path0, 0, 0, 4, fragById, fragPathAlignsById);

    createFragmentWithAlignment(
        "frag2",
        "ATTCGATTCGACAGCAGCAG", // Full coverage
        "ATGTCGATGTCG", graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]", path0, 0, 0, 4, fragById, fragPathAlignsById);

    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    ASSERT_EQ(1u, result.alleleConsensuses.size());
    const AlleleConsensus& consensus = result.alleleConsensuses[0];

    // Verify ConsensusSequences format: should be a simple DNA string
    std::string seqStr = consensus.toString();
    EXPECT_EQ(9u, seqStr.length());
    // Should contain only valid DNA bases (A, C, G, T) or N for unknown
    for (char c : seqStr)
    {
        EXPECT_TRUE(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N')
            << "Invalid character '" << c << "' in consensus sequence";
    }

    // Verify ConsensusSequencesReadSupport format: should be digits only
    std::string supportStr = consensus.toReadSupportString();
    EXPECT_EQ(seqStr.length(), supportStr.length());
    for (char c : supportStr)
    {
        EXPECT_TRUE(c >= '0' && c <= '9') << "Invalid character '" << c << "' in read support string";
    }

    // Verify the actual values
    EXPECT_EQ("CAGCAGCAG", seqStr);
    EXPECT_EQ("222222222", supportStr); // 2 reads at each position
}

// Integration test: Deletion mapping with homozygous genotype
// This tests that deletions are handled correctly in homozygous mode
TEST(IntegrationTest, DeletionMappingInHomozygousGenotype)
{
    // Create graph: "ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAGCAG)*ATGTCGATGTCG"));

    // Create homozygous diplotype (two identical paths with 1 repeat unit = 6bp)
    Path path0(&graph, 0, { 0, 1, 2 }, 12);
    Path path1(&graph, 0, { 0, 1, 2 }, 12);
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Fragment 1: perfect alignment on allele 0
    createFragmentWithAlignment(
        "frag_perfect_a0",
        "ATTCGATTCGACAGCAG", // Full left flank + full repeat unit
        "ATGTCGATGTCG", graph, "0[11M]1[6M]", "2[12M]", path0, 0, 0, 2, fragById, fragPathAlignsById);

    // Fragment 2: alignment with deletion on allele 1
    {
        std::string fragId = "frag_deletion_a1";
        std::string readSeq = "ATTCGATTCGACACAG"; // 11bp flank + 5bp repeat (deletion of position 2)
        std::string mateSeq = "ATGTCGATGTCG";

        ehunter::ReadId readIdR(fragId, MateNumber::kFirstMate);
        ehunter::Read readR(readIdR, readSeq, false);
        GraphAlignment readGraphAlign = decodeGraphAlignment(0, "0[11M]1[2M1D3M]", &graph);

        ehunter::ReadId mateIdR(fragId, MateNumber::kSecondMate);
        ehunter::Read mateR(mateIdR, mateSeq, false);
        GraphAlignment mateGraphAlign = decodeGraphAlignment(0, "2[12M]", &graph);

        ReadWithAlign readWithAlign(std::move(readR), std::move(readGraphAlign));
        ReadWithAlign mateWithAlign(std::move(mateR), std::move(mateGraphAlign));

        Frag frag(std::move(readWithAlign), std::move(mateWithAlign));
        fragById.emplace(fragId, std::move(frag));

        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[2M1D3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        // Aligned to path1 (allele 1) - but should still contribute in homozygous mode
        ReadPathAlign readPathAlign(path1, 1, 0, readAlign);
        ReadPathAlign matePathAlign(path1, 1, 2, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Should produce only ONE consensus (homozygous)
    ASSERT_EQ(1u, result.alleleConsensuses.size());

    const AlleleConsensus& consensus = result.alleleConsensuses[0];

    // Both fragments should have contributed
    EXPECT_EQ(2, result.totalAnchors);

    // Consensus should still be correct (CAGCAG)
    std::string consensusStr = consensus.toString();
    EXPECT_EQ("CAGCAG", consensusStr);

    // Read support should show:
    // - Position 2 has support from only 1 read (the one without deletion)
    // - All other positions have support from 2 reads
    std::string supportStr = consensus.toReadSupportString();
    EXPECT_EQ(6u, supportStr.length());
    EXPECT_EQ('2', supportStr[0]); // Position 0: 2 reads
    EXPECT_EQ('2', supportStr[1]); // Position 1: 2 reads
    EXPECT_EQ('1', supportStr[2]); // Position 2: 1 read (deletion in other)
    EXPECT_EQ('2', supportStr[3]); // Position 3: 2 reads
    EXPECT_EQ('2', supportStr[4]); // Position 4: 2 reads
    EXPECT_EQ('2', supportStr[5]); // Position 5: 2 reads
}

// Integration test: Heterozygous genotype with one allele having ambiguous reads
// Verifies that ambiguous reads don't contaminate the other allele's consensus
TEST(IntegrationTest, HeterozygousWithAmbiguousAndUnambiguousReads)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Heterozygous: 3 repeats vs 5 repeats
    Path path0(&graph, 0, { 0, 1, 1, 1, 2 }, 12);              // 3 repeats (9bp)
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 1, 2 }, 12);        // 5 repeats (15bp)
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Unambiguous anchor for allele 0 (only aligns to shorter allele)
    createFragmentWithAlignment(
        "frag_unambig_a0",
        "ATTCGATTCGACAGCAGCAG", // Covers all 3 repeats + right flank
        "ATGTCGATGTCG", graph, "0[11M]1[3M]1[3M]1[3M]", "2[12M]", path0, 0, 0, 4, fragById, fragPathAlignsById);

    // Unambiguous anchor for allele 1 (length-discriminating, extends beyond shorter allele)
    {
        std::string fragId = "frag_unambig_a1";
        std::string readSeq = "ATTCGATTCGACAGCAGCAGCAGCAG"; // 11bp flank + 15bp repeat (5 units)
        std::string mateSeq = "ATGTCGATGTCG";

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]1[3M]1[3M]1[3M]1[3M]1[3M]", "2[12M]");
        fragById.emplace(fragId, std::move(frag));

        auto readAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]1[3M]1[3M]1[3M]1[3M]1[3M]", &graph));
        auto mateAlign = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "2[12M]", &graph));

        ReadPathAlign readPathAlign(path1, 1, 0, readAlign);
        ReadPathAlign matePathAlign(path1, 1, 6, mateAlign);
        FragPathAlign fragPathAlign(std::move(readPathAlign), std::move(matePathAlign));

        fragPathAlignsById[fragId].push_back(std::move(fragPathAlign));
    }

    // Ambiguous fragment (aligns to both alleles, covers only 2 repeats)
    {
        std::string fragId = "frag_ambiguous";
        std::string readSeq = "ATTCGATTCGA"; // Left flank only
        std::string mateSeq = "CAGCAG";       // 2 repeats - fits on both alleles

        Frag frag = makeFragment(fragId, readSeq, mateSeq, graph, "0[11M]", "1[3M]1[3M]");
        fragById.emplace(fragId, std::move(frag));

        // Alignment to allele 0
        auto readAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign0 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));
        ReadPathAlign rpa0(path0, 0, 0, readAlign0);
        ReadPathAlign mpa0(path0, 0, 1, mateAlign0);
        FragPathAlign fpa0(std::move(rpa0), std::move(mpa0));

        // Alignment to allele 1
        auto readAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "0[11M]", &graph));
        auto mateAlign1 = std::make_shared<GraphAlignment>(decodeGraphAlignment(0, "1[3M]1[3M]", &graph));
        ReadPathAlign rpa1(path1, 1, 0, readAlign1);
        ReadPathAlign mpa1(path1, 1, 1, mateAlign1);
        FragPathAlign fpa1(std::move(rpa1), std::move(mpa1));

        fragPathAlignsById[fragId].push_back(std::move(fpa0));
        fragPathAlignsById[fragId].push_back(std::move(fpa1));
    }

    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    // Should have 2 consensus sequences (heterozygous)
    ASSERT_EQ(2u, result.alleleConsensuses.size());

    // Only 2 anchors (the two unambiguous ones), NOT 3
    EXPECT_EQ(2, result.totalAnchors);
    EXPECT_EQ(3, result.totalFragments); // 3 fragments examined

    // Verify allele 0 consensus (3 repeats)
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(0, consensus0.alleleIndex);
    EXPECT_EQ(9, consensus0.repeatLength);
    EXPECT_EQ("CAGCAGCAG", consensus0.toString());
    EXPECT_EQ("111111111", consensus0.toReadSupportString()); // 1 unambiguous anchor

    // Verify allele 1 consensus (5 repeats)
    const AlleleConsensus& consensus1 = result.alleleConsensuses[1];
    EXPECT_EQ(1, consensus1.alleleIndex);
    EXPECT_EQ(15, consensus1.repeatLength);
    EXPECT_EQ("CAGCAGCAGCAGCAG", consensus1.toString());
    EXPECT_EQ("111111111111111", consensus1.toReadSupportString()); // 1 unambiguous anchor
}

// Integration test: Read support capping at 9
// This tests that coverage values >= 10 are displayed as '9' in the support string
TEST(IntegrationTest, ReadSupportCappedAtNine)
{
    AlleleConsensus consensus;
    consensus.alleleIndex = 0;
    consensus.repeatLength = 3;
    consensus.positions.resize(3);

    // Position 0: 5 reads
    consensus.positions[0] = PositionVotes();
    for (int i = 0; i < 5; ++i)
    {
        consensus.positions[0]->addVote('C', 1.0);
    }

    // Position 1: 10 reads (should cap at '9')
    consensus.positions[1] = PositionVotes();
    for (int i = 0; i < 10; ++i)
    {
        consensus.positions[1]->addVote('A', 1.0);
    }

    // Position 2: 15 reads (should cap at '9')
    consensus.positions[2] = PositionVotes();
    for (int i = 0; i < 15; ++i)
    {
        consensus.positions[2]->addVote('G', 1.0);
    }

    std::string seqStr = consensus.toString();
    EXPECT_EQ("CAG", seqStr);

    std::string supportStr = consensus.toReadSupportString();
    EXPECT_EQ(3u, supportStr.length());
    EXPECT_EQ('5', supportStr[0]); // 5 reads
    EXPECT_EQ('9', supportStr[1]); // 10 reads -> capped to 9
    EXPECT_EQ('9', supportStr[2]); // 15 reads -> capped to 9
}

// Integration test: Empty/uncovered positions in heterozygous genotype
// Verifies that uncovered positions produce 'N' in sequence and '0' in support
TEST(IntegrationTest, PartialCoverageProducesCorrectOutput)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGATTCGA(CAG)*ATGTCGATGTCG"));

    // Heterozygous: 2 repeats vs 4 repeats
    Path path0(&graph, 0, { 0, 1, 1, 2 }, 12);              // 2 repeats (6bp)
    Path path1(&graph, 0, { 0, 1, 1, 1, 1, 2 }, 12);        // 4 repeats (12bp)
    std::vector<Path> diplotype = { path0, path1 };

    FragById fragById;
    FragPathAlignsById fragPathAlignsById;

    // Anchor for shorter allele (allele 0)
    createFragmentWithAlignment(
        "frag_short",
        "ATTCGATTCGACAGCAG", // 11bp flank + 6bp repeat
        "ATGTCGATGTCG", graph, "0[11M]1[3M]1[3M]", "2[12M]", path0, 0, 0, 3, fragById, fragPathAlignsById);

    // No anchor for longer allele - it should have all N's

    ConsensusResult result = buildConsensusFromAnchors(diplotype, fragById, fragPathAlignsById, createFragAssignment(fragPathAlignsById));

    ASSERT_EQ(2u, result.alleleConsensuses.size());

    // Allele 0 should have full coverage
    const AlleleConsensus& consensus0 = result.alleleConsensuses[0];
    EXPECT_EQ(6, consensus0.repeatLength);
    EXPECT_EQ("CAGCAG", consensus0.toString());
    EXPECT_EQ("111111", consensus0.toReadSupportString());

    // Allele 1 should have no coverage (all N's, all 0's)
    const AlleleConsensus& consensus1 = result.alleleConsensuses[1];
    EXPECT_EQ(12, consensus1.repeatLength);
    EXPECT_EQ("NNNNNNNNNNNN", consensus1.toString());
    EXPECT_EQ("000000000000", consensus1.toReadSupportString());
}
