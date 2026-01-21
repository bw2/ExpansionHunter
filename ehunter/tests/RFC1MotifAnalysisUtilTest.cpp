//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

#include "locus/RFC1MotifAnalysisUtil.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using namespace ehunter;

TEST(RFC1MotifAnalysisTests, MeanTest)
{
    std::vector<double> x { 10, 3, 4, 5, 10 };
    EXPECT_DOUBLE_EQ(4.0, mean(x.begin() + 1, x.end() - 1));
}

TEST(RFC1MotifAnalysisTests, MinRotationTest)
{
    EXPECT_EQ("AAGGC", getMinRotation("GGCAA"));
    EXPECT_EQ("GGGGT", getMinRotation("GGGGT"));
}

TEST(RFC1MotifAnalysisTests, FindUsableBaseRange_Test)
{
    // Test fwd orientation
    //
    {
        std::vector<uint8_t> binaryQuals(15, 1);
        auto baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 14u), *baseRange);

        binaryQuals[13] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 12u), *baseRange);

        binaryQuals[1] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, false);
        EXPECT_FALSE(baseRange);
    }

    // Test rev orientation
    //
    {
        std::vector<uint8_t> binaryQuals(15, 1);
        auto baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(0u, 14u), *baseRange);

        binaryQuals[1] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_TRUE(baseRange);
        EXPECT_EQ(std::make_pair(2u, 14u), *baseRange);

        binaryQuals[14] = 0;
        baseRange = findUsableReadBaseRange(binaryQuals, true);
        EXPECT_FALSE(baseRange);

        // Test that method didn't alter binaryQuals:
        EXPECT_EQ(0, binaryQuals[1]);
    }
}

// Tests for extractRFC1ReadAlignments
//
// Creates a graph with regex "TAAT(CCG)*CCTT" where:
// - Node 0: left flank "TAAT"
// - Node 1: repeat unit "CCG"
// - Node 2: right flank "CCTT"

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_EmptyBuffer_ReturnsEmpty)
{
    locus::AlignmentBuffer buffer;
    auto result = extractRFC1ReadAlignments(buffer);
    EXPECT_TRUE(result.empty());
}

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_SingleFragmentOverlappingNode1_ReturnsOneAlignment)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    locus::AlignmentBuffer buffer;
    locus::AlignedFragment fragment{
        "frag1",
        "TAATCCGCCG",
        decodeGraphAlignment(0, "0[4M]1[3M]1[3M]", &graph),
        true,
        "CCGCCTT",
        decodeGraphAlignment(0, "1[3M]2[4M]", &graph),
        false
    };

    buffer.tryInsertFragment(fragment);

    auto result = extractRFC1ReadAlignments(buffer);
    ASSERT_EQ(1u, result.size());
    EXPECT_EQ("TAATCCGCCG", result[0].read);
    EXPECT_FALSE(result[0].isReversed);  // isReversed = !readIsForwardStrand = !true = false
}

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_NeitherReadNorMateOverlapsNode1_ReturnsEmpty)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    locus::AlignmentBuffer buffer;
    locus::AlignedFragment fragment{
        "frag1",
        "TAAT",
        decodeGraphAlignment(0, "0[4M]", &graph),  // Only node 0
        true,
        "CCTT",
        decodeGraphAlignment(0, "2[4M]", &graph),  // Only node 2
        false
    };

    buffer.tryInsertFragment(fragment);

    auto result = extractRFC1ReadAlignments(buffer);
    EXPECT_TRUE(result.empty());
}

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_OnlyMateOverlapsNode1_ReturnsEmpty)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    locus::AlignmentBuffer buffer;
    locus::AlignedFragment fragment{
        "frag1",
        "TAAT",
        decodeGraphAlignment(0, "0[4M]", &graph),  // Only node 0
        true,
        "CCGCCTT",
        decodeGraphAlignment(0, "1[3M]2[4M]", &graph),  // Overlaps node 1
        false
    };

    buffer.tryInsertFragment(fragment);

    // Only the read is extracted, not the mate, so if read doesn't overlap node 1, nothing returned
    auto result = extractRFC1ReadAlignments(buffer);
    EXPECT_TRUE(result.empty());
}

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_MultipleFragments_OnlyOverlappingReturned)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    locus::AlignmentBuffer buffer;

    // Fragment 1: read overlaps node 1
    locus::AlignedFragment frag1{
        "frag1",
        "TAATCCG",
        decodeGraphAlignment(0, "0[4M]1[3M]", &graph),
        true,
        "CCTT",
        decodeGraphAlignment(0, "2[4M]", &graph),
        false
    };
    buffer.tryInsertFragment(frag1);

    // Fragment 2: read does not overlap node 1
    locus::AlignedFragment frag2{
        "frag2",
        "TAAT",
        decodeGraphAlignment(0, "0[4M]", &graph),
        false,
        "CCGCCTT",
        decodeGraphAlignment(0, "1[3M]2[4M]", &graph),
        true
    };
    buffer.tryInsertFragment(frag2);

    // Fragment 3: read overlaps node 1
    locus::AlignedFragment frag3{
        "frag3",
        "CCGCCGCCTT",
        decodeGraphAlignment(0, "1[3M]1[3M]2[4M]", &graph),
        false,
        "TAAT",
        decodeGraphAlignment(0, "0[4M]", &graph),
        true
    };
    buffer.tryInsertFragment(frag3);

    auto result = extractRFC1ReadAlignments(buffer);
    ASSERT_EQ(2u, result.size());

    // Results are ordered by fragment ID (map ordering)
    EXPECT_EQ("TAATCCG", result[0].read);
    EXPECT_EQ("CCGCCGCCTT", result[1].read);
}

TEST(RFC1MotifAnalysisTests, ExtractRFC1ReadAlignments_IsReversedSetCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("TAAT(CCG)*CCTT"));

    locus::AlignmentBuffer buffer;

    // Fragment with forward strand read
    locus::AlignedFragment fragFwd{
        "fragA",
        "TAATCCG",
        decodeGraphAlignment(0, "0[4M]1[3M]", &graph),
        true,
        "CCTT",
        decodeGraphAlignment(0, "2[4M]", &graph),
        false
    };
    buffer.tryInsertFragment(fragFwd);

    // Fragment with reverse strand read
    locus::AlignedFragment fragRev{
        "fragB",
        "CCGCCTT",
        decodeGraphAlignment(0, "1[3M]2[4M]", &graph),
        false,
        "TAAT",
        decodeGraphAlignment(0, "0[4M]", &graph),
        true
    };
    buffer.tryInsertFragment(fragRev);

    auto result = extractRFC1ReadAlignments(buffer);
    ASSERT_EQ(2u, result.size());

    // fragA: readIsForwardStrand=true -> isReversed=false
    EXPECT_EQ("TAATCCG", result[0].read);
    EXPECT_FALSE(result[0].isReversed);

    // fragB: readIsForwardStrand=false -> isReversed=true
    EXPECT_EQ("CCGCCTT", result[1].read);
    EXPECT_TRUE(result[1].isReversed);
}
