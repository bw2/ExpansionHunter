//
// Expansion Hunter
// Copyright 2016-2021 Illumina, Inc.
// All rights reserved.
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

#include "locus/AlignmentBuffer.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/Graph.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;

using namespace ehunter;
using namespace ehunter::locus;

namespace
{

// Helper function to create a simple graph for testing
Graph makeTestGraph()
{
    return makeRegionGraph(decodeFeaturesFromRegex("ATAT(CCG)*ATTT"));
}

// Helper function to create an AlignedFragment for testing
AlignedFragment makeTestFragment(
    const std::string& fragmentId, const Graph& graph, const std::string& readCigar, const std::string& mateCigar,
    const std::string& readBases = "ATATCCGCCGAT", const std::string& mateBases = "CCGCCGATTT",
    bool readIsForward = true, bool mateIsForward = false)
{
    return AlignedFragment{
        fragmentId,
        readBases,
        decodeGraphAlignment(0, readCigar, &graph),
        readIsForward,
        mateBases,
        decodeGraphAlignment(0, mateCigar, &graph),
        mateIsForward
    };
}

} // namespace

// Test AlignedFragment struct construction (all fields set correctly)
TEST(AlignedFragmentTest, ConstructionSetsAllFields)
{
    Graph graph = makeTestGraph();

    AlignedFragment fragment{
        "test_fragment_1",
        "ATATCCG",
        decodeGraphAlignment(0, "0[4M]1[3M]", &graph),
        true,
        "CCGATTT",
        decodeGraphAlignment(0, "1[3M]2[4M]", &graph),
        false
    };

    EXPECT_EQ("test_fragment_1", fragment.fragmentId);
    EXPECT_EQ("ATATCCG", fragment.readBases);
    EXPECT_TRUE(fragment.readIsForwardStrand);
    EXPECT_EQ("CCGATTT", fragment.mateBases);
    EXPECT_FALSE(fragment.mateIsForwardStrand);
}

// Test tryInsertFragment() - successful insertion
TEST(AlignmentBufferTest, TryInsertFragment_SuccessfulInsertion)
{
    Graph graph = makeTestGraph();
    AlignmentBuffer buffer;

    AlignedFragment fragment = makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]");

    bool inserted = buffer.tryInsertFragment(std::move(fragment));

    EXPECT_TRUE(inserted);
    EXPECT_EQ(1u, buffer.size());
    EXPECT_FALSE(buffer.empty());
}

// Test tryInsertFragment() - duplicate fragment ID (map emplace does not replace)
TEST(AlignmentBufferTest, TryInsertFragment_DuplicateFragmentId)
{
    Graph graph = makeTestGraph();
    AlignmentBuffer buffer;

    // Insert first fragment
    AlignedFragment fragment1 = makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]", "FIRST_READ");
    buffer.tryInsertFragment(std::move(fragment1));

    // Insert second fragment with same ID
    AlignedFragment fragment2 = makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]", "SECOND_READ");
    buffer.tryInsertFragment(std::move(fragment2));

    // Map should still have 1 entry (emplace does not replace existing keys)
    EXPECT_EQ(1u, buffer.size());

    // The first fragment should still be there (emplace doesn't overwrite)
    const auto& storedFragment = buffer.getBuffer().at("frag1");
    EXPECT_EQ("FIRST_READ", storedFragment.readBases);
}

// Test getBuffer() - empty buffer returns empty map
TEST(AlignmentBufferTest, GetBuffer_EmptyBufferReturnsEmptyMap)
{
    AlignmentBuffer buffer;

    const auto& bufferMap = buffer.getBuffer();

    EXPECT_TRUE(bufferMap.empty());
    EXPECT_EQ(0u, bufferMap.size());
}

// Test getBuffer() - multiple fragments all retrievable
TEST(AlignmentBufferTest, GetBuffer_MultipleFragmentsAllRetrievable)
{
    Graph graph = makeTestGraph();
    AlignmentBuffer buffer;

    // Insert multiple fragments
    AlignedFragment fragment1 = makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]");
    AlignedFragment fragment2 = makeTestFragment("frag2", graph, "0[4M]1[3M]", "1[3M]2[4M]");
    AlignedFragment fragment3 = makeTestFragment("frag3", graph, "0[4M]1[3M]", "1[3M]2[4M]");

    buffer.tryInsertFragment(std::move(fragment1));
    buffer.tryInsertFragment(std::move(fragment2));
    buffer.tryInsertFragment(std::move(fragment3));

    const auto& bufferMap = buffer.getBuffer();

    EXPECT_EQ(3u, bufferMap.size());
    EXPECT_TRUE(bufferMap.find("frag1") != bufferMap.end());
    EXPECT_TRUE(bufferMap.find("frag2") != bufferMap.end());
    EXPECT_TRUE(bufferMap.find("frag3") != bufferMap.end());

    // Verify fragment IDs are correct
    EXPECT_EQ("frag1", bufferMap.at("frag1").fragmentId);
    EXPECT_EQ("frag2", bufferMap.at("frag2").fragmentId);
    EXPECT_EQ("frag3", bufferMap.at("frag3").fragmentId);
}

// Test size() - returns correct count
TEST(AlignmentBufferTest, Size_ReturnsCorrectCount)
{
    Graph graph = makeTestGraph();
    AlignmentBuffer buffer;

    EXPECT_EQ(0u, buffer.size());

    buffer.tryInsertFragment(makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]"));
    EXPECT_EQ(1u, buffer.size());

    buffer.tryInsertFragment(makeTestFragment("frag2", graph, "0[4M]1[3M]", "1[3M]2[4M]"));
    EXPECT_EQ(2u, buffer.size());

    buffer.tryInsertFragment(makeTestFragment("frag3", graph, "0[4M]1[3M]", "1[3M]2[4M]"));
    EXPECT_EQ(3u, buffer.size());
}

// Test empty() - true when empty, false after insert
TEST(AlignmentBufferTest, Empty_TrueWhenEmpty_FalseAfterInsert)
{
    Graph graph = makeTestGraph();
    AlignmentBuffer buffer;

    EXPECT_TRUE(buffer.empty());

    buffer.tryInsertFragment(makeTestFragment("frag1", graph, "0[4M]1[3M]", "1[3M]2[4M]"));

    EXPECT_FALSE(buffer.empty());
}
