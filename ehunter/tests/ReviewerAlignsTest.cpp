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

#include "reviewer/Aligns.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using namespace ehunter;
using namespace ehunter::reviewer;

// Test ReadPathAlign constructor sets all fields correctly
TEST(ReadPathAlignConstruction, AllFieldsSetCorrectly)
{
    // Create a simple graph: "ATTCGA(C)*ATGTCG"
    // Node 0: ATTCGA (6 bases)
    // Node 1: C (1 base, repeat)
    // Node 2: ATGTCG (6 bases)
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create a haplotype path spanning nodes 0, 1, 1, 1, 2
    Path hapPath(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    // Create a graph alignment starting at position 3 on node 0
    GraphAlignment align = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph);
    auto alignPtr = std::make_shared<GraphAlign>(align);

    int pathIndex = 0;
    int startIndexOnPath = 0;  // alignment starts at node 0

    ReadPathAlign readPathAlign(hapPath, pathIndex, startIndexOnPath, alignPtr);

    EXPECT_EQ(0, readPathAlign.pathIndex);
    EXPECT_EQ(0, readPathAlign.startIndexOnPath);
    EXPECT_NE(nullptr, readPathAlign.align);
}

// Test ReadPathAlign begin/end calculation with alignment starting at first node
TEST(ReadPathAlignBeginEnd, AlignmentAtFirstNode_CorrectPositions)
{
    // Graph: ATTCGA(C)*ATGTCG
    // Node 0: ATTCGA (6 bases)
    // Node 1: C (1 base)
    // Node 2: ATGTCG (6 bases)
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Haplotype path: 0 -> 1 -> 1 -> 1 -> 2
    Path hapPath(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    // Alignment starting at position 3 on node 0, spanning 3M on node 0 and 1M on node 1, and 1M on node 1
    // Total reference length: 3 + 1 + 1 = 5
    GraphAlignment align = decodeGraphAlignment(3, "0[3M]1[1M]1[1M]", &graph);
    auto alignPtr = std::make_shared<GraphAlign>(align);

    ReadPathAlign readPathAlign(hapPath, 0, 0, alignPtr);

    // startIndexOnPath = 0, so prefix length = 0
    // begin = 0 + align.path().startPosition() = 0 + 3 = 3
    // end = begin + referenceLength = 3 + 5 = 8
    EXPECT_EQ(3, readPathAlign.begin);
    EXPECT_EQ(8, readPathAlign.end);
}

// Test ReadPathAlign begin/end calculation with alignment starting at later node
TEST(ReadPathAlignBeginEnd, AlignmentAtLaterNode_CorrectPositions)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Haplotype path: 0 -> 1 -> 1 -> 1 -> 2
    Path hapPath(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    // Alignment starting at position 0 on node 1 (first repeat unit)
    GraphAlignment align = decodeGraphAlignment(0, "1[1M]1[1M]2[4M]", &graph);
    auto alignPtr = std::make_shared<GraphAlign>(align);

    // startIndexOnPath = 1, meaning alignment path starts at node index 1 of hapPath
    // Prefix length = length of node 0 = 6
    ReadPathAlign readPathAlign(hapPath, 0, 1, alignPtr);

    // begin = 6 (prefix) + 0 (startPosition) = 6
    // referenceLength = 1 + 1 + 4 = 6
    // end = 6 + 6 = 12
    EXPECT_EQ(6, readPathAlign.begin);
    EXPECT_EQ(12, readPathAlign.end);
}

// Test ReadPathAlign with alignment starting at right flank
TEST(ReadPathAlignBeginEnd, AlignmentAtRightFlank_CorrectPositions)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Haplotype path: 0 -> 1 -> 1 -> 2
    Path hapPath(&graph, 0, { 0, 1, 1, 2 }, 6);

    // Alignment starting at position 1 on node 2 (right flank), spanning 5 bases
    GraphAlignment align = decodeGraphAlignment(1, "2[5M]", &graph);
    auto alignPtr = std::make_shared<GraphAlign>(align);

    // startIndexOnPath = 3 (node index of node 2 in the path)
    // Prefix length = 6 (node 0) + 1 (node 1) + 1 (node 1) = 8
    ReadPathAlign readPathAlign(hapPath, 0, 3, alignPtr);

    // begin = 8 (prefix) + 1 (startPosition) = 9
    // end = 9 + 5 = 14
    EXPECT_EQ(9, readPathAlign.begin);
    EXPECT_EQ(14, readPathAlign.end);
}

// Test FragAssignment constructor sets fields correctly
TEST(FragAssignmentConstruction, FieldsSetCorrectly)
{
    std::vector<std::string> fragIds = { "frag1", "frag2", "frag3" };
    std::vector<int> alignIndexByFrag = { 0, 1, 0 };

    FragAssignment assignment(fragIds, alignIndexByFrag);

    EXPECT_EQ(3u, assignment.fragIds.size());
    EXPECT_EQ("frag1", assignment.fragIds[0]);
    EXPECT_EQ("frag2", assignment.fragIds[1]);
    EXPECT_EQ("frag3", assignment.fragIds[2]);

    EXPECT_EQ(3u, assignment.alignIndexByFrag.size());
    EXPECT_EQ(0, assignment.alignIndexByFrag[0]);
    EXPECT_EQ(1, assignment.alignIndexByFrag[1]);
    EXPECT_EQ(0, assignment.alignIndexByFrag[2]);
}

// Test FragAssignment with empty vectors
TEST(FragAssignmentConstruction, EmptyVectors)
{
    std::vector<std::string> fragIds;
    std::vector<int> alignIndexByFrag;

    FragAssignment assignment(fragIds, alignIndexByFrag);

    EXPECT_TRUE(assignment.fragIds.empty());
    EXPECT_TRUE(assignment.alignIndexByFrag.empty());
}

// Test FragPathAlign construction
TEST(FragPathAlignConstruction, BothReadsOnSamePath)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Haplotype path: 0 -> 1 -> 1 -> 2
    Path hapPath(&graph, 0, { 0, 1, 1, 2 }, 6);

    // Read alignment
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[6M]1[1M]", &graph);
    auto readAlignPtr = std::make_shared<GraphAlign>(readAlign);
    ReadPathAlign readPathAlign(hapPath, 0, 0, readAlignPtr);

    // Mate alignment (on same path index)
    GraphAlignment mateAlign = decodeGraphAlignment(0, "1[1M]2[4M]", &graph);
    auto mateAlignPtr = std::make_shared<GraphAlign>(mateAlign);
    ReadPathAlign matePathAlign(hapPath, 0, 1, mateAlignPtr);

    FragPathAlign fragPathAlign(readPathAlign, matePathAlign);

    EXPECT_EQ(0, fragPathAlign.readAlign.pathIndex);
    EXPECT_EQ(0, fragPathAlign.mateAlign.pathIndex);
}

// Test ReadWithAlign struct
TEST(ReadWithAlignConstruction, FieldsSetCorrectly)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    GraphAlignment align = decodeGraphAlignment(0, "0[6M]", &graph);

    ReadId readId("frag1", MateNumber::kFirstMate);
    ehunter::Read read(readId, "ATTCGA", false);

    ReadWithAlign readWithAlign(read, align);

    EXPECT_EQ("ATTCGA", readWithAlign.bases());
    EXPECT_EQ("frag1", readWithAlign.read.fragmentId());
}

// Test Frag struct construction
TEST(FragConstruction, BothReadsStored)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));
    GraphAlignment readAlign = decodeGraphAlignment(0, "0[6M]", &graph);
    GraphAlignment mateAlign = decodeGraphAlignment(0, "2[6M]", &graph);

    ReadId readId("frag1", MateNumber::kFirstMate);
    ehunter::Read read(readId, "ATTCGA", false);
    ReadWithAlign readWithAlignObj(read, readAlign);

    ReadId mateId("frag1", MateNumber::kSecondMate);
    ehunter::Read mate(mateId, "ATGTCG", true);
    ReadWithAlign mateWithAlignObj(mate, mateAlign);

    Frag frag(readWithAlignObj, mateWithAlignObj);

    EXPECT_EQ("ATTCGA", frag.read.bases());
    EXPECT_EQ("ATGTCG", frag.mate.bases());
    EXPECT_EQ("frag1", frag.read.read.fragmentId());
    EXPECT_EQ("frag1", frag.mate.read.fragmentId());
}
