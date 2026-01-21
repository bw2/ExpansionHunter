//
// ExpansionHunter
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

#include "reviewer/Origin.hh"

#include "gtest/gtest.h"

#include <vector>
#include <string>
#include <set>

#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"
#include "graphalign/GraphAlignment.hh"

using namespace ehunter::reviewer;
using graphtools::Graph;
using graphtools::NodeId;
using graphtools::Path;
using graphtools::GraphAlignment;
using std::string;
using std::vector;

namespace
{

// Helper to create a simple graph for testing
Graph makeSimpleGraph()
{
    // Create a simple graph: node0 -> node1 -> node2
    Graph graph(3);
    graph.setNodeSeq(0, "AAAA");
    graph.setNodeSeq(1, "CCCC");
    graph.setNodeSeq(2, "GGGG");
    graph.addEdge(0, 1);
    graph.addEdge(1, 2);
    return graph;
}

// Helper to create a GraphAlignPtr
GraphAlignPtr makeGraphAlign(Graph* graph, int startNode)
{
    // Create a simple alignment: matches starting at startNode
    vector<NodeId> nodes = { static_cast<NodeId>(startNode) };
    Path path(graph, 0, nodes, 4);
    vector<graphtools::Alignment> alignments = { graphtools::Alignment(0, "4M") };
    return std::make_shared<GraphAlign>(path, alignments);
}

// Helper to create a ReadPathAlign
ReadPathAlign makeReadPathAlign(const Path& hapPath, int pathIndex, int startIndexOnPath, GraphAlignPtr align)
{
    return ReadPathAlign(hapPath, pathIndex, startIndexOnPath, align);
}

// Helper to create a FragPathAlign
FragPathAlign makeFragPathAlign(const Path& hapPath, int pathIndex, GraphAlignPtr readAlign, GraphAlignPtr mateAlign)
{
    ReadPathAlign readPathAlign = makeReadPathAlign(hapPath, pathIndex, 0, readAlign);
    ReadPathAlign matePathAlign = makeReadPathAlign(hapPath, pathIndex, 1, mateAlign);
    return FragPathAlign(readPathAlign, matePathAlign);
}

}  // namespace

TEST(ReviewerOrigin, EmptyInput_ReturnsEmptyAssignment)
{
    vector<Path> hapPaths;
    FragPathAlignsById fragPathAlignsById;

    FragAssignment assignment = getBestFragAssignment(hapPaths, fragPathAlignsById);

    EXPECT_TRUE(assignment.fragIds.empty());
    EXPECT_TRUE(assignment.alignIndexByFrag.empty());
}

TEST(ReviewerOrigin, SingleFragment_AssignsToOneHaplotype)
{
    Graph graph = makeSimpleGraph();
    vector<NodeId> nodeIds = { 0, 1, 2 };
    Path hapPath(&graph, 0, nodeIds, 4);
    vector<Path> hapPaths = { hapPath };

    GraphAlignPtr align = makeGraphAlign(&graph, 0);
    FragPathAlign fragAlign = makeFragPathAlign(hapPath, 0, align, align);

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"] = { fragAlign };

    FragAssignment assignment = getBestFragAssignment(hapPaths, fragPathAlignsById);

    ASSERT_EQ(1u, assignment.fragIds.size());
    ASSERT_EQ(1u, assignment.alignIndexByFrag.size());
    EXPECT_EQ("frag1", assignment.fragIds[0]);
    EXPECT_EQ(0, assignment.alignIndexByFrag[0]);
}

TEST(ReviewerOrigin, MultipleFragments_EachGetsAssignment)
{
    Graph graph = makeSimpleGraph();
    vector<NodeId> nodeIds = { 0, 1, 2 };
    Path hapPath(&graph, 0, nodeIds, 4);
    vector<Path> hapPaths = { hapPath };

    GraphAlignPtr align = makeGraphAlign(&graph, 0);
    FragPathAlign fragAlign = makeFragPathAlign(hapPath, 0, align, align);

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"] = { fragAlign };
    fragPathAlignsById["frag2"] = { fragAlign };
    fragPathAlignsById["frag3"] = { fragAlign };

    FragAssignment assignment = getBestFragAssignment(hapPaths, fragPathAlignsById);

    ASSERT_EQ(3u, assignment.fragIds.size());
    ASSERT_EQ(3u, assignment.alignIndexByFrag.size());

    // Verify all fragment IDs are present
    std::set<string> expectedFragIds = { "frag1", "frag2", "frag3" };
    std::set<string> actualFragIds(assignment.fragIds.begin(), assignment.fragIds.end());
    EXPECT_EQ(expectedFragIds, actualFragIds);

    // Each fragment should have a valid assignment index
    for (size_t i = 0; i < assignment.alignIndexByFrag.size(); ++i)
    {
        EXPECT_EQ(0, assignment.alignIndexByFrag[i]);  // Only one alignment option
    }
}

TEST(ReviewerOrigin, AmbiguousFragment_PicksOneHaplotype)
{
    Graph graph = makeSimpleGraph();
    vector<NodeId> nodeIds = { 0, 1, 2 };
    Path hapPath0(&graph, 0, nodeIds, 4);
    Path hapPath1(&graph, 0, nodeIds, 4);
    vector<Path> hapPaths = { hapPath0, hapPath1 };

    GraphAlignPtr align = makeGraphAlign(&graph, 0);

    // Create two different alignments for the same fragment (ambiguous)
    FragPathAlign fragAlign0 = makeFragPathAlign(hapPath0, 0, align, align);
    FragPathAlign fragAlign1 = makeFragPathAlign(hapPath1, 1, align, align);

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["ambiguous_frag"] = { fragAlign0, fragAlign1 };

    FragAssignment assignment = getBestFragAssignment(hapPaths, fragPathAlignsById);

    ASSERT_EQ(1u, assignment.fragIds.size());
    ASSERT_EQ(1u, assignment.alignIndexByFrag.size());
    EXPECT_EQ("ambiguous_frag", assignment.fragIds[0]);

    // The assignment index should be either 0 or 1 (one of the two options)
    EXPECT_TRUE(assignment.alignIndexByFrag[0] == 0 || assignment.alignIndexByFrag[0] == 1);
}

TEST(ReviewerOrigin, MultipleFragmentsWithMultipleAlignments_AllGetValidAssignments)
{
    Graph graph = makeSimpleGraph();
    vector<NodeId> nodeIds = { 0, 1, 2 };
    Path hapPath0(&graph, 0, nodeIds, 4);
    Path hapPath1(&graph, 0, nodeIds, 4);
    vector<Path> hapPaths = { hapPath0, hapPath1 };

    GraphAlignPtr align = makeGraphAlign(&graph, 0);

    // Frag1 has 2 alignment options
    FragPathAlign fragAlign0 = makeFragPathAlign(hapPath0, 0, align, align);
    FragPathAlign fragAlign1 = makeFragPathAlign(hapPath1, 1, align, align);

    // Frag2 has 3 alignment options
    FragPathAlign fragAlign2 = makeFragPathAlign(hapPath0, 0, align, align);

    FragPathAlignsById fragPathAlignsById;
    fragPathAlignsById["frag1"] = { fragAlign0, fragAlign1 };
    fragPathAlignsById["frag2"] = { fragAlign2, fragAlign0, fragAlign1 };

    FragAssignment assignment = getBestFragAssignment(hapPaths, fragPathAlignsById);

    ASSERT_EQ(2u, assignment.fragIds.size());
    ASSERT_EQ(2u, assignment.alignIndexByFrag.size());

    // Each fragment should have a valid assignment
    for (size_t i = 0; i < assignment.fragIds.size(); ++i)
    {
        const string& fragId = assignment.fragIds[i];
        int numAligns = static_cast<int>(fragPathAlignsById.at(fragId).size());
        EXPECT_GE(assignment.alignIndexByFrag[i], 0);
        EXPECT_LT(assignment.alignIndexByFrag[i], numAligns);
    }
}
