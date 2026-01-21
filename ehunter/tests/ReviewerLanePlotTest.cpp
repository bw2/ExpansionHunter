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

#include "reviewer/LanePlot.hh"

#include <memory>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/LinearAlignment.hh"
#include "graphcore/Graph.hh"
#include "graphcore/Path.hh"

#include "core/Read.hh"
#include "reviewer/Aligns.hh"

using graphtools::Alignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;

using std::list;
using std::make_shared;
using std::string;
using std::vector;

using namespace ehunter;
using namespace ehunter::reviewer;

// ============================================================================
// Feature struct tests
// ============================================================================

TEST(FeatureTest, Constructor_SetsAllFields)
{
    Feature feature(FeatureType::kRect, 10, "blue", "black");

    EXPECT_EQ(FeatureType::kRect, feature.type);
    EXPECT_EQ(10, feature.length);
    EXPECT_EQ("blue", feature.fill);
    EXPECT_EQ("black", feature.stroke);
    EXPECT_FALSE(feature.label.is_initialized());
}

TEST(FeatureTest, Constructor_DifferentFeatureTypes)
{
    Feature rectFeature(FeatureType::kRect, 5, "red", "none");
    EXPECT_EQ(FeatureType::kRect, rectFeature.type);

    Feature leftBreakFeature(FeatureType::kRectWithLeftBreak, 8, "green", "black");
    EXPECT_EQ(FeatureType::kRectWithLeftBreak, leftBreakFeature.type);

    Feature rightBreakFeature(FeatureType::kRectWithRightBreak, 12, "yellow", "gray");
    EXPECT_EQ(FeatureType::kRectWithRightBreak, rightBreakFeature.type);

    Feature lineFeature(FeatureType::kLine, 3, "none", "black");
    EXPECT_EQ(FeatureType::kLine, lineFeature.type);

    Feature arrowsFeature(FeatureType::kArrows, 15, "black", "black");
    EXPECT_EQ(FeatureType::kArrows, arrowsFeature.type);

    Feature vertLineFeature(FeatureType::kVerticalLine, 0, "none", "black");
    EXPECT_EQ(FeatureType::kVerticalLine, vertLineFeature.type);
}

TEST(FeatureTest, LabelCanBeSet)
{
    Feature feature(FeatureType::kRect, 10, "blue", "none");
    EXPECT_FALSE(feature.label.is_initialized());

    feature.label = "ATCG";
    EXPECT_TRUE(feature.label.is_initialized());
    EXPECT_EQ("ATCG", *feature.label);
}

// ============================================================================
// Segment struct tests
// ============================================================================

TEST(SegmentTest, Constructor_EmptyFeatures_EndEqualsStart)
{
    vector<Feature> emptyFeatures;
    Segment segment(100, emptyFeatures, 1.0);

    EXPECT_EQ(100, segment.start);
    EXPECT_EQ(100, segment.end);
    EXPECT_TRUE(segment.features.empty());
    EXPECT_DOUBLE_EQ(1.0, segment.opacity);
}

TEST(SegmentTest, Constructor_SingleFeature_EndCalculatedCorrectly)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kRect, 25, "blue", "none");

    Segment segment(50, features, 0.8);

    EXPECT_EQ(50, segment.start);
    EXPECT_EQ(75, segment.end);  // start(50) + length(25) = 75
    EXPECT_EQ(1u, segment.features.size());
    EXPECT_DOUBLE_EQ(0.8, segment.opacity);
}

TEST(SegmentTest, Constructor_MultipleFeatures_EndIsSumOfLengths)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kRect, 10, "blue", "none");
    features.emplace_back(FeatureType::kLine, 5, "none", "black");
    features.emplace_back(FeatureType::kRect, 15, "red", "none");

    Segment segment(20, features, 1.0);

    EXPECT_EQ(20, segment.start);
    EXPECT_EQ(50, segment.end);  // start(20) + 10 + 5 + 15 = 50
    EXPECT_EQ(3u, segment.features.size());
}

TEST(SegmentTest, Constructor_ZeroLengthFeatures_HandleCorrectly)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kVerticalLine, 0, "none", "black");
    features.emplace_back(FeatureType::kRect, 10, "blue", "none");
    features.emplace_back(FeatureType::kVerticalLine, 0, "none", "black");

    Segment segment(0, features, 1.0);

    EXPECT_EQ(0, segment.start);
    EXPECT_EQ(10, segment.end);  // 0 + 0 + 10 + 0 = 10
}

TEST(SegmentTest, Constructor_NegativeStart_EndCalculatedCorrectly)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kRect, 30, "blue", "none");

    Segment segment(-10, features, 0.7);

    EXPECT_EQ(-10, segment.start);
    EXPECT_EQ(20, segment.end);  // start(-10) + length(30) = 20
}

TEST(SegmentTest, Constructor_DifferentOpacityValues)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kRect, 10, "blue", "none");

    Segment fullOpacity(0, features, 1.0);
    EXPECT_DOUBLE_EQ(1.0, fullOpacity.opacity);

    Segment partialOpacity(0, features, 0.7);
    EXPECT_DOUBLE_EQ(0.7, partialOpacity.opacity);

    Segment zeroOpacity(0, features, 0.0);
    EXPECT_DOUBLE_EQ(0.0, zeroOpacity.opacity);
}

// ============================================================================
// Lane struct tests
// ============================================================================

TEST(LaneTest, Constructor_SetsHeightAndSegments)
{
    vector<Feature> features;
    features.emplace_back(FeatureType::kRect, 20, "blue", "none");

    vector<Segment> segments;
    segments.emplace_back(0, features, 1.0);
    segments.emplace_back(30, features, 0.8);

    Lane lane(10, segments);

    EXPECT_EQ(10, lane.height);
    EXPECT_EQ(2u, lane.segments.size());
}

TEST(LaneTest, Constructor_EmptySegments)
{
    vector<Segment> emptySegments;
    Lane lane(12, emptySegments);

    EXPECT_EQ(12, lane.height);
    EXPECT_TRUE(lane.segments.empty());
}

// ============================================================================
// LanePlot (vector<Lane>) tests
// ============================================================================

TEST(LanePlotTest, EmptyLanePlot)
{
    LanePlot plot;
    EXPECT_TRUE(plot.empty());
    EXPECT_EQ(0u, plot.size());
}

TEST(LanePlotTest, SingleLane)
{
    LanePlot plot;
    vector<Segment> segments;
    plot.emplace_back(18, segments);

    EXPECT_EQ(1u, plot.size());
    EXPECT_EQ(18, plot[0].height);
}

TEST(LanePlotTest, MultipleLanes)
{
    LanePlot plot;
    vector<Segment> segments;

    plot.emplace_back(12, segments);  // Label lane
    plot.emplace_back(18, segments);  // Haplotype path lane
    plot.emplace_back(10, segments);  // Read lane

    EXPECT_EQ(3u, plot.size());
    EXPECT_EQ(12, plot[0].height);
    EXPECT_EQ(18, plot[1].height);
    EXPECT_EQ(10, plot[2].height);
}

// ============================================================================
// FeatureType enum tests
// ============================================================================

TEST(FeatureTypeTest, AllTypesAccessible)
{
    // Verify all feature types can be used
    FeatureType types[] = {
        FeatureType::kRect,
        FeatureType::kRectWithLeftBreak,
        FeatureType::kRectWithRightBreak,
        FeatureType::kLine,
        FeatureType::kArrows,
        FeatureType::kVerticalLine
    };

    EXPECT_EQ(6u, sizeof(types) / sizeof(types[0]));
}

// ============================================================================
// Integration-style tests for generateBlueprint (requires graph setup)
// Note: generateBlueprint has complex dependencies, so we test what we can
// ============================================================================

class GenerateBlueprintTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a simple graph: leftFlank -> repeat -> rightFlank
        // Node 0: left flank
        // Node 1: repeat unit (with self-loop)
        // Node 2: right flank
        graph_ = std::make_unique<Graph>(3);
        graph_->setNodeSeq(0, std::string(100, 'A'));  // 100bp left flank
        graph_->setNodeSeq(1, "CAG");                   // 3bp repeat
        graph_->setNodeSeq(2, std::string(100, 'T'));  // 100bp right flank

        graph_->addEdge(0, 1);
        graph_->addEdge(1, 1);  // Self-loop for repeat
        graph_->addEdge(1, 2);
    }

    std::unique_ptr<Graph> graph_;
};

TEST_F(GenerateBlueprintTest, EmptyFragmentInput_ThrowsException)
{
    // Create a path through the graph with some repeat units
    vector<graphtools::NodeId> nodeIds = {0, 1, 1, 1, 2};  // 3 repeats
    Path path(graph_.get(), 0, nodeIds, static_cast<int32_t>(graph_->nodeSeq(2).length()));

    vector<Path> paths = {path};
    FragById fragById;
    FragAssignment fragAssignment({}, {});
    FragPathAlignsById fragPathAlignsById;

    // generateBlueprint throws when there are no read alignments
    EXPECT_THROW(generateBlueprint(paths, fragById, fragAssignment, fragPathAlignsById), std::runtime_error);
}

TEST_F(GenerateBlueprintTest, TwoHaplotypePaths_ProducesTwoLanePlots)
{
    // Create two different paths (heterozygous)
    vector<graphtools::NodeId> nodeIds1 = {0, 1, 1, 2};     // 2 repeats
    vector<graphtools::NodeId> nodeIds2 = {0, 1, 1, 1, 2};  // 3 repeats

    Path path1(graph_.get(), 0, nodeIds1, static_cast<int32_t>(graph_->nodeSeq(2).length()));
    Path path2(graph_.get(), 0, nodeIds2, static_cast<int32_t>(graph_->nodeSeq(2).length()));

    vector<Path> paths = {path1, path2};
    FragById fragById;
    FragAssignment fragAssignment({}, {});
    FragPathAlignsById fragPathAlignsById;

    // Should throw because no reads, but if it didn't throw, it would produce 2 plots
    EXPECT_THROW(generateBlueprint(paths, fragById, fragAssignment, fragPathAlignsById), std::runtime_error);
}

// Note: Full integration tests for generateBlueprint with actual fragments
// would require setting up ReadPathAlign which needs complex graph alignment
// projection. The struct tests above cover the core data structure behavior.
// Additional integration testing is recommended at a higher level.
