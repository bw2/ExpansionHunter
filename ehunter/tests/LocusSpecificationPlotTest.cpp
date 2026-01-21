//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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

#include "locus/LocusSpecification.hh"

#include "gtest/gtest.h"

#include "graphcore/Graph.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using namespace ehunter;

namespace
{

// Helper function to create a minimal LocusSpecification for testing
LocusSpecification createTestLocusSpec(
    bool useRFC1MotifAnalysis = false,
    std::vector<PlotReadVisualization> plotConditions = {})
{
    graphtools::Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(CAG)*ATGTCG"));
    std::vector<GenomicRegion> targetRegions = { GenomicRegion(1, 100, 200) };
    NodeToRegionAssociation referenceRegions;
    GenotyperParameters genotyperParams(10);

    return LocusSpecification(
        "TestLocus",
        ChromType::kAutosome,
        targetRegions,
        graph,
        referenceRegions,
        genotyperParams,
        useRFC1MotifAnalysis,
        plotConditions);
}

} // namespace

// Test: Default PlotPolicy is kConditional
TEST(LocusSpecificationPlotTest, DefaultPlotPolicyIsConditional)
{
    LocusSpecification locusSpec = createTestLocusSpec();
    EXPECT_EQ(PlotPolicy::kConditional, locusSpec.plotPolicy());
}

// Test: setPlotPolicy(kAll) and verify getter
TEST(LocusSpecificationPlotTest, SetPlotPolicyAll)
{
    LocusSpecification locusSpec = createTestLocusSpec();
    locusSpec.setPlotPolicy(PlotPolicy::kAll);
    EXPECT_EQ(PlotPolicy::kAll, locusSpec.plotPolicy());
}

// Test: setPlotPolicy(kNone) and verify getter
TEST(LocusSpecificationPlotTest, SetPlotPolicyNone)
{
    LocusSpecification locusSpec = createTestLocusSpec();
    locusSpec.setPlotPolicy(PlotPolicy::kNone);
    EXPECT_EQ(PlotPolicy::kNone, locusSpec.plotPolicy());
}

// Test: hasPlotConditions() false when empty
TEST(LocusSpecificationPlotTest, HasPlotConditionsFalseWhenEmpty)
{
    LocusSpecification locusSpec = createTestLocusSpec(false, {});
    EXPECT_FALSE(locusSpec.hasPlotConditions());
}

// Test: hasPlotConditions() true when conditions present
TEST(LocusSpecificationPlotTest, HasPlotConditionsTrueWhenPresent)
{
    std::vector<PlotReadVisualization> conditions = {
        { PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 50 }
    };
    LocusSpecification locusSpec = createTestLocusSpec(false, conditions);
    EXPECT_TRUE(locusSpec.hasPlotConditions());
}

// Test: plotConditions() returns correct conditions
TEST(LocusSpecificationPlotTest, PlotConditionsReturnsCorrectConditions)
{
    std::vector<PlotReadVisualization> conditions = {
        { PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 50 },
        { PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kLessThan, 10 }
    };
    LocusSpecification locusSpec = createTestLocusSpec(false, conditions);

    const auto& returnedConditions = locusSpec.plotConditions();
    ASSERT_EQ(2u, returnedConditions.size());

    EXPECT_EQ(PlotThresholdAppliedTo::kLongAllele, returnedConditions[0].appliedTo);
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThanOrEqual, returnedConditions[0].op);
    EXPECT_EQ(50, returnedConditions[0].threshold);

    EXPECT_EQ(PlotThresholdAppliedTo::kShortAllele, returnedConditions[1].appliedTo);
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThan, returnedConditions[1].op);
    EXPECT_EQ(10, returnedConditions[1].threshold);
}

// Test: requiresAlignmentBuffer() true when RFC1=true
TEST(LocusSpecificationPlotTest, RequiresAlignmentBufferTrueWhenRFC1True)
{
    LocusSpecification locusSpec = createTestLocusSpec(true, {});
    // Ensure we're testing the RFC1 path (not other conditions)
    EXPECT_EQ(PlotPolicy::kConditional, locusSpec.plotPolicy());
    EXPECT_FALSE(locusSpec.hasPlotConditions());
    EXPECT_TRUE(locusSpec.requiresAlignmentBuffer());
}

// Test: requiresAlignmentBuffer() true when plotPolicy=kAll
TEST(LocusSpecificationPlotTest, RequiresAlignmentBufferTrueWhenPlotPolicyAll)
{
    LocusSpecification locusSpec = createTestLocusSpec(false, {});
    locusSpec.setPlotPolicy(PlotPolicy::kAll);
    // Ensure we're testing the plotPolicy path (not other conditions)
    EXPECT_FALSE(locusSpec.useRFC1MotifAnalysis());
    EXPECT_FALSE(locusSpec.hasPlotConditions());
    EXPECT_TRUE(locusSpec.requiresAlignmentBuffer());
}

// Test: requiresAlignmentBuffer() true when hasPlotConditions()
TEST(LocusSpecificationPlotTest, RequiresAlignmentBufferTrueWhenHasPlotConditions)
{
    std::vector<PlotReadVisualization> conditions = {
        { PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 50 }
    };
    LocusSpecification locusSpec = createTestLocusSpec(false, conditions);
    // Ensure we're testing the hasPlotConditions path (not other conditions)
    EXPECT_FALSE(locusSpec.useRFC1MotifAnalysis());
    EXPECT_EQ(PlotPolicy::kConditional, locusSpec.plotPolicy());
    EXPECT_TRUE(locusSpec.requiresAlignmentBuffer());
}

// Test: requiresAlignmentBuffer() false when none of the above
TEST(LocusSpecificationPlotTest, RequiresAlignmentBufferFalseWhenNoConditions)
{
    LocusSpecification locusSpec = createTestLocusSpec(false, {});
    // Ensure none of the conditions are true
    EXPECT_FALSE(locusSpec.useRFC1MotifAnalysis());
    EXPECT_EQ(PlotPolicy::kConditional, locusSpec.plotPolicy());
    EXPECT_FALSE(locusSpec.hasPlotConditions());
    EXPECT_FALSE(locusSpec.requiresAlignmentBuffer());
}
