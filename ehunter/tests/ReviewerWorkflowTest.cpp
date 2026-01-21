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

#include "reviewer/ReviewerWorkflow.hh"

#include "gtest/gtest.h"

#include "core/CountTable.hh"
#include "core/Parameters.hh"
#include "genotyping/RepeatGenotype.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantFindings.hh"

using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Helper to create a minimal LocusSpecification for testing
LocusSpecification createTestLocusSpec(
    PlotPolicy policy,
    std::vector<PlotReadVisualization> conditions = {})
{
    // Create minimal graph with single node
    graphtools::Graph graph(1);
    graph.setNodeSeq(0, "ACGT");

    // Create minimal parameters (minLocusCoverage = 10)
    GenotyperParameters params(10);

    // Create LocusSpecification
    LocusSpecification spec(
        "TEST_LOCUS",
        ChromType::kAutosome,
        {GenomicRegion(0, 100, 200)},
        graph,
        {},  // referenceRegions
        params,
        false,  // useRFC1MotifAnalysis
        conditions);

    spec.setPlotPolicy(policy);
    return spec;
}

// Helper to create LocusFindings with a given genotype
LocusFindings createFindingsWithGenotype(int shortAllele, int longAllele)
{
    LocusFindings findings;

    // Create empty CountTables
    CountTable spanningCounts;
    CountTable flankingCounts;
    CountTable inrepeatCounts;

    // Create genotype with the specified allele sizes (3-bp repeat unit)
    RepeatGenotype genotype(3, {shortAllele, longAllele});

    // Create RepeatFindings with the genotype
    auto repeatFindings = std::make_unique<RepeatFindings>(
        spanningCounts,
        flankingCounts,
        inrepeatCounts,
        AlleleCount::kTwo,
        genotype,
        static_cast<GenotypeFilter>(0));

    findings.findingsForEachVariant["TEST_VARIANT"] = std::move(repeatFindings);
    return findings;
}

// Helper to create LocusFindings without a genotype
LocusFindings createFindingsWithoutGenotype()
{
    LocusFindings findings;
    return findings;
}

}  // namespace

//
// Tests for shouldPlotReadVisualization() with different plot policies
//

TEST(ShouldPlotReadVisualization, PolicyNone_AlwaysReturnsFalse)
{
    auto spec = createTestLocusSpec(PlotPolicy::kNone);
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

TEST(ShouldPlotReadVisualization, PolicyAll_AlwaysReturnsTrue)
{
    auto spec = createTestLocusSpec(PlotPolicy::kAll);
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ShouldPlotReadVisualization, PolicyConditionalNoConditions_ReturnsFalse)
{
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, {});
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

TEST(ShouldPlotReadVisualization, PolicyConditionalConditionMet_ReturnsTrue)
{
    // Condition: longAllele >= 15
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 15}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ShouldPlotReadVisualization, PolicyConditionalConditionNotMet_ReturnsFalse)
{
    // Condition: longAllele >= 50 (not met since longAllele is 20)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 50}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

//
// Tests for condition evaluation via shouldPlotReadVisualization (indirect testing of evaluateCondition)
//

TEST(ConditionEvaluation, LessThan_True)
{
    // Condition: shortAllele < 20 (shortAllele is 10)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kLessThan, 20}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, LessThan_False)
{
    // Condition: shortAllele < 20 (shortAllele is 30)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kLessThan, 20}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(30, 40);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, LessThanOrEqual_True)
{
    // Condition: longAllele <= 20 (longAllele is 20)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kLessThanOrEqual, 20}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 20);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, Equal_True)
{
    // Condition: shortAllele == 20 (shortAllele is 20)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kEqual, 20}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(20, 30);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, NotEqual_True)
{
    // Condition: shortAllele != 20 (shortAllele is 10)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kNotEqual, 20}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, GreaterThanOrEqual_True)
{
    // Condition: longAllele >= 400 (longAllele is 400)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 400}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(100, 400);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, GreaterThan_True)
{
    // Condition: longAllele > 400 (longAllele is 500)
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThan, 400}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(100, 500);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

//
// Tests for appliedTo (kShortAllele vs kLongAllele)
//

TEST(ConditionEvaluation, AppliedToShortAllele_UsesShortAllele)
{
    // Condition: shortAllele >= 15, shortAllele is 10, longAllele is 30
    // This should be false because shortAllele (10) < 15
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 15}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, AppliedToLongAllele_UsesLongAllele)
{
    // Condition: longAllele >= 15, shortAllele is 10, longAllele is 30
    // This should be true because longAllele (30) >= 15
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThanOrEqual, 15}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

//
// Tests for multiple conditions (OR logic)
//

TEST(ConditionEvaluation, MultipleConditionsOneMet_ReturnsTrue)
{
    // Conditions: shortAllele < 5 OR longAllele > 25
    // shortAllele is 10 (not < 5), longAllele is 30 (> 25) -> second condition met
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kLessThan, 5},
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThan, 25}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}

TEST(ConditionEvaluation, MultipleConditionsNoneMet_ReturnsFalse)
{
    // Conditions: shortAllele < 5 OR longAllele > 50
    // shortAllele is 10 (not < 5), longAllele is 30 (not > 50) -> neither condition met
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kShortAllele, PlotThresholdComparisonOp::kLessThan, 5},
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kGreaterThan, 50}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithGenotype(10, 30);

    EXPECT_FALSE(shouldPlotReadVisualization(spec, findings));
}

//
// Edge cases
//

TEST(ShouldPlotReadVisualization, NoGenotype_ConditionsUseZeroForAlleles)
{
    // When there's no genotype, evaluateCondition uses default allele values (0)
    // Condition: longAllele < 10 should be true since 0 < 10
    std::vector<PlotReadVisualization> conditions = {
        {PlotThresholdAppliedTo::kLongAllele, PlotThresholdComparisonOp::kLessThan, 10}
    };
    auto spec = createTestLocusSpec(PlotPolicy::kConditional, conditions);
    auto findings = createFindingsWithoutGenotype();

    EXPECT_TRUE(shouldPlotReadVisualization(spec, findings));
}
