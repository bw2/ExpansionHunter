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

#include "reviewer/GenotypePaths.hh"

#include "gtest/gtest.h"

#include "core/CountTable.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantFindings.hh"

using graphtools::Graph;
using graphtools::Path;
using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Helper to create a LocusSpecification with a single STR for testing
// Graph structure: "ATTCGA(CAG)*ATGTCG"
// Node 0: ATTCGA (left flank, 6 bases)
// Node 1: CAG (repeat unit, 3 bases)
// Node 2: ATGTCG (right flank, 6 bases)
LocusSpecification createTestLocusSpec(const std::string& variantId = "TEST_VARIANT")
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(CAG)*ATGTCG"));
    GenotyperParameters params(10);

    LocusSpecification spec(
        "TEST_LOCUS",
        ChromType::kAutosome,
        {GenomicRegion(0, 100, 200)},
        graph,
        {},  // referenceRegions
        params,
        false,  // useRFC1MotifAnalysis
        {});    // plotConditions

    // Add variant specification for the repeat node
    spec.addVariantSpecification(
        variantId,
        VariantClassification(VariantType::kRepeat, VariantSubtype::kCommonRepeat),
        GenomicRegion(0, 106, 106),  // Reference locus for the repeat
        {1},  // Node 1 is the repeat node
        boost::none);

    return spec;
}

// Helper to create LocusFindings with a given genotype
LocusFindings createFindingsWithGenotype(
    const std::string& variantId,
    int shortAllele,
    int longAllele,
    int repeatUnitLen = 3)
{
    LocusFindings findings;

    CountTable spanningCounts;
    CountTable flankingCounts;
    CountTable inrepeatCounts;

    RepeatGenotype genotype(repeatUnitLen, {shortAllele, longAllele});

    auto repeatFindings = std::make_unique<RepeatFindings>(
        spanningCounts,
        flankingCounts,
        inrepeatCounts,
        AlleleCount::kTwo,
        genotype,
        static_cast<GenotypeFilter>(0));

    findings.findingsForEachVariant[variantId] = std::move(repeatFindings);
    return findings;
}

// Helper to create LocusFindings with a haploid (single-allele) genotype
LocusFindings createFindingsWithHaploidGenotype(
    const std::string& variantId,
    int alleleSize,
    int repeatUnitLen = 3)
{
    LocusFindings findings;

    CountTable spanningCounts;
    CountTable flankingCounts;
    CountTable inrepeatCounts;

    RepeatGenotype genotype(repeatUnitLen, {alleleSize});

    auto repeatFindings = std::make_unique<RepeatFindings>(
        spanningCounts,
        flankingCounts,
        inrepeatCounts,
        AlleleCount::kOne,
        genotype,
        static_cast<GenotypeFilter>(0));

    findings.findingsForEachVariant[variantId] = std::move(repeatFindings);
    return findings;
}

// Helper to create LocusFindings without a genotype for a variant
LocusFindings createFindingsWithoutGenotype(const std::string& variantId)
{
    LocusFindings findings;

    CountTable spanningCounts;
    CountTable flankingCounts;
    CountTable inrepeatCounts;

    auto repeatFindings = std::make_unique<RepeatFindings>(
        spanningCounts,
        flankingCounts,
        inrepeatCounts,
        AlleleCount::kTwo,
        boost::none,  // No genotype
        static_cast<GenotypeFilter>(0));

    findings.findingsForEachVariant[variantId] = std::move(repeatFindings);
    return findings;
}

}  // namespace

//
// Tests for getCandidateDiplotypes with homozygous genotype
//

TEST(GetCandidateDiplotypes, HomozygousGenotype_SingleDiplotypeWithSamePathTwice)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 3, 3);  // Homozygous 3/3
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);

    // Should return exactly one diplotype for homozygous genotype
    ASSERT_EQ(1u, diplotypes.size());

    // The diplotype should have two paths that are identical
    const Diplotype& diplotype = diplotypes[0];
    ASSERT_EQ(2u, diplotype.size());
    EXPECT_EQ(diplotype[0], diplotype[1]);

    // Each path should have nodes: 0 -> 1 -> 1 -> 1 -> 2 (3 repeat units)
    const Path& path = diplotype[0];
    const auto& nodeIds = path.nodeIds();
    // Expected: node 0 (left flank), node 1 x3 (repeats), node 2 (right flank)
    EXPECT_EQ(5u, nodeIds.size());
    EXPECT_EQ(0u, nodeIds[0]);
    EXPECT_EQ(1u, nodeIds[1]);
    EXPECT_EQ(1u, nodeIds[2]);
    EXPECT_EQ(1u, nodeIds[3]);
    EXPECT_EQ(2u, nodeIds[4]);
}

//
// Tests for getCandidateDiplotypes with heterozygous genotype
//

TEST(GetCandidateDiplotypes, HeterozygousGenotype_SingleDiplotypeWithTwoDifferentPaths)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 2, 4);  // Heterozygous 2/4
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);

    // Should return exactly one diplotype for heterozygous genotype (paths are ordered)
    ASSERT_EQ(1u, diplotypes.size());

    const Diplotype& diplotype = diplotypes[0];
    ASSERT_EQ(2u, diplotype.size());

    // The two paths should be different (compare by node IDs since Path lacks != operator)
    EXPECT_TRUE(diplotype[0].nodeIds() != diplotype[1].nodeIds());

    // One path should have 2 repeat units, the other should have 4
    const Path& path0 = diplotype[0];
    const Path& path1 = diplotype[1];

    // Count repeat node occurrences in each path
    int repeatCount0 = 0;
    int repeatCount1 = 0;
    for (auto nodeId : path0.nodeIds())
    {
        if (nodeId == 1)
            ++repeatCount0;
    }
    for (auto nodeId : path1.nodeIds())
    {
        if (nodeId == 1)
            ++repeatCount1;
    }

    // One should have 2, the other should have 4
    EXPECT_TRUE((repeatCount0 == 2 && repeatCount1 == 4) || (repeatCount0 == 4 && repeatCount1 == 2));
}

//
// Tests for getCandidateDiplotypes with zero repeat count
//

TEST(GetCandidateDiplotypes, ZeroRepeatCount_PathWithNoRepeatNodes)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 0, 2);  // 0/2 genotype
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);

    ASSERT_EQ(1u, diplotypes.size());

    const Diplotype& diplotype = diplotypes[0];
    ASSERT_EQ(2u, diplotype.size());

    // Check that one path has 0 repeat nodes and one has 2
    int repeatCount0 = 0;
    int repeatCount1 = 0;
    for (auto nodeId : diplotype[0].nodeIds())
    {
        if (nodeId == 1)
            ++repeatCount0;
    }
    for (auto nodeId : diplotype[1].nodeIds())
    {
        if (nodeId == 1)
            ++repeatCount1;
    }

    EXPECT_TRUE((repeatCount0 == 0 && repeatCount1 == 2) || (repeatCount0 == 2 && repeatCount1 == 0));
}

//
// Tests for getCandidateDiplotypes with repeat sizes capped by fragment length
//

TEST(GetCandidateDiplotypes, LargeRepeatSize_CappedByFragmentLength)
{
    auto spec = createTestLocusSpec();
    // Genotype with very large repeat count (500), should be capped by meanFragLen
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 10, 500);
    int meanFragLen = 50;  // Small fragment length to cap the repeat

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);

    ASSERT_EQ(1u, diplotypes.size());

    const Diplotype& diplotype = diplotypes[0];

    // Find the maximum repeat count in the diplotype
    int maxRepeatCount = 0;
    for (const auto& path : diplotype)
    {
        int repeatCount = 0;
        for (auto nodeId : path.nodeIds())
        {
            if (nodeId == 1)
                ++repeatCount;
        }
        maxRepeatCount = std::max(maxRepeatCount, repeatCount);
    }

    // The larger repeat should be capped to meanFragLen (50)
    EXPECT_LE(maxRepeatCount, meanFragLen);
}

//
// Tests for getCandidateDiplotypes with haploid genotype
//

TEST(GetCandidateDiplotypes, HaploidGenotype_SinglePathInDiplotype)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithHaploidGenotype("TEST_VARIANT", 5);
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);

    ASSERT_EQ(1u, diplotypes.size());

    const Diplotype& diplotype = diplotypes[0];
    // Haploid genotype produces a diplotype with a single path
    ASSERT_EQ(1u, diplotype.size());

    // Count repeat nodes
    int repeatCount = 0;
    for (auto nodeId : diplotype[0].nodeIds())
    {
        if (nodeId == 1)
            ++repeatCount;
    }
    EXPECT_EQ(5, repeatCount);
}

//
// Tests for getCandidateDiplotypes with missing genotype
//

TEST(GetCandidateDiplotypes, NoGenotype_ThrowsException)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithoutGenotype("TEST_VARIANT");
    int meanFragLen = 300;

    // The function should throw when genotype is missing
    EXPECT_THROW(getCandidateDiplotypes(meanFragLen, spec, findings), std::runtime_error);
}

//
// Tests for getCandidateDiplotypes with missing variant findings
//

TEST(GetCandidateDiplotypes, MissingVariantFindings_ThrowsException)
{
    auto spec = createTestLocusSpec();
    LocusFindings findings;  // Empty findings, no variant data
    int meanFragLen = 300;

    // Should throw because variant findings are not found
    EXPECT_THROW(getCandidateDiplotypes(meanFragLen, spec, findings), std::runtime_error);
}

//
// Tests for Diplotype output operator
//

TEST(DiplotypeOutputStream, HomozygousDiplotype_OutputsCorrectFormat)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 2, 2);
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);
    ASSERT_EQ(1u, diplotypes.size());

    std::ostringstream oss;
    oss << diplotypes[0];
    std::string output = oss.str();

    // Output should contain the repeat motif and count
    EXPECT_FALSE(output.empty());
    // The output format includes (LF)(CAG){2}(RF)/(LF)(CAG){2}(RF) for homozygous 2/2
    EXPECT_NE(std::string::npos, output.find("CAG"));
}

TEST(DiplotypeOutputStream, HeterozygousDiplotype_OutputsCorrectFormat)
{
    auto spec = createTestLocusSpec();
    auto findings = createFindingsWithGenotype("TEST_VARIANT", 1, 3);
    int meanFragLen = 300;

    auto diplotypes = getCandidateDiplotypes(meanFragLen, spec, findings);
    ASSERT_EQ(1u, diplotypes.size());

    std::ostringstream oss;
    oss << diplotypes[0];
    std::string output = oss.str();

    // Output should contain the repeat motif
    EXPECT_FALSE(output.empty());
    EXPECT_NE(std::string::npos, output.find("CAG"));
    // Heterozygous format includes a '/' separator
    EXPECT_NE(std::string::npos, output.find("/"));
}
