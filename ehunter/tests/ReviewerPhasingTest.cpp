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

#include "reviewer/Phasing.hh"

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "core/CountTable.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/GraphBlueprint.hh"
#include "io/JsonWriter.hh"
#include "io/RegionGraph.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "locus/VariantFindings.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

using graphtools::decodeGraphAlignment;
using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Helper to create a Read with given ID and sequence
Read makeRead(const std::string& fragId, MateNumber mate, const std::string& seq)
{
    ReadId readId(fragId, mate);
    return Read(readId, seq, false);
}

// Helper to create a Frag from a graph alignment
Frag makeFrag(
    const std::string& fragId, const Graph& graph, const std::string& readSeq, const std::string& readAlignStr,
    const std::string& mateSeq, const std::string& mateAlignStr)
{
    Read read = makeRead(fragId, MateNumber::kFirstMate, readSeq);
    GraphAlignment readAlign = decodeGraphAlignment(0, readAlignStr, &graph);
    ReadWithAlign readWithAlign(std::move(read), std::move(readAlign));

    Read mate = makeRead(fragId, MateNumber::kSecondMate, mateSeq);
    GraphAlignment mateAlign = decodeGraphAlignment(0, mateAlignStr, &graph);
    ReadWithAlign mateWithAlign(std::move(mate), std::move(mateAlign));

    return Frag(std::move(readWithAlign), std::move(mateWithAlign));
}

// Locus with two adjacent repeats, like HTT's CAG and CCG: "ATTCGA(C)*TT(G)*ATGTCG"
// Node 0: ATTCGA (left flank), node 1: C repeat, node 2: TT spacer, node 3: G repeat, node 4: ATGTCG (right flank)
LocusSpecification makeTwoRepeatLocusSpec()
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*TT(G)*ATGTCG"));
    LocusSpecification spec(
        "TWO_REPEAT_LOCUS", ChromType::kAutosome, { GenomicRegion(0, 100, 200) }, graph, {}, GenotyperParameters(10),
        false, {});
    const VariantClassification repeatClassification(VariantType::kRepeat, VariantSubtype::kCommonRepeat);
    spec.addVariantSpecification("C_REPEAT", repeatClassification, GenomicRegion(0, 106, 107), { 1 }, boost::none);
    spec.addVariantSpecification("G_REPEAT", repeatClassification, GenomicRegion(0, 109, 110), { 3 }, boost::none);
    return spec;
}

LocusFindings makeTwoRepeatFindings(
    std::vector<int> cRepeatAlleleSizes, std::vector<int> gRepeatAlleleSizes)
{
    LocusFindings findings;
    for (const auto& variantIdAndAlleleSizes :
         { std::make_pair("C_REPEAT", cRepeatAlleleSizes), std::make_pair("G_REPEAT", gRepeatAlleleSizes) })
    {
        findings.findingsForEachVariant[variantIdAndAlleleSizes.first] = std::make_unique<RepeatFindings>(
            CountTable(), CountTable(), CountTable(), AlleleCount::kTwo,
            RepeatGenotype(1, variantIdAndAlleleSizes.second), static_cast<GenotypeFilter>(0));
    }
    return findings;
}

// Runs the same pairing steps as runReviewerWorkflow and summarizes the chosen pairing
std::optional<RepeatAllelePhasing> phaseRepeatAlleles(
    const LocusSpecification& spec, const LocusFindings& findings, const FragById& fragById, int meanFragLen = 300)
{
    const ScoredDiplotypes scoredDiplotypes
        = scoreDiplotypes(fragById, getCandidateDiplotypes(meanFragLen, spec, findings));
    return summarizeRepeatAllelePhasing(
        spec, findings, scoredDiplotypes.front().first, compareTopTwoDiplotypes(fragById, scoredDiplotypes));
}

}  // namespace

// Test scoreDiplotypes with a single diplotype
TEST(ReviewerPhasing_ScoreDiplotypes, SingleDiplotype_ReturnsScored)
{
    // Create a simple graph: ATTCGA(C)*ATGTCG
    // Node 0: ATTCGA (left flank), Node 1: C (repeat), Node 2: ATGTCG (right flank)
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create a diplotype with a single path through nodes 0, 1, 1, 2 (3 copies of repeat)
    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Create fragments that align to this path
    FragById fragById;
    // Read spanning left flank and repeat
    Frag frag1 = makeFrag("frag1", graph, "ATTCGAC", "0[6M]1[1M]", "CATGTCG", "1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // The diplotype should have a positive score since fragments align
    EXPECT_GT(result[0].second, 0);
    EXPECT_EQ(diplotype.size(), result[0].first.size());
}

// Test scoreDiplotypes with multiple diplotypes sorted by score
TEST(ReviewerPhasing_ScoreDiplotypes, MultipleDiplotypes_SortedByScoreDescending)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create two diplotypes with different repeat counts
    // Path with 2 repeats
    Path path2repeats(&graph, 0, { 0, 1, 1, 2 }, 6);
    // Path with 3 repeats
    Path path3repeats(&graph, 0, { 0, 1, 1, 1, 2 }, 6);

    Diplotype diplotype2 = { path2repeats };
    Diplotype diplotype3 = { path3repeats };
    std::vector<Diplotype> diplotypes = { diplotype2, diplotype3 };

    // Create a fragment that aligns better to the 2-repeat path
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(2u, result.size());
    // Results should be sorted by score in descending order
    EXPECT_GE(result[0].second, result[1].second);
}

// Test scoreDiplotypes with tied scores maintains stable ordering
TEST(ReviewerPhasing_ScoreDiplotypes, TiedScores_StableOrdering)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create identical paths - they should have the same score
    Path path1(&graph, 0, { 0, 1, 1, 2 }, 6);
    Path path2(&graph, 0, { 0, 1, 1, 2 }, 6);

    Diplotype diplotype1 = { path1 };
    Diplotype diplotype2 = { path2 };
    std::vector<Diplotype> diplotypes = { diplotype1, diplotype2 };

    // Create a fragment
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(2u, result.size());
    // Both should have the same score
    EXPECT_EQ(result[0].second, result[1].second);
    // Note: scoreDiplotypes uses std::stable_sort, so tied diplotypes keep their input order
}

// Test scoreDiplotypes with no fragments - all diplotypes should score 0
TEST(ReviewerPhasing_ScoreDiplotypes, NoFragments_AllDiplotypesScoreZero)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Empty fragment map
    FragById fragById;

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // With no fragments, score should be 0
    EXPECT_EQ(0, result[0].second);
}

// Test scoreDiplotypes with heterozygous diplotype (two different haplotypes)
TEST(ReviewerPhasing_ScoreDiplotypes, HeterozygousDiplotype_BothHaplotypesScored)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    // Create a heterozygous diplotype with paths of different repeat counts
    Path pathShort(&graph, 0, { 0, 1, 2 }, 6);     // 1 repeat
    Path pathLong(&graph, 0, { 0, 1, 1, 1, 2 }, 6);  // 3 repeats

    Diplotype diplotype = { pathShort, pathLong };
    std::vector<Diplotype> diplotypes = { diplotype };

    // Create fragments
    FragById fragById;
    Frag frag1 = makeFrag("frag1", graph, "ATTCGAC", "0[6M]1[1M]", "CATGTCG", "1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    ScoredDiplotypes result = scoreDiplotypes(fragById, diplotypes);

    ASSERT_EQ(1u, result.size());
    // The diplotype should have a positive score
    EXPECT_GT(result[0].second, 0);
    // Verify it's still a heterozygous diplotype
    EXPECT_EQ(2u, result[0].first.size());
}

// Test scoreDiplotypes with multiple fragments
TEST(ReviewerPhasing_ScoreDiplotypes, MultipleFragments_ScoresAccumulate)
{
    Graph graph = makeRegionGraph(decodeFeaturesFromRegex("ATTCGA(C)*ATGTCG"));

    Path path(&graph, 0, { 0, 1, 1, 2 }, 6);
    Diplotype diplotype = { path };
    std::vector<Diplotype> diplotypes = { diplotype };

    FragById fragById;
    // Add first fragment
    Frag frag1 = makeFrag("frag1", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag1", std::move(frag1));

    // Get score with one fragment
    ScoredDiplotypes result1 = scoreDiplotypes(fragById, diplotypes);
    int score1 = result1[0].second;

    // Add second fragment
    Frag frag2 = makeFrag("frag2", graph, "ATTCGACC", "0[6M]1[1M]1[1M]", "CCATGTCG", "1[1M]1[1M]2[6M]");
    fragById.emplace("frag2", std::move(frag2));

    // Get score with two fragments
    ScoredDiplotypes result2 = scoreDiplotypes(fragById, diplotypes);
    int score2 = result2[0].second;

    // Score with two fragments should be greater than with one
    EXPECT_GT(score2, score1);
}

// A fragment spanning both repeats on the C=1, G=1 haplotype places the other alleles (C=3, G=2) together
TEST(ReviewerPhasing_RepeatAllelePhasing, FragmentSpanningBothRepeats_PairsShortAllelesTogether)
{
    const LocusSpecification spec = makeTwoRepeatLocusSpec();
    const LocusFindings findings = makeTwoRepeatFindings({ 1, 3 }, { 1, 2 });
    const Graph& graph = spec.regionGraph();

    FragById fragById;
    fragById.emplace(
        "frag1",
        makeFrag("frag1", graph, "ATTCGACTTGATGTCG", "0[6M]1[1M]2[2M]3[1M]4[6M]", "ATGTCG", "4[6M]"));

    const auto phasing = phaseRepeatAlleles(spec, findings, fragById);
    ASSERT_TRUE(phasing);
    const std::vector<std::map<std::string, int>> expectedHaplotypes
        = { { { "C_REPEAT", 1 }, { "G_REPEAT", 1 } }, { { "C_REPEAT", 3 }, { "G_REPEAT", 2 } } };
    EXPECT_EQ(expectedHaplotypes, phasing->repeatSizeByVariantIdOnEachHaplotype);
    EXPECT_EQ(2, phasing->numberOfPossiblePairings);
    EXPECT_EQ(1, phasing->fragmentsSupportingChosenPairingOverNextBest);
    EXPECT_EQ(0, phasing->fragmentsSupportingNextBestPairingOverChosen);
    EXPECT_FALSE(phasing->chosenPairingIsTiedWithNextBest);
}

// The opposite phase: the short C allele sits with the long G allele, which sorting alleles by size would hide
TEST(ReviewerPhasing_RepeatAllelePhasing, FragmentSpanningBothRepeats_PairsShortAlleleWithLongAllele)
{
    const LocusSpecification spec = makeTwoRepeatLocusSpec();
    const LocusFindings findings = makeTwoRepeatFindings({ 1, 3 }, { 1, 2 });
    const Graph& graph = spec.regionGraph();

    FragById fragById;
    fragById.emplace(
        "frag1",
        makeFrag("frag1", graph, "ATTCGACTTGGATGTCG", "0[6M]1[1M]2[2M]3[1M]3[1M]4[6M]", "ATGTCG", "4[6M]"));

    const auto phasing = phaseRepeatAlleles(spec, findings, fragById);
    ASSERT_TRUE(phasing);
    const std::vector<std::map<std::string, int>> expectedHaplotypes
        = { { { "C_REPEAT", 1 }, { "G_REPEAT", 2 } }, { { "C_REPEAT", 3 }, { "G_REPEAT", 1 } } };
    EXPECT_EQ(expectedHaplotypes, phasing->repeatSizeByVariantIdOnEachHaplotype);
    EXPECT_EQ(1, phasing->fragmentsSupportingChosenPairingOverNextBest);
    EXPECT_EQ(0, phasing->fragmentsSupportingNextBestPairingOverChosen);
    EXPECT_FALSE(phasing->chosenPairingIsTiedWithNextBest);
}

// Fragments that fit every haplotype equally well cannot tell the pairings apart, so the choice is a tie
TEST(ReviewerPhasing_RepeatAllelePhasing, NoFragmentSpansBothRepeats_ReportsTie)
{
    const LocusSpecification spec = makeTwoRepeatLocusSpec();
    const LocusFindings findings = makeTwoRepeatFindings({ 1, 3 }, { 1, 2 });
    const Graph& graph = spec.regionGraph();

    FragById fragById;
    fragById.emplace("frag1", makeFrag("frag1", graph, "ATTCGA", "0[6M]", "ATGTCG", "4[6M]"));

    const auto phasing = phaseRepeatAlleles(spec, findings, fragById);
    ASSERT_TRUE(phasing);
    EXPECT_EQ(2, phasing->numberOfPossiblePairings);
    EXPECT_EQ(0, phasing->fragmentsSupportingChosenPairingOverNextBest);
    EXPECT_EQ(0, phasing->fragmentsSupportingNextBestPairingOverChosen);
    EXPECT_TRUE(phasing->chosenPairingIsTiedWithNextBest);
}

// With one homozygous repeat there is only one way to pair the alleles
TEST(ReviewerPhasing_RepeatAllelePhasing, OneHomozygousRepeat_SinglePossiblePairing)
{
    const LocusSpecification spec = makeTwoRepeatLocusSpec();
    const LocusFindings findings = makeTwoRepeatFindings({ 2, 2 }, { 1, 2 });

    const auto phasing = phaseRepeatAlleles(spec, findings, FragById());
    ASSERT_TRUE(phasing);
    const std::vector<std::map<std::string, int>> expectedHaplotypes
        = { { { "C_REPEAT", 2 }, { "G_REPEAT", 1 } }, { { "C_REPEAT", 2 }, { "G_REPEAT", 2 } } };
    EXPECT_EQ(expectedHaplotypes, phasing->repeatSizeByVariantIdOnEachHaplotype);
    EXPECT_EQ(1, phasing->numberOfPossiblePairings);
    EXPECT_FALSE(phasing->chosenPairingIsTiedWithNextBest);
}

// Both C alleles exceed the path length cap, so the paths no longer show which haplotype carries which
TEST(ReviewerPhasing_RepeatAllelePhasing, BothAllelesCappedToSamePathLength_ReturnsEmpty)
{
    const LocusSpecification spec = makeTwoRepeatLocusSpec();
    const LocusFindings findings = makeTwoRepeatFindings({ 3, 4 }, { 1, 2 });

    EXPECT_FALSE(phaseRepeatAlleles(spec, findings, FragById(), 2));
}

TEST(ReviewerPhasing_RepeatAllelePhasing, EncodeJson_OmitsComparisonFieldsWhenOnlyOnePairing)
{
    RepeatAllelePhasing phasing;
    phasing.repeatSizeByVariantIdOnEachHaplotype
        = { { { "C_REPEAT", 2 }, { "G_REPEAT", 1 } }, { { "C_REPEAT", 2 }, { "G_REPEAT", 2 } } };
    phasing.numberOfPossiblePairings = 1;

    const nlohmann::json record = encodeRepeatAllelePhasing(phasing);
    EXPECT_EQ(1, record["NumberOfPossiblePairings"]);
    EXPECT_EQ(2, record["RepeatSizesOnEachHaplotype"][1]["G_REPEAT"]);
    EXPECT_FALSE(record.contains("FragmentsSupportingChosenPairingOverNextBest"));
    EXPECT_FALSE(record.contains("ChosenPairingIsTiedWithNextBest"));

    phasing.numberOfPossiblePairings = 2;
    phasing.fragmentsSupportingChosenPairingOverNextBest = 5;
    const nlohmann::json recordWithComparison = encodeRepeatAllelePhasing(phasing);
    EXPECT_EQ(5, recordWithComparison["FragmentsSupportingChosenPairingOverNextBest"]);
    EXPECT_EQ(0, recordWithComparison["FragmentsSupportingNextBestPairingOverChosen"]);
    EXPECT_EQ(false, recordWithComparison["ChosenPairingIsTiedWithNextBest"]);
}
