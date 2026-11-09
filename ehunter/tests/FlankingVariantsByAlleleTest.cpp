//
// ExpansionHunter
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

#include "reviewer/FlankingVariantsByAllele.hh"

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/GraphAlignmentOperations.hh"

#include "io/GraphBlueprint.hh"
#include "io/RegionGraph.hh"

using graphtools::Graph;
using graphtools::GraphAlignment;
using graphtools::Path;
using std::string;
using std::vector;
using namespace ehunter;
using namespace ehunter::reviewer;

namespace
{

// Left flank, a CAG repeat and right flank, at chr index 0: left flank [90, 100), right flank starting at
// 100 + 3 * <number of reference repeat units>
const string kLeftFlank = "TTTCGATGCA";
const string kRightFlank = "GAAAACTTGC";
const int64_t kLeftFlankStart = 90;
const int kReferenceRepeatUnits = 2;
const int64_t kRightFlankStart = kLeftFlankStart + 10 + 3 * kReferenceRepeatUnits;

Graph makeLocusGraph() { return makeRegionGraph(decodeFeaturesFromRegex(kLeftFlank + "(CAG)*" + kRightFlank)); }

NodeToRegionAssociation makeNodeRegions()
{
    NodeToRegionAssociation regions;
    regions.emplace(0, GenomicRegion(0, kLeftFlankStart, kLeftFlankStart + 10));
    regions.emplace(1, GenomicRegion(0, kLeftFlankStart + 10, kRightFlankStart));
    regions.emplace(2, GenomicRegion(0, kRightFlankStart, kRightFlankStart + 10));
    return regions;
}

Path makeHaplotypePath(const Graph& graph, int repeatUnits)
{
    vector<graphtools::NodeId> nodes = { 0 };
    nodes.insert(nodes.end(), repeatUnits, 1);
    nodes.push_back(2);
    return Path(&graph, 0, nodes, static_cast<int32_t>(kRightFlank.size()));
}

// One read of a read pair: its bases and its alignment, given by the offset in its first node and a graph CIGAR
struct TestRead
{
    string bases;
    int startOffset;
    string graphCigar;
};

// Read pairs and how each one is placed on the haplotype paths, built up by the tests
struct TestReadPairs
{
    FragById fragById;
    FragPathAlignsById fragPathAlignsById;
    std::set<string> fragIdsPlacedOnOnePathByAlignmentScore;

    // Adds a read pair whose alignments fall on the given paths (one alignment per path, all with the same reads).
    // A pair placed on one path got there by its alignment score unless `placedByFragmentLength` is set, in which
    // case it scored equally on every path and the fragment-length filter picked this one.
    void add(
        const string& fragId, const TestRead& read, const TestRead& mate, const Graph& graph,
        const vector<std::pair<const Path*, int>>& pathsAndIndexes, bool placedByFragmentLength = false)
    {
        if (pathsAndIndexes.size() == 1 && !placedByFragmentLength)
        {
            fragIdsPlacedOnOnePathByAlignmentScore.insert(fragId);
        }
        auto makeReadWithAlign = [&](const TestRead& testRead, MateNumber mateNumber) {
            return ReadWithAlign(
                Read(ReadId(fragId, mateNumber), testRead.bases, false),
                decodeGraphAlignment(testRead.startOffset, testRead.graphCigar, &graph));
        };
        fragById.emplace(
            fragId, Frag(makeReadWithAlign(read, MateNumber::kFirstMate), makeReadWithAlign(mate, MateNumber::kSecondMate)));

        for (const auto& pathAndIndex : pathsAndIndexes)
        {
            const Path& path = *pathAndIndex.first;
            // Where the read's first node first appears on the path; enough for these tests, whose reads start in a
            // flank or in the first repeat unit
            auto startIndexOnPath = [&path](const GraphAlignment& align) {
                const auto& pathNodes = path.nodeIds();
                const auto firstNode = align.path().nodeIds().front();
                return static_cast<int>(std::find(pathNodes.begin(), pathNodes.end(), firstNode) - pathNodes.begin());
            };
            auto readAlign = std::make_shared<GraphAlignment>(
                decodeGraphAlignment(read.startOffset, read.graphCigar, &graph));
            auto mateAlign = std::make_shared<GraphAlignment>(
                decodeGraphAlignment(mate.startOffset, mate.graphCigar, &graph));
            ReadPathAlign readPathAlign(path, pathAndIndex.second, startIndexOnPath(*readAlign), readAlign);
            ReadPathAlign matePathAlign(path, pathAndIndex.second, startIndexOnPath(*mateAlign), mateAlign);
            fragPathAlignsById[fragId].emplace_back(std::move(readPathAlign), std::move(matePathAlign));
        }
    }

    FragAssignment assignFirstAlignment() const
    {
        vector<string> fragIds;
        vector<int> alignIndexes;
        for (const auto& fragIdAndAligns : fragPathAlignsById)
        {
            fragIds.push_back(fragIdAndAligns.first);
            alignIndexes.push_back(0);
        }
        return FragAssignment(std::move(fragIds), std::move(alignIndexes));
    }
};

// A read crossing the right flank's 4 A's, with nothing else in it
const TestRead kRightFlankRead { kRightFlank, 0, "2[10M]" };

// By default the reference bases next to the flanks are those of the CAG repeat: C after the left flank, G before the
// right flank
std::optional<vector<AlleleFlankingVariants>> findVariants(
    const Graph& graph, const RepeatGenotype& genotype, const Diplotype& diplotype, const TestReadPairs& readPairs,
    char referenceBaseBeforeRightFlank = 'G', char referenceBaseAfterLeftFlank = 'C')
{
    return findFlankingVariantsByAllele(
        graph, makeNodeRegions(), referenceBaseAfterLeftFlank, referenceBaseBeforeRightFlank, 1, genotype, diplotype,
        readPairs.fragById, readPairs.fragPathAlignsById, readPairs.assignFirstAlignment(),
        readPairs.fragIdsPlacedOnOnePathByAlignmentScore);
}

const FlankingVariantReadCounts* findCounts(const AlleleFlankingVariants& allele, const FlankingVariant& variant)
{
    for (const auto& counts : allele.variants)
    {
        if (counts.variant == variant)
        {
            return &counts;
        }
    }
    return nullptr;
}

}  // namespace

TEST(ShiftDeletionLeftTest, MovesDeletionToStartOfHomopolymer)
{
    EXPECT_EQ(1, shiftDeletionLeft("CAAAAT", 4, 5));
    EXPECT_EQ(1, shiftDeletionLeft("CAAAAT", 1, 2));
}

TEST(ShiftDeletionLeftTest, MovesDeletionAcrossDinucleotideRepeat)
{
    // Deleting "AC" at [5, 7) of GACACACT is the same as deleting it at [1, 3)
    EXPECT_EQ(1, shiftDeletionLeft("GACACACT", 5, 7));
}

TEST(ShiftDeletionLeftTest, LeavesDeletionWithoutRepeatedBasesInPlace)
{
    EXPECT_EQ(3, shiftDeletionLeft("ACGTAC", 3, 4));
}

TEST(ShiftInsertionLeftTest, MovesInsertionToStartOfHomopolymer)
{
    const auto shifted = shiftInsertionLeft("CAAAT", 4, "A");
    EXPECT_EQ(1, shifted.first);
    EXPECT_EQ("A", shifted.second);
}

TEST(ShiftInsertionLeftTest, RotatesInsertedBasesWhileShifting)
{
    // Inserting "CA" before offset 4 of GACACT gives GACACACT, the same as inserting "AC" before offset 1
    const auto shifted = shiftInsertionLeft("GACACT", 4, "CA");
    EXPECT_EQ(1, shifted.first);
    EXPECT_EQ("AC", shifted.second);
}

TEST(FindFlankingVariantsByAlleleTest, PhasesVariantsWithHeterozygousRepeatAlleles)
{
    const Graph graph = makeLocusGraph();
    const Path shortPath = makeHaplotypePath(graph, 2);
    const Path longPath = makeHaplotypePath(graph, 4);
    // The long allele's path comes first, to check that entries are ordered by allele size
    const Diplotype diplotype = { longPath, shortPath };
    const RepeatGenotype genotype(3, { 2, 4 });

    TestReadPairs readPairs;
    // Short allele: an SNV at left flank offset 3 (C>A), in two read pairs whose mates cover the right flank
    const TestRead shortAlleleRead { "TTTAGATGCA" "CAGCAG" "GAAAA", 0, "0[3M1X6M]1[3M]1[3M]2[5M]" };
    readPairs.add("short1", shortAlleleRead, kRightFlankRead, graph, { { &shortPath, 1 } });
    readPairs.add("short2", shortAlleleRead, kRightFlankRead, graph, { { &shortPath, 1 } });
    // Long allele: reads that start at left flank offset 4 (so don't cover the SNV), and mates with a T inserted
    // after right flank offset 5
    const TestRead longAlleleRead { "GATGCA" "CAGCAGCAGCAG" "GA", 4, "0[6M]1[3M]1[3M]1[3M]1[3M]2[2M]" };
    const TestRead longAlleleMate { "GAAAACTTTGC", 0, "2[6M1I4M]" };
    readPairs.add("long1", longAlleleRead, longAlleleMate, graph, { { &longPath, 0 } });
    readPairs.add("long2", longAlleleRead, longAlleleMate, graph, { { &longPath, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(2u, alleles->size());
    const AlleleFlankingVariants& shortAllele = alleles->at(0);
    const AlleleFlankingVariants& longAllele = alleles->at(1);
    EXPECT_EQ(2, shortAllele.alleleSize);
    EXPECT_EQ(4, longAllele.alleleSize);

    const FlankingVariant snv { kLeftFlankStart + 3, "C", "A" };
    // The T is inserted between the C at right flank offset 5 and the TT at offsets 6-7. It can't shift further
    // left, so its anchor base is that C.
    const FlankingVariant insertion { kRightFlankStart + 5, "C", "CT" };

    // The same variants, in the same order, under both alleles
    ASSERT_EQ(2u, shortAllele.variants.size());
    ASSERT_EQ(2u, longAllele.variants.size());
    EXPECT_EQ(snv, shortAllele.variants[0].variant);
    EXPECT_EQ(insertion, shortAllele.variants[1].variant);
    EXPECT_EQ(snv, longAllele.variants[0].variant);
    EXPECT_EQ(insertion, longAllele.variants[1].variant);

    EXPECT_EQ(0, shortAllele.variants[0].refReads);
    EXPECT_EQ(2, shortAllele.variants[0].altReads);
    // The long allele's reads don't reach the SNV
    EXPECT_EQ(0, longAllele.variants[0].refReads);
    EXPECT_EQ(0, longAllele.variants[0].altReads);

    // Only the short allele's mates cross the insertion point; its reads end at right flank offset 4
    EXPECT_EQ(2, shortAllele.variants[1].refReads);
    EXPECT_EQ(0, shortAllele.variants[1].altReads);
    EXPECT_EQ(0, longAllele.variants[1].refReads);
    EXPECT_EQ(2, longAllele.variants[1].altReads);

    ASSERT_TRUE(shortAllele.haplotypeRegion.has_value());
    EXPECT_EQ(GenomicRegion(0, kLeftFlankStart, kRightFlankStart + 10), *shortAllele.haplotypeRegion);
    ASSERT_TRUE(longAllele.haplotypeRegion.has_value());
    EXPECT_EQ(GenomicRegion(0, kLeftFlankStart + 4, kRightFlankStart + 10), *longAllele.haplotypeRegion);
}

TEST(FindFlankingVariantsByAlleleTest, IgnoresReadPairsThatFitBothAlleles)
{
    const Graph graph = makeLocusGraph();
    const Path shortPath = makeHaplotypePath(graph, 2);
    const Path longPath = makeHaplotypePath(graph, 4);
    const Diplotype diplotype = { shortPath, longPath };
    const RepeatGenotype genotype(3, { 2, 4 });

    TestReadPairs readPairs;
    // Two read pairs that lie in the right flank and fit both paths, with an SNV at right flank offset 6 (T>C)
    const TestRead flankOnlyRead { "GAAAACCTGC", 0, "2[6M1X3M]" };
    readPairs.add("both1", flankOnlyRead, flankOnlyRead, graph, { { &shortPath, 0 }, { &longPath, 1 } });
    readPairs.add("both2", flankOnlyRead, flankOnlyRead, graph, { { &shortPath, 0 }, { &longPath, 1 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(2u, alleles->size());
    EXPECT_TRUE(alleles->at(0).variants.empty());
    EXPECT_TRUE(alleles->at(1).variants.empty());
    EXPECT_FALSE(alleles->at(0).haplotypeRegion.has_value());
    EXPECT_FALSE(alleles->at(1).haplotypeRegion.has_value());
}

TEST(FindFlankingVariantsByAlleleTest, IgnoresReadPairsPlacedOnAnAlleleByFragmentLengthAlone)
{
    const Graph graph = makeLocusGraph();
    const Path shortPath = makeHaplotypePath(graph, 2);
    const Path longPath = makeHaplotypePath(graph, 4);
    const Diplotype diplotype = { shortPath, longPath };
    const RepeatGenotype genotype(3, { 2, 4 });

    TestReadPairs readPairs;
    // Pairs with one read in each flank fit both paths equally well; the fragment-length filter put these two on the
    // long allele's path. They show an SNV at right flank offset 6 (T>C).
    const TestRead leftFlankRead { kLeftFlank, 0, "0[10M]" };
    const TestRead snvRead { "GAAAACCTGC", 0, "2[6M1X3M]" };
    readPairs.add("byLength1", leftFlankRead, snvRead, graph, { { &longPath, 1 } }, true);
    readPairs.add("byLength2", leftFlankRead, snvRead, graph, { { &longPath, 1 } }, true);

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    EXPECT_TRUE(alleles->at(0).variants.empty());
    EXPECT_FALSE(alleles->at(1).haplotypeRegion.has_value());
}

TEST(FindFlankingVariantsByAlleleTest, HomozygousGenotypeHasOneEntryWithEveryReadPair)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path, path };
    const RepeatGenotype genotype(3, { 3, 3 });

    TestReadPairs readPairs;
    const TestRead snvRead { "GAAAACCTGC", 0, "2[6M1X3M]" };
    // The pairs fit both (identical) paths, which is how REViewer places them for a homozygous genotype
    readPairs.add("alt1", snvRead, snvRead, graph, { { &path, 0 }, { &path, 1 } });
    readPairs.add("ref1", kRightFlankRead, kRightFlankRead, graph, { { &path, 0 }, { &path, 1 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->size());
    EXPECT_EQ(3, alleles->at(0).alleleSize);
    ASSERT_EQ(1u, alleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 6, "T", "C" }), alleles->at(0).variants[0].variant);
    EXPECT_EQ(2, alleles->at(0).variants[0].refReads);
    EXPECT_EQ(2, alleles->at(0).variants[0].altReads);
}

TEST(FindFlankingVariantsByAlleleTest, ReadsCantBeSplitWhenBothAllelePathsHaveTheSameLength)
{
    const Graph graph = makeLocusGraph();
    // Both alleles were capped to the same path length
    const Path cappedPath = makeHaplotypePath(graph, 5);
    const Diplotype diplotype = { cappedPath, cappedPath };
    const RepeatGenotype genotype(3, { 200, 300 });

    TestReadPairs readPairs;
    readPairs.add("frag1", kRightFlankRead, kRightFlankRead, graph, { { &cappedPath, 0 } });

    EXPECT_FALSE(findVariants(graph, genotype, diplotype, readPairs).has_value());
}

TEST(FindFlankingVariantsByAlleleTest, ListsOnlyVariantsShownByAtLeastTwoReadsOfOneAllele)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    // A single read shows T>C at right flank offset 6
    readPairs.add("frag1", TestRead { "GAAAACCTGC", 0, "2[6M1X3M]" }, kRightFlankRead, graph, { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->size());
    EXPECT_TRUE(alleles->at(0).variants.empty());
}

TEST(FindFlankingVariantsByAlleleTest, LowQualityBasesShowNeitherAllele)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    const TestRead highQualityAltRead { "GAAAACCTGC", 0, "2[6M1X3M]" };
    // The same C, but low quality (lower case), and a low-quality reference T
    const TestRead lowQualityAltRead { "GAAAACcTGC", 0, "2[6M1X3M]" };
    const TestRead lowQualityRefRead { "GAAAACtTGC", 0, "2[10M]" };
    readPairs.add("frag1", highQualityAltRead, highQualityAltRead, graph, { { &path, 0 } });
    readPairs.add("frag2", lowQualityAltRead, lowQualityRefRead, graph, { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->at(0).variants.size());
    EXPECT_EQ(0, alleles->at(0).variants[0].refReads);
    EXPECT_EQ(2, alleles->at(0).variants[0].altReads);
}

TEST(FindFlankingVariantsByAlleleTest, MergesADeletionPlacedAnywhereInAHomopolymer)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    // One A deleted from the right flank's AAAA (offsets 1-4), placed at offset 4 by one read and offset 1 by another
    readPairs.add(
        "frag1", TestRead { "GAAACTTGC", 0, "2[4M1D5M]" }, TestRead { "GAAACTTGC", 0, "2[1M1D8M]" }, graph,
        { { &path, 0 } });
    // A read without the deletion shows the REF allele; a read with two A's deleted shows neither allele
    readPairs.add(
        "frag2", kRightFlankRead, TestRead { "GAACTTGC", 0, "2[3M2D5M]" }, graph, { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    const FlankingVariant oneADeleted { kRightFlankStart, "GA", "G" };
    const FlankingVariantReadCounts* counts = findCounts(alleles->at(0), oneADeleted);
    ASSERT_NE(nullptr, counts);
    EXPECT_EQ(1, counts->refReads);
    EXPECT_EQ(2, counts->altReads);
    // The two-A deletion has a single read, so it is not listed
    EXPECT_EQ(1u, alleles->at(0).variants.size());
}

TEST(FindFlankingVariantsByAlleleTest, MergesAChangeThatReadsSplitIntoDifferentPieces)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    // Both reads replace the AAC at right flank offsets 3-5 (GAA[AAC]TTGC) with TTGG. One read is aligned as 3
    // mismatches followed by an inserted G, the other as an inserted T followed by 3 mismatches.
    readPairs.add(
        "frag1", TestRead { "GAATTGGTTGC", 0, "2[3M3X1I4M]" }, TestRead { "GAATTGGTTGC", 0, "2[3M1I3X4M]" }, graph,
        { { &path, 0 } });
    readPairs.add("frag2", kRightFlankRead, kRightFlankRead, graph, { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 3, "AAC", "TTGG" }), alleles->at(0).variants[0].variant);
    EXPECT_EQ(2, alleles->at(0).variants[0].refReads);
    EXPECT_EQ(2, alleles->at(0).variants[0].altReads);
}

TEST(FindFlankingVariantsByAlleleTest, MergesSnvsLessThanFiveBasesApart)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    // SNVs at right flank offsets 1 and 5 (4 bases apart) become one variant
    const TestRead closeSnvsRead { "GTAAAGTTGC", 0, "2[1M1X3M1X4M]" };
    readPairs.add("close", closeSnvsRead, closeSnvsRead, graph, { { &path, 0 } });
    const auto closeAlleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(closeAlleles.has_value());
    ASSERT_EQ(1u, closeAlleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 1, "AAAAC", "TAAAG" }), closeAlleles->at(0).variants[0].variant);

    // SNVs at right flank offsets 1 and 6 (5 bases apart) stay separate
    TestReadPairs farReadPairs;
    const TestRead farSnvsRead { "GTAAACCTGC", 0, "2[1M1X4M1X3M]" };
    farReadPairs.add("far", farSnvsRead, farSnvsRead, graph, { { &path, 0 } });
    const auto farAlleles = findVariants(graph, genotype, diplotype, farReadPairs);
    ASSERT_TRUE(farAlleles.has_value());
    ASSERT_EQ(2u, farAlleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 1, "A", "T" }), farAlleles->at(0).variants[0].variant);
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 6, "T", "C" }), farAlleles->at(0).variants[1].variant);
}

TEST(FindFlankingVariantsByAlleleTest, LowQualityMismatchNextToAnSnvDoesNotMergeWithIt)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 3);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 3 });

    TestReadPairs readPairs;
    const TestRead snvRead { "GAAAACCTGC", 0, "2[6M1X3M]" };
    // The same SNV, with a low-quality mismatch 2 bases away (offset 8, G read as t)
    const TestRead snvWithLowQualityNeighborRead { "GAAAACCTtC", 0, "2[6M1X1M1X1M]" };
    readPairs.add("frag1", snvRead, snvRead, graph, { { &path, 0 } });
    readPairs.add("frag2", snvWithLowQualityNeighborRead, snvWithLowQualityNeighborRead, graph, { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs);
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart + 6, "T", "C" }), alleles->at(0).variants[0].variant);
    EXPECT_EQ(0, alleles->at(0).variants[0].refReads);
    EXPECT_EQ(4, alleles->at(0).variants[0].altReads);
}

TEST(FindFlankingVariantsByAlleleTest, LeavesOutIndelsThatCouldEquallyWellBeInTheRepeat)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 2);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 2 });

    // The right flank's first G deleted; the repeat ends in G, so the same deletion could be the repeat's last G
    TestReadPairs deletionReadPairs;
    const TestRead deletionRead { "CAGCAG" "AAAACTTGC", 0, "1[3M]1[3M]2[1D9M]" };
    deletionReadPairs.add("frag1", deletionRead, deletionRead, graph, { { &path, 0 } });
    const auto deletionAlleles = findVariants(graph, genotype, diplotype, deletionReadPairs, 'G');
    ASSERT_TRUE(deletionAlleles.has_value());
    EXPECT_TRUE(deletionAlleles->at(0).variants.empty());
    // If the base before the right flank were a T, the deletion would be reported
    const auto deletionAllelesNextToT = findVariants(graph, genotype, diplotype, deletionReadPairs, 'T');
    ASSERT_TRUE(deletionAllelesNextToT.has_value());
    ASSERT_EQ(1u, deletionAllelesNextToT->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart - 1, "TG", "T" }), deletionAllelesNextToT->at(0).variants[0].variant);

    // A C inserted at the end of the left flank; the repeat starts with C, so it could be the repeat's first base
    TestReadPairs insertionReadPairs;
    const TestRead insertionRead { kLeftFlank + "C" + "CAGCAG", 0, "0[10M1I]1[3M]1[3M]" };
    insertionReadPairs.add("frag1", insertionRead, insertionRead, graph, { { &path, 0 } });
    const auto insertionAlleles = findVariants(graph, genotype, diplotype, insertionReadPairs, 'G', 'C');
    ASSERT_TRUE(insertionAlleles.has_value());
    EXPECT_TRUE(insertionAlleles->at(0).variants.empty());
    // If the base after the left flank were a G, the insertion would be reported, anchored on the flank's last A
    const auto insertionAllelesNextToG = findVariants(graph, genotype, diplotype, insertionReadPairs, 'G', 'G');
    ASSERT_TRUE(insertionAllelesNextToG.has_value());
    ASSERT_EQ(1u, insertionAllelesNextToG->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kLeftFlankStart + 9, "A", "AC" }), insertionAllelesNextToG->at(0).variants[0].variant);
}

TEST(FindFlankingVariantsByAlleleTest, AnchorsInsertionAtStartOfRightFlankOnBaseBeforeIt)
{
    const Graph graph = makeLocusGraph();
    const Path path = makeHaplotypePath(graph, 2);
    const Diplotype diplotype = { path };
    const RepeatGenotype genotype(3, { 2 });

    TestReadPairs readPairs;
    // A T inserted between the repeat and the right flank, placed at the start of the right flank's alignment
    const TestRead insertionRead { "CAGCAG" "T" "GAAAACTTGC", 0, "1[3M]1[3M]2[1I10M]" };
    readPairs.add("frag1", insertionRead, insertionRead, graph, { { &path, 0 } });
    // A read crossing from the repeat into the right flank shows the REF allele; a read that starts in the right
    // flank doesn't, since it can't tell whether something was inserted before it
    readPairs.add(
        "frag2", TestRead { "CAGCAG" "GAAAACTTGC", 0, "1[3M]1[3M]2[10M]" }, kRightFlankRead, graph,
        { { &path, 0 } });

    const auto alleles = findVariants(graph, genotype, diplotype, readPairs, 'G');
    ASSERT_TRUE(alleles.has_value());
    ASSERT_EQ(1u, alleles->at(0).variants.size());
    EXPECT_EQ((FlankingVariant { kRightFlankStart - 1, "G", "GT" }), alleles->at(0).variants[0].variant);
    EXPECT_EQ(1, alleles->at(0).variants[0].refReads);
    EXPECT_EQ(2, alleles->at(0).variants[0].altReads);

    // Without the base before the right flank, the insertion can't be written in VCF form, so it is left out
    const auto allelesWithoutAnchor = findVariants(graph, genotype, diplotype, readPairs, 'N');
    ASSERT_TRUE(allelesWithoutAnchor.has_value());
    EXPECT_TRUE(allelesWithoutAnchor->at(0).variants.empty());
}
