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

#include "locus/MotifComposition.hh"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <htslib/sam.h>

#include "gtest/gtest.h"

#include "core/Read.hh"
#include "graphutils/SequenceOperations.hh"

using namespace ehunter;
using namespace ehunter::motifcomposition;
using std::string;
using std::vector;

namespace
{

uint32_t cigarOp(int length, int op)
{
    return (static_cast<uint32_t>(length) << BAM_CIGAR_SHIFT) | static_cast<uint32_t>(op);
}

string repeatMotif(const string& motif, int count)
{
    string sequence;
    for (int index = 0; index != count; ++index)
    {
        sequence += motif;
    }
    return sequence;
}

// A CAG x 10 repeat at [100, 130) on contig 0, with non-repetitive flanks.
const string kLeftFlank = "ATCGATTGCATGCAATGCCG"; // reference [80, 100)
const string kRightFlank = "TTAGGCTAACGTTGACCTAG"; // reference [130, 150)
const int kReadLength = 150;
const int kRegionExtensionLength = 1000;

MotifCompositionLocus makeCagLocus()
{
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.referenceRepeatStart = 100;
    locus.referenceRepeatEnd = 130;
    locus.catalogMotif = "CAG";
    locus.referenceRepeatSequence = repeatMotif("CAG", 10);
    locus.meanFragmentLength = 300;
    return locus;
}

FullRead makeMappedRead(
    const string& name, MateNumber mateNumber, const string& sequence, int64_t pos, const vector<uint32_t>& cigar,
    bool isReversed = false, int mapq = 60)
{
    Read read(ReadId(name, mateNumber), sequence, isReversed);
    LinearAlignmentStats stats;
    stats.chromId = 0;
    stats.pos = static_cast<int32_t>(pos);
    stats.mapq = mapq;
    stats.isPaired = true;
    stats.isMapped = true;
    stats.isMateMapped = true;
    stats.cigar = cigar;
    return FullRead(std::move(read), std::move(stats));
}

FullRead makeUnmappedRead(const string& name, MateNumber mateNumber, const string& sequence, int64_t pos)
{
    Read read(ReadId(name, mateNumber), sequence, false);
    LinearAlignmentStats stats;
    stats.chromId = 0;
    stats.pos = static_cast<int32_t>(pos);
    stats.mapq = 0;
    stats.isPaired = true;
    stats.isMapped = false;
    stats.isMateMapped = true;
    return FullRead(std::move(read), std::move(stats));
}

// A spanning read over the whole repeat carrying `repeat` (any length) between the two 20 bp flanks. A repeat
// longer or shorter than the reference is represented as a whole-motif insertion or deletion at S.
FullReadPair makeSpanningRead(const string& name, const string& repeat)
{
    const int lengthChange = static_cast<int>(repeat.size()) - 30;
    vector<uint32_t> cigar;
    if (lengthChange > 0)
    {
        cigar = { cigarOp(20, BAM_CMATCH), cigarOp(lengthChange, BAM_CINS), cigarOp(50, BAM_CMATCH) };
    }
    else if (lengthChange < 0)
    {
        cigar = { cigarOp(20, BAM_CMATCH), cigarOp(-lengthChange, BAM_CDEL), cigarOp(50 + lengthChange, BAM_CMATCH) };
    }
    else
    {
        cigar = { cigarOp(70, BAM_CMATCH) };
    }
    FullReadPair pair;
    pair.firstMate = makeMappedRead(name, MateNumber::kFirstMate, kLeftFlank + repeat + kRightFlank, 80, cigar);
    return pair;
}

struct ReadSet
{
    vector<FullReadPair> pairs;

    vector<const FullReadPair*> pointers() const
    {
        vector<const FullReadPair*> result;
        for (const FullReadPair& pair : pairs)
        {
            result.push_back(&pair);
        }
        return result;
    }
};

// CAG x 10 with a CAA as the third unit.
const string kInterruptedRepeat = "CAGCAGCAA" + repeatMotif("CAG", 7);

} // namespace

TEST(MotifCompositionEligibility, MotifLengthBetweenTwoAndAThirdOfTheReadLength)
{
    EXPECT_FALSE(isEligibleForMotifComposition(1, 150));
    EXPECT_TRUE(isEligibleForMotifComposition(2, 150));
    EXPECT_TRUE(isEligibleForMotifComposition(50, 150));
    EXPECT_FALSE(isEligibleForMotifComposition(51, 150));
}

TEST(MotifCompositionPeriod, RepeatOfTheMotifLengthPasses)
{
    EXPECT_DOUBLE_EQ(1.0, computePeriodScore(repeatMotif("CAG", 10), 3));
    EXPECT_TRUE(passesPeriodTest(repeatMotif("CAG", 10), 3, 0.75));
    // Low-quality (lower case) bases compare case-insensitively.
    EXPECT_TRUE(passesPeriodTest("CAGcagCAGCAGCAGCAG", 3, 0.75));
}

TEST(MotifCompositionPeriod, HomopolymersAndShorterPeriodsAreRejected)
{
    // A homopolymer scores 1.0 at lag 3 but also at lag 1.
    EXPECT_FALSE(passesPeriodTest(string(30, 'G'), 3, 0.75));
    // (CA)n scores 1.0 at lag 4 but also at lag 2.
    EXPECT_FALSE(passesPeriodTest(repeatMotif("CA", 15), 4, 0.75));
    EXPECT_TRUE(passesPeriodTest(repeatMotif("CA", 15), 2, 0.75));
    // Sequence without the period.
    EXPECT_FALSE(passesPeriodTest("ATCGATTGCATGCAATGCCGTTAGGCTAACG", 3, 0.75));
    // N never matches.
    EXPECT_DOUBLE_EQ(0.0, computePeriodScore(string(12, 'N'), 3));
}

TEST(MotifCompositionFrame, OffsetOfTheBestTiling)
{
    EXPECT_EQ(0, computeFrameOffset(repeatMotif("CAG", 5), "CAG"));
    EXPECT_EQ(1, computeFrameOffset("A" + repeatMotif("CAG", 5), "CAG"));
    EXPECT_EQ(2, computeFrameOffset("AG" + repeatMotif("CAG", 5), "CAG"));
    // IUPAC catalog motif.
    EXPECT_EQ(0, computeFrameOffset(repeatMotif("AAGGG", 4), "AARRG"));
}

TEST(MotifCompositionSplitting, WorkedExampleFromThePlan)
{
    // A CAA interruption and a 1 bp deletion: CAG CAG CAA CAG | gap | CAG CAG.
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAACAGCACAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(7u, sequenceSubstrings.size());
    const vector<int> starts = { 0, 3, 6, 9, 12, 14, 17 };
    const vector<SequenceSubstringType> types = { SequenceSubstringType::kKnownMotif,
                                                  SequenceSubstringType::kKnownMotif,
                                                  SequenceSubstringType::kNewMotifCandidate,
                                                  SequenceSubstringType::kKnownMotif,
                                                  SequenceSubstringType::kGap,
                                                  SequenceSubstringType::kKnownMotif,
                                                  SequenceSubstringType::kKnownMotif };
    for (size_t index = 0; index != sequenceSubstrings.size(); ++index)
    {
        EXPECT_EQ(starts[index], sequenceSubstrings[index].offsetWithinRepeatTract) << index;
        EXPECT_EQ(types[index], sequenceSubstrings[index].type) << index;
    }
    EXPECT_EQ(2, sequenceSubstrings[4].length);
}

TEST(MotifCompositionSplitting, CandidateNextToAGapIsRemoved)
{
    // "CAA" followed by a 1-base shift is a substring cut from a longer unit, not a real motif.
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAAGCAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(6, sequenceSubstrings[2].offsetWithinRepeatTract);
    EXPECT_EQ(4, sequenceSubstrings[2].length);
}

TEST(MotifCompositionSplitting, ConsecutiveCandidatesBetweenKnownMotifsAreKept)
{
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAACAACAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSplitting, CandidateAtAReadEndNeedsAMatchingPartialUnit)
{
    // Neither end is at the repeat's edge. The partial unit "AG" before CAA matches the end of CAG, so the read start is trusted.
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("AGCAACAGCAG", "CAG", { "CAG" }, 2, false, false, string());
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);

    // "TG" does not match the end of CAG, so CAA could be a frame-shifted substring and is removed.
    sequenceSubstrings = splitIntoMotifs("TGCAACAGCAG", "CAG", { "CAG" }, 2, false, false, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[0].type);
    EXPECT_EQ(5, sequenceSubstrings[0].length);

    // No partial unit at all at a start that is not at the repeat's edge.
    sequenceSubstrings = splitIntoMotifs("CAACAGCAG", "CAG", { "CAG" }, 0, false, false, string());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[0].type);
}

TEST(MotifCompositionSplitting, EndAtRepeatEdgeMustCarryTheReferencePartialUnit)
{
    // Reference repeat CAG x 3 + "CA": the partial last unit is "CA".
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAACA", "CAG", { "CAG" }, 0, true, true, string("CA"));
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);

    // A 1-base leftover instead of 2 shows a frame shift, so the last substring is not trusted.
    sequenceSubstrings = splitIntoMotifs("CAGCAGCAAC", "CAG", { "CAG" }, 0, true, true, string("CA"));
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSplitting, HomopolymerSequenceSubstringsAreNeverCandidates)
{
    // GGG at a CAG locus is a homopolymer substring.
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGGGGCAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    // ATAT at an ATCT locus repeats a 2 bp unit but is not a homopolymer, so it can be a new motif.
    sequenceSubstrings = splitIntoMotifs("ATCTATATATCT", "ATCT", { "ATCT" }, 0, true, true, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
}

TEST(MotifCompositionSplitting, IupacCatalogMotifMatchesAreRecordedAsTheirOwnBases)
{
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("AAAAGAAGGGAAAAG", "AARRG", { "AAAAG" }, 0, true, true, boost::none);
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kKnownMotif, sequenceSubstrings[0].type);
    EXPECT_EQ(SequenceSubstringType::kCatalogMotifMatch, sequenceSubstrings[1].type);
    EXPECT_EQ(SequenceSubstringType::kKnownMotif, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSoftClip, KeepsRepeatAndDropsWhatFollowsIt)
{
    // Aligned repeat "CAGCAGCAG", then a clip of 9 more repeat bases followed by 15 non-repeat bases.
    const string read = "CAGCAGCAG" + string("CAGCAGCAG") + "TTGACCTAGGATCCA";
    const int clipStart = 9;
    const int kept = computeSoftClipBasesToKeep(read, 0, clipStart, read.size(), 0, true, 3, 0.75);
    EXPECT_GE(kept, 9);
    EXPECT_LE(kept, 12);

    // All repeat: the whole clip is kept.
    const string allRepeat = repeatMotif("CAG", 8);
    EXPECT_EQ(15, computeSoftClipBasesToKeep(allRepeat, 0, 9, allRepeat.size(), 0, true, 3, 0.75));

    // A clip on the left, compared toward the aligned part on its right.
    const string leftClipped = "TTGACCTAGGATCCA" + repeatMotif("CAG", 6);
    const int keptOnLeft = computeSoftClipBasesToKeep(leftClipped, 0, 0, 24, leftClipped.size(), false, 3, 0.75);
    EXPECT_GE(keptOnLeft, 9);
    EXPECT_LE(keptOnLeft, 12);
}

TEST(MotifCompositionPoisson, UpperTail)
{
    EXPECT_DOUBLE_EQ(1.0, computePoissonUpperTail(2.0, 0));
    EXPECT_DOUBLE_EQ(0.0, computePoissonUpperTail(0.0, 1));
    // Values from scipy.stats.poisson.sf; the first two are the plan's CAG x 20 example.
    EXPECT_NEAR(3.2246e-5, computePoissonUpperTail(0.58, 6), 1e-8);
    EXPECT_NEAR(2.6429e-6, computePoissonUpperTail(0.58, 7), 1e-9);
    EXPECT_NEAR(1.10249e-3, computePoissonUpperTail(3.0, 10), 1e-7);
}

TEST(MotifComposition, ReferenceReadsGiveOneMotif)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("read" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(100, 10), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(90, 10), composition->locus.motifPairs.at({ 1, 1 }));
    EXPECT_FALSE(composition->hasAlleleBlocks);

    EXPECT_FALSE(computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}

TEST(MotifComposition, InterruptionSeenInSeveralReadsIsAccepted)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    for (int index = 0; index != 6; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("alt" + std::to_string(index), kInterruptedRepeat));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(154, 16), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifs.at(2));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifPairs.at({ 1, 2 }));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifPairs.at({ 2, 1 }));
    EXPECT_EQ(std::make_pair(90 + 6 * 7, 16), composition->locus.motifPairs.at({ 1, 1 }));

    EXPECT_TRUE(computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}

TEST(MotifComposition, SingleReadWithAnInterruptionIsTreatedAsAnError)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt", kInterruptedRepeat));
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    // The rejected substring becomes a gap: it is not counted as CAG, and it breaks the pairs on both sides.
    EXPECT_EQ(std::make_pair(109, 11), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(90 + 7, 11), composition->locus.motifPairs.at({ 1, 1 }));
}

TEST(MotifComposition, InterruptionWithALowQualityDistinguishingBaseIsNotCounted)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    // The base that tells CAA from CAG is low quality (lower case) in every read.
    const string lowQualityInterruption = "CAGCAGCAa" + repeatMotif("CAG", 7);
    for (int index = 0; index != 6; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("alt" + std::to_string(index), lowQualityInterruption));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
}

TEST(MotifComposition, HeterozygousAllelesGetTheirOwnCounts)
{
    ReadSet reads;
    for (int index = 0; index != 5; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("short" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    const string longAllele = "CAGCAGCAA" + repeatMotif("CAG", 11);
    for (int index = 0; index != 5; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("long" + std::to_string(index), longAllele));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 14 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    ASSERT_TRUE(composition->hasAlleleBlocks);
    EXPECT_EQ(std::make_pair(50, 5), composition->allele1.motifs.at(1));
    EXPECT_EQ(0u, composition->allele1.motifs.count(2));
    EXPECT_EQ(std::make_pair(45, 5), composition->allele1.motifPairs.at({ 1, 1 }));
    EXPECT_EQ(std::make_pair(65, 5), composition->allele2.motifs.at(1));
    EXPECT_EQ(std::make_pair(5, 5), composition->allele2.motifs.at(2));
    EXPECT_EQ(std::make_pair(55, 5), composition->allele2.motifPairs.at({ 1, 1 }));
    EXPECT_EQ(std::make_pair(5, 5), composition->allele2.motifPairs.at({ 1, 2 }));

    // Alleles one unit apart: stutter reads cannot be told apart, so only locus totals are written.
    const auto closeAlleles = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 13, 14 }), false);
    ASSERT_TRUE(closeAlleles);
    EXPECT_FALSE(closeAlleles->hasAlleleBlocks);
}

TEST(MotifComposition, InrepeatReadTakesItsOrientationFromTheAnchoredMate)
{
    ReadSet reads;
    for (int index = 0; index != 4; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    // A forward mate anchored in the left flank; its unmapped mate was sequenced from the other strand, so it is
    // stored as the reverse complement of the repeat.
    for (int index = 0; index != 3; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "irr" + std::to_string(index), MateNumber::kFirstMate, string(30, 'T'), 40, { cigarOp(30, BAM_CMATCH) });
        pair.secondMate = makeUnmappedRead(
            "irr" + std::to_string(index), MateNumber::kSecondMate,
            graphtools::reverseComplement(repeatMotif("CAG", 20)), 40);
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(4 * 10 + 3 * 20, 7), composition->locus.motifs.at(1));

    // Without an anchored mate the read is not used.
    for (FullReadPair& pair : reads.pairs)
    {
        if (pair.secondMate)
        {
            pair.firstMate->s.pos = 2000;
        }
    }
    const auto unanchored = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(unanchored);
    EXPECT_EQ(std::make_pair(4 * 10, 4), unanchored->locus.motifs.at(1));
}

TEST(MotifComposition, NoCallInAFlippedInrepeatReadIsNotHighQuality)
{
    // In-repeat reads carry CAN, where N is a low-quality no-call ('n'). Reverse-complementing the read to the
    // reference orientation turns 'n' into 'N', which must still not count as a high-quality distinguishing base.
    ReadSet reads;
    for (int index = 0; index != 4; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    const string referenceOrientation = repeatMotif("CAG", 3) + "CAN" + repeatMotif("CAG", 16);
    string storedBases = graphtools::reverseComplement(referenceOrientation);
    std::replace(storedBases.begin(), storedBases.end(), 'N', 'n');
    for (int index = 0; index != 6; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "irr" + std::to_string(index), MateNumber::kFirstMate, string(30, 'T'), 40, { cigarOp(30, BAM_CMATCH) });
        pair.secondMate = makeUnmappedRead("irr" + std::to_string(index), MateNumber::kSecondMate, storedBases, 40);
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
}

TEST(MotifComposition, ReadConfidentlyPlacedNearTheRepeatIsNotAnInrepeatRead)
{
    // A repeat-like read that BWA placed with high MAPQ 300 bp to the right of the repeat, whose mate is anchored
    // facing the repeat, is flank sequence (for example a neighbouring repeat), not an in-repeat read.
    ReadSet reads;
    for (int index = 0; index != 4; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    FullReadPair pair;
    pair.firstMate = makeMappedRead("near", MateNumber::kFirstMate, string(30, 'T'), 40, { cigarOp(30, BAM_CMATCH) });
    pair.secondMate = makeMappedRead(
        "near", MateNumber::kSecondMate, repeatMotif("CAG", 20), 430, { cigarOp(60, BAM_CMATCH) }, true, 60);
    reads.pairs.push_back(std::move(pair));
    auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(40, 4), composition->locus.motifs.at(1));

    // With low MAPQ, BWA could not place it, so it is used as an in-repeat read.
    reads.pairs.back().secondMate->s.mapq = 0;
    composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(40 + 20, 5), composition->locus.motifs.at(1));
}

TEST(MotifComposition, NewUnitInsertedAtTheRepeatEdgeIsNotTrusted)
{
    // Every read carries CAG x 10 followed by AAG, which the aligner wrote as a 3 bp insertion exactly at E. Where
    // the aligner put an edge insertion does not show whether it belongs to the repeat, so AAG is not reported.
    ReadSet insertedAtEdge;
    for (int index = 0; index != 8; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "ins" + std::to_string(index), MateNumber::kFirstMate,
            kLeftFlank + repeatMotif("CAG", 10) + "AAG" + kRightFlank, 80,
            { cigarOp(50, BAM_CMATCH), cigarOp(3, BAM_CINS), cigarOp(20, BAM_CMATCH) });
        insertedAtEdge.pairs.push_back(std::move(pair));
    }
    auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, insertedAtEdge.pointers(), RepeatGenotype(3, { 11, 11 }),
        false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    // The same AAG as the repeat's own last unit (a substitution, no insertion) is at an edge the alignment placed,
    // so it is reported.
    ReadSet substitutedLastUnit;
    for (int index = 0; index != 8; ++index)
    {
        substitutedLastUnit.pairs.push_back(
            makeSpanningRead("sub" + std::to_string(index), repeatMotif("CAG", 9) + "AAG"));
    }
    composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, substitutedLastUnit.pointers(),
        RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "AAG" }), composition->motifs);
}

TEST(MotifComposition, FlankBasesAlignedOntoTheRepeatAreCutAtTheFlankAnchor)
{
    // A right flank that starts out resembling the repeat. Reads from a CAG x 6 allele, aligned without the
    // deletion, put 12 bases of that flank on the reference repeat. With the reference flanks known, the right
    // flank's anchor marks where the repeat ends.
    const string repeatLikeRightFlank = "CAACAGCAT" + kRightFlank;
    ReadSet reads;
    for (int index = 0; index != 8; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "short" + std::to_string(index), MateNumber::kFirstMate,
            kLeftFlank + repeatMotif("CAG", 6) + repeatLikeRightFlank.substr(0, 20), 80, { cigarOp(58, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    MotifCompositionLocus locus = makeCagLocus();
    locus.leftFlankSequence = kLeftFlank;
    locus.rightFlankSequence = repeatLikeRightFlank;
    auto composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 6, 6 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(48, 8), composition->locus.motifs.at(1));

    // Without the flanks, the flank bases are split into motif-sized substrings and reported as new motifs.
    composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 6, 6 }), false);
    ASSERT_TRUE(composition);
    EXPECT_GT(composition->motifs.size(), 1u);
}

TEST(MotifComposition, SoftClippedRepeatExtendsAFlankingRead)
{
    ReadSet reads;
    // Aligned through the left flank and 15 bases of repeat, then 15 more repeat bases soft-clipped.
    FullReadPair pair;
    pair.firstMate = makeMappedRead(
        "clipped", MateNumber::kFirstMate, kLeftFlank + repeatMotif("CAG", 10), 80,
        { cigarOp(35, BAM_CMATCH), cigarOp(15, BAM_CSOFT_CLIP) });
    reads.pairs.push_back(std::move(pair));
    const auto composition = computeMotifComposition(
        makeCagLocus(), kReadLength, kRegionExtensionLength, reads.pointers(), boost::none, false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(std::make_pair(10, 1), composition->locus.motifs.at(1));
}

TEST(MotifComposition, NoReadsGiveAnEmptyComposition)
{
    const auto composition
        = computeMotifComposition(makeCagLocus(), kReadLength, kRegionExtensionLength, { }, boost::none, false);
    ASSERT_TRUE(composition);
    EXPECT_TRUE(composition->motifs.empty());
    EXPECT_TRUE(composition->locus.motifs.empty());
}

TEST(MotifCompositionCatalogKnownMotifs, UnusableEntriesAreSkipped)
{
    EXPECT_EQ(
        vector<string>({ "CAA", "CAG", "CCG" }),
        selectCatalogKnownMotifs({ "caa", "CAG", "CAGG", "CNG", "CAA", "CcG", "" }, 3, "locus"));
    EXPECT_TRUE(selectCatalogKnownMotifs({}, 3, "locus").empty());
    // AAC and ACA are CAA written from other starting bases, so they repeat it.
    EXPECT_EQ(vector<string>({ "AAC", "GCA" }), selectCatalogKnownMotifs({ "AAC", "GCA", "ACA", "CAA" }, 3, "locus"));
}

TEST(MotifComposition, CatalogKnownMotifInTwoReadPairsIsCounted)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt0", kInterruptedRepeat));
    MotifCompositionLocus locus = makeCagLocus();
    // CAT is never seen, so it gets no ID.
    locus.catalogKnownMotifs = { "CAG", "CAA", "CAT" };

    // One read pair is not enough, even for a listed motif: it could be a sequencing error.
    auto composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    // Two are: without the catalog list, two CAA read pairs among ten CAG ones still fail the error test.
    reads.pairs.push_back(makeSpanningRead("alt1", kInterruptedRepeat));
    locus.catalogKnownMotifs.clear();
    composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    locus.catalogKnownMotifs = { "CAG", "CAA", "CAT" };
    composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(118, 12), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifs.at(2));
    EXPECT_EQ(std::make_pair(90 + 2 * 7, 12), composition->locus.motifPairs.at({ 1, 1 }));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifPairs.at({ 1, 2 }));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifPairs.at({ 2, 1 }));

    // CAA is not in the reference repeat, so it counts as a non-reference motif.
    EXPECT_TRUE(computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}

TEST(MotifComposition, CatalogKnownMotifIsMatchedInTheRotationTheReadsShow)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt0", kInterruptedRepeat));
    reads.pairs.push_back(makeSpanningRead("alt1", kInterruptedRepeat));
    MotifCompositionLocus locus = makeCagLocus();
    // ACA is CAA written from another starting base; the reads are cut in line with CAG, so they show CAA.
    locus.catalogKnownMotifs = { "ACA" };
    const auto composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifs.at(2));
}

TEST(MotifComposition, UnseenCatalogKnownMotifIsNotANonReferenceMotif)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    MotifCompositionLocus locus = makeCagLocus();
    locus.catalogKnownMotifs = { "CAA" };
    const auto composition = computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    EXPECT_FALSE(computeMotifComposition(
        locus, kReadLength, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}
