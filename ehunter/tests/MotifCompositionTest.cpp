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

// A CAG x 10 locus at [100, 130) on contig 0, with non-repetitive flanks.
const string kLeftFlank = "ATCGATTGCATGCAATGCCG"; // reference [80, 100)
const string kRightFlank = "TTAGGCTAACGTTGACCTAG"; // reference [130, 150)
const int kRegionExtensionLength = 1000;

MotifCompositionLocus makeCagLocus()
{
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.locusStart = 100;
    locus.locusEnd = 130;
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

// A spanning read over the whole locus carrying `repeat` (any length) between the two 20 bp flanks. A repeat
// longer or shorter than the reference repeat sequence is represented as a whole-motif insertion or deletion at S.
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

// CAG x 10 with a CAA as the third repeat unit.
const string kInterruptedRepeatSequence = "CAGCAGCAA" + repeatMotif("CAG", 7);

} // namespace

TEST(MotifCompositionEligibility, MotifSizeBetweenTwoAndAThirdOfTheReadLength)
{
    EXPECT_FALSE(isEligibleForMotifComposition(1, 150));
    EXPECT_TRUE(isEligibleForMotifComposition(2, 150));
    EXPECT_TRUE(isEligibleForMotifComposition(50, 150));
    EXPECT_FALSE(isEligibleForMotifComposition(51, 150));
}

TEST(MotifCompositionPeriod, RepeatSequenceOfTheMotifSizePasses)
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
    EXPECT_EQ(0, computeReferenceRepeatSequenceFrame(repeatMotif("CAG", 5), "CAG"));
    EXPECT_EQ(1, computeReferenceRepeatSequenceFrame("A" + repeatMotif("CAG", 5), "CAG"));
    EXPECT_EQ(2, computeReferenceRepeatSequenceFrame("AG" + repeatMotif("CAG", 5), "CAG"));
    // IUPAC catalog motif.
    EXPECT_EQ(0, computeReferenceRepeatSequenceFrame(repeatMotif("AAGGG", 4), "AARRG"));
}

TEST(MotifCompositionSplitting, InterruptionAndDeletion)
{
    // A CAA interruption and a 1 bp deletion: CAG CAG CAA CAG | gap | CAG CAG.
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAACAGCACAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(7u, sequenceSubstrings.size());
    const vector<int> starts = { 0, 3, 6, 9, 12, 14, 17 };
    const vector<SequenceSubstringType> types = { SequenceSubstringType::kAcceptedMotif,
                                                  SequenceSubstringType::kAcceptedMotif,
                                                  SequenceSubstringType::kNewMotifCandidate,
                                                  SequenceSubstringType::kAcceptedMotif,
                                                  SequenceSubstringType::kGap,
                                                  SequenceSubstringType::kAcceptedMotif,
                                                  SequenceSubstringType::kAcceptedMotif };
    for (size_t index = 0; index != sequenceSubstrings.size(); ++index)
    {
        EXPECT_EQ(starts[index], sequenceSubstrings[index].offsetWithinRepeatTract) << index;
        EXPECT_EQ(types[index], sequenceSubstrings[index].type) << index;
    }
    EXPECT_EQ(2, sequenceSubstrings[4].length);
}

TEST(MotifCompositionSplitting, MotifCandidateNextToAGapIsRemoved)
{
    // "CAA" followed by a 1-base shift is a repeat unit cut out of frame, not a real motif.
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAAGCAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(6, sequenceSubstrings[2].offsetWithinRepeatTract);
    EXPECT_EQ(4, sequenceSubstrings[2].length);
}

TEST(MotifCompositionSplitting, ConsecutiveMotifCandidatesBetweenAcceptedMotifsAreKept)
{
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAACAACAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSplitting, MotifCandidateAtAReadEndNeedsAMatchingPartialRepeatUnit)
{
    // Neither end is at the edge of the repeat sequence. The partial repeat unit "AG" before CAA matches the end of
    // CAG, so the read start is trusted.
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("AGCAACAGCAG", "CAG", { "CAG" }, 2, false, false, string());
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);

    // "TG" does not match the end of CAG, so CAA could be a frame-shifted repeat unit and is removed.
    sequenceSubstrings = splitIntoMotifs("TGCAACAGCAG", "CAG", { "CAG" }, 2, false, false, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[0].type);
    EXPECT_EQ(5, sequenceSubstrings[0].length);

    // No partial repeat unit at all at a start that is not at the edge of the repeat sequence.
    sequenceSubstrings = splitIntoMotifs("CAACAGCAG", "CAG", { "CAG" }, 0, false, false, string());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[0].type);
}

TEST(MotifCompositionSplitting, PartialRepeatUnitMayHaveOneMismatchOnlyFromSixBases)
{
    // A 5-base partial repeat unit must match the end of the motif exactly: CGTTG passes, CGTTA does not.
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CGTTGACCTTGACGTTG", "ACGTTG", { "ACGTTG" }, 5, false, false, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);

    sequenceSubstrings = splitIntoMotifs("CGTTAACCTTGACGTTG", "ACGTTG", { "ACGTTG" }, 5, false, false, string());
    ASSERT_EQ(2u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[0].type);
    EXPECT_EQ(11, sequenceSubstrings[0].length);

    // A 6-base partial repeat unit may have one mismatch: GTTGCT differs from the motif's end GTTGCA at one base.
    sequenceSubstrings
        = splitIntoMotifs("GTTGCTACCTTGCAACGTTGCA", "ACGTTGCA", { "ACGTTGCA" }, 6, false, false, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
}

TEST(MotifCompositionSplitting, EndAtEdgeOfRepeatSequenceMustCarryTheReferencePartialRepeatUnit)
{
    // Reference repeat sequence CAG x 3 + "CA": the partial last repeat unit is "CA".
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGCAACA", "CAG", { "CAG" }, 0, true, true, string("CA"));
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);

    // A 1-base leftover instead of 2 shows a frame shift, so the last repeat unit is not trusted.
    sequenceSubstrings = splitIntoMotifs("CAGCAGCAAC", "CAG", { "CAG" }, 0, true, true, string("CA"));
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSplitting, HomopolymerRepeatUnitsAreNeverMotifCandidates)
{
    // GGG at a CAG locus is a homopolymer repeat unit.
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("CAGCAGGGGCAGCAG", "CAG", { "CAG" }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    // ATAT at an ATCT locus repeats a 2 bp repeat unit but is not a homopolymer, so it can be a new motif.
    sequenceSubstrings = splitIntoMotifs("ATCTATATATCT", "ATCT", { "ATCT" }, 0, true, true, string());
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
}

TEST(MotifCompositionSplitting, IupacCatalogMotifMatchesAreRecordedAsTheirOwnBases)
{
    const vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs("AAAAGAAGGGAAAAG", "AARRG", { "AAAAG" }, 0, true, true, boost::none);
    ASSERT_EQ(3u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kAcceptedMotif, sequenceSubstrings[0].type);
    EXPECT_EQ(SequenceSubstringType::kCatalogMotifMatch, sequenceSubstrings[1].type);
    EXPECT_EQ(SequenceSubstringType::kAcceptedMotif, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSplitting, AcceptedMotifOfAnotherLengthIsRejected)
{
    EXPECT_THROW(
        splitIntoMotifs("CAGCAGCAGCAG", "CAG", { "CAG", "CAGCAG" }, 0, true, true, string()), std::logic_error);
    EXPECT_THROW(splitIntoMotifs("CAGCAGCAGCAG", "CAG", { "CA" }, 0, true, true, string()), std::logic_error);
}

TEST(MotifCompositionSplitting, UnanchoredMotifCandidateThatPassesTheInFrameTestIsKept)
{
    // ACGTTGCAAT is 1 base from the 10 bp motif in frame and at least 7 from each of its rotations. It is followed
    // by a 2-base insertion, so its right side is not anchored.
    const string motif = "ACGTTGCAAG";
    const string tract = motif + motif + "ACGTTGCAAT" + "GG" + motif + motif;
    vector<SequenceSubstring> sequenceSubstrings = splitIntoMotifs(tract, motif, { motif }, 0, true, true, string());
    ASSERT_EQ(6u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[3].type);

    // The closest motif in frame can be an accepted motif rather than the catalog motif.
    sequenceSubstrings = splitIntoMotifs(tract, "TTGACCATGA", { motif }, 0, true, true, string());
    ASSERT_EQ(6u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);

    // Neither side needs to be anchored: here both are untrusted read ends.
    sequenceSubstrings = splitIntoMotifs("ACGTTGCAAT" + string("TT"), motif, { motif }, 0, false, false, string());
    ASSERT_EQ(2u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[0].type);

    // A motif candidate kept this way anchors its neighbors: ACCATGCAAG, 2 bases from the motif, fails the in-frame
    // test but touches an accepted motif on its left and the kept ACGTTGCAAT on its right.
    sequenceSubstrings = splitIntoMotifs(
        motif + motif + "ACCATGCAAG" + "ACGTTGCAAT" + "GG" + motif + motif, motif, { motif }, 0, true, true, string());
    ASSERT_EQ(7u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[3].type);
}

TEST(MotifCompositionSplitting, UnanchoredMotifCandidateThatFailsTheInFrameTestIsRemoved)
{
    // ACGTTGCATT is 2 bases from the 10 bp motif, and the in-frame test allows 1 at 10 bp.
    const string motif = "ACGTTGCAAG";
    vector<SequenceSubstring> sequenceSubstrings = splitIntoMotifs(
        motif + motif + "ACGTTGCATT" + "GG" + motif + motif, motif, { motif }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(12, sequenceSubstrings[2].length);

    // AGAGAGAGAG is 1 base from AGAGAGAGAC in frame, but also from its rotation AGAGAGACAG.
    const string dinucleotideLikeMotif = "AGAGAGAGAC";
    sequenceSubstrings = splitIntoMotifs(
        repeatMotif(dinucleotideLikeMotif, 2) + "AGAGAGAGAG" + "TT" + repeatMotif(dinucleotideLikeMotif, 2),
        dinucleotideLikeMotif, { dinucleotideLikeMotif }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(12, sequenceSubstrings[2].length);
}

TEST(MotifCompositionSplitting, InFrameTestNeverKeepsAMotifCandidateBelowTenBases)
{
    // ACGTTGCAT is 1 base from the 9 bp motif in frame, but below 10 bp the in-frame test allows no mismatch.
    const string motif = "ACGTTGCAG";
    const vector<SequenceSubstring> sequenceSubstrings = splitIntoMotifs(
        motif + motif + "ACGTTGCAT" + "GG" + motif + motif, motif, { motif }, 0, true, true, string());
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(11, sequenceSubstrings[2].length);
}

TEST(MotifCompositionSplitting, CatalogKnownMotifsReplaceTheInFrameTest)
{
    // With a KnownMotifs list, a motif candidate the anchoring rule would remove is kept if it is a rotation of a
    // listed motif, even below 10 bp, where the in-frame test never keeps one.
    const string motif = "ACGTTGCAG";
    const string tract = motif + motif + "ACGTTGCAT" + "GG" + motif + motif;
    vector<SequenceSubstring> sequenceSubstrings
        = splitIntoMotifs(tract, motif, { motif }, 0, true, true, string(), { motif, "ACGTTGCAT" });
    ASSERT_EQ(6u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[3].type);

    // The list may write the motif in any rotation: TACGTTGCA is ACGTTGCAT.
    sequenceSubstrings = splitIntoMotifs(tract, motif, { motif }, 0, true, true, string(), { "TACGTTGCA" });
    ASSERT_EQ(6u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);

    // A motif candidate that is not listed is removed, even one the in-frame test would keep: ACGTTGCAAT is 1 base
    // from the 10 bp motif in frame.
    const string longerMotif = "ACGTTGCAAG";
    sequenceSubstrings = splitIntoMotifs(
        longerMotif + longerMotif + "ACGTTGCAAT" + "GG" + longerMotif + longerMotif, longerMotif, { longerMotif }, 0,
        true, true, string(), { longerMotif });
    ASSERT_EQ(5u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kGap, sequenceSubstrings[2].type);
    EXPECT_EQ(12, sequenceSubstrings[2].length);

    // Motif candidates anchored on both sides are kept whether or not they are listed.
    sequenceSubstrings = splitIntoMotifs("CAGCAACAACAG", "CAG", { "CAG" }, 0, true, true, string(), { "CAG" });
    ASSERT_EQ(4u, sequenceSubstrings.size());
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[1].type);
    EXPECT_EQ(SequenceSubstringType::kNewMotifCandidate, sequenceSubstrings[2].type);
}

TEST(MotifCompositionSoftClip, KeepsRepeatSequenceAndDropsWhatFollowsIt)
{
    // Aligned repeat "CAGCAGCAG", then a clip of 9 more repeat bases followed by 15 non-repeat bases.
    const string read = "CAGCAGCAG" + string("CAGCAGCAG") + "TTGACCTAGGATCCA";
    const int clipStart = 9;
    const int kept = computeSoftClipBasesToKeep(read, 0, clipStart, read.size(), 0, true, 3, 0.75);
    EXPECT_GE(kept, 9);
    EXPECT_LE(kept, 12);

    // All repeat: the whole clip is kept.
    const string allRepeatSequence = repeatMotif("CAG", 8);
    EXPECT_EQ(15, computeSoftClipBasesToKeep(allRepeatSequence, 0, 9, allRepeatSequence.size(), 0, true, 3, 0.75));

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
    // Values from scipy.stats.poisson.sf.
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(100, 10), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(90, 10), composition->locus.motifPairs.at({ 1, 1 }));
    EXPECT_FALSE(composition->hasAlleleBlocks);

    EXPECT_FALSE(computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
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
        reads.pairs.push_back(makeSpanningRead("alt" + std::to_string(index), kInterruptedRepeatSequence));
    }
    const auto composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(154, 16), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifs.at(2));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifPairs.at({ 1, 2 }));
    EXPECT_EQ(std::make_pair(6, 6), composition->locus.motifPairs.at({ 2, 1 }));
    EXPECT_EQ(std::make_pair(90 + 6 * 7, 16), composition->locus.motifPairs.at({ 1, 1 }));

    EXPECT_TRUE(computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}

TEST(MotifComposition, SingleReadWithAnInterruptionIsTreatedAsAnError)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt", kInterruptedRepeatSequence));
    const auto composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    // The rejected repeat unit becomes a gap: it is not counted as CAG, and it breaks the pairs on both sides.
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 14 }), false);
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

    // Alleles one repeat unit apart: stutter reads cannot be told apart, so only locus totals are written.
    const auto closeAlleles = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 13, 14 }), false);
    ASSERT_TRUE(closeAlleles);
    EXPECT_FALSE(closeAlleles->hasAlleleBlocks);
}

TEST(MotifComposition, InsideReadIsAssignedToAnAlleleByTheRepeatSequenceItKeeps)
{
    // Reads of the short allele (10 repeat units), plus a 48-base read mapped inside the locus whose last 30 bases
    // are soft-clipped flank sequence. Only its 18 aligned bases are repeat sequence, which fit in the short allele
    // as well as in the long one (60 repeat units), so the read is counted for the locus but for neither allele.
    ReadSet reads;
    for (int index = 0; index != 5; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("short" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    FullReadPair pair;
    pair.firstMate = makeMappedRead(
        "inside", MateNumber::kFirstMate, repeatMotif("CAG", 6) + "TTAGGCTAACGTTGACCTAGATCGATTGCA", 106,
        { cigarOp(18, BAM_CMATCH), cigarOp(30, BAM_CSOFT_CLIP) });
    reads.pairs.push_back(std::move(pair));
    const auto composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 60 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(std::make_pair(56, 6), composition->locus.motifs.at(1));
    ASSERT_TRUE(composition->hasAlleleBlocks);
    EXPECT_EQ(std::make_pair(50, 5), composition->allele1.motifs.at(1));
    EXPECT_TRUE(composition->allele2.motifs.empty());
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(unanchored);
    EXPECT_EQ(std::make_pair(4 * 10, 4), unanchored->locus.motifs.at(1));
}

TEST(MotifComposition, InrepeatReadIsCutAtTheInnermostFlankAnchorMatch)
{
    // HG002's chr22:23881827-23881841 (AC x 7) lies in a tandem array of ~120bp units, each holding an AC repeat and
    // an AG repeat, so the locus's flanks recur one unit away, and every read there has MAPQ 0. A read placed in
    // the next unit is an in-repeat read; one that starts in an AG repeat and runs on through the AC repeat, the
    // right flank, the next AG repeat and the next AC repeat matches each flank's anchor twice. It is cut at the
    // matches nearest its middle, so that only the AC units between them are counted and none of the neighbouring
    // AG repeat; cut at the outermost matches it would hold a whole array unit.
    const string leftFlank = "TGAGAGAGACAGAGGCAGAGAGAGAGAGAA";
    const string rightFlank = "AGACACAGAGACAGAGAGATTGAGAGAGAC";
    const string repeat = repeatMotif("AC", 7);
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.locusStart = 23881827;
    locus.locusEnd = 23881841;
    locus.catalogMotif = "AC";
    locus.referenceRepeatSequence = repeat;
    locus.leftFlankSequence = leftFlank;
    locus.rightFlankSequence = rightFlank;
    locus.meanFragmentLength = 300;
    ReadSet reads;
    for (int index = 0; index != 4; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "span" + std::to_string(index), MateNumber::kFirstMate, leftFlank + repeat + rightFlank,
            locus.locusStart - 30, { cigarOp(74, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    const string readThroughTheArrayUnit
        = leftFlank.substr(18) + repeat + rightFlank + "AGAGGCAGAGAGAGAGAA" + repeat + rightFlank.substr(0, 12);
    for (int index = 0; index != 3; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "irr" + std::to_string(index), MateNumber::kFirstMate, string(30, 'T'), locus.locusStart - 60,
            { cigarOp(30, BAM_CMATCH) });
        pair.secondMate = makeUnmappedRead(
            "irr" + std::to_string(index), MateNumber::kSecondMate,
            graphtools::reverseComplement(readThroughTheArrayUnit), locus.locusStart - 60);
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(2, { 7, 7 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "AC" }), composition->motifs);
    EXPECT_EQ(std::make_pair(4 * 7 + 3 * 7, 7), composition->locus.motifs.at(1));
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
}

TEST(MotifComposition, ReadConfidentlyPlacedNearTheLocusIsNotAnInrepeatRead)
{
    // A repeat-like read that BWA placed with high MAPQ 300 bp to the right of the locus is flank sequence (for
    // example a neighbouring repeat), not an in-repeat read, even though its mate is on the forward strand in the
    // left flank, where the pair's other read would lie over the locus.
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
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(40, 4), composition->locus.motifs.at(1));

    // With low MAPQ, BWA could not place it, so it is used as an in-repeat read.
    reads.pairs.back().secondMate->s.mapq = 0;
    composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(40 + 20, 5), composition->locus.motifs.at(1));
}

TEST(MotifComposition, NewRepeatUnitInsertedAtTheEdgeOfRepeatSequenceIsNotTrusted)
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
        makeCagLocus(), kRegionExtensionLength, insertedAtEdge.pointers(), RepeatGenotype(3, { 11, 11 }),
        false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    // The same AAG as the repeat's own last repeat unit (a substitution, no insertion) is at an edge the alignment
    // placed, so it is reported.
    ReadSet substitutedLastRepeatUnit;
    for (int index = 0; index != 8; ++index)
    {
        substitutedLastRepeatUnit.pairs.push_back(
            makeSpanningRead("sub" + std::to_string(index), repeatMotif("CAG", 9) + "AAG"));
    }
    composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, substitutedLastRepeatUnit.pointers(),
        RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "AAG" }), composition->motifs);
}

TEST(MotifComposition, EndAtEdgeOfRepeatSequenceIsTrustedPastAnImpurityInTheReferenceRepeatSequence)
{
    // The 19-base reference repeat sequence CAG CAG CAG C CAG CAG CAG has an extra C. The splitter shifts its frame by
    // one base there, so the sequence ends with no partial repeat unit, and so do reads ending at E. Their end is
    // therefore trusted, and a CAA in their last repeat unit is anchored on both sides.
    const string referenceRepeatSequence = "CAGCAGCAGCCAGCAGCAG";
    MotifCompositionLocus locus = makeCagLocus();
    locus.locusEnd = 119;
    locus.referenceRepeatSequence = referenceRepeatSequence;
    ReadSet reads;
    for (int index = 0; index != 8; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "ref" + std::to_string(index), MateNumber::kFirstMate, kLeftFlank + referenceRepeatSequence + kRightFlank,
            80, { cigarOp(59, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    for (int index = 0; index != 5; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "alt" + std::to_string(index), MateNumber::kFirstMate, kLeftFlank + "CAGCAGCAGCCAGCAGCAA" + kRightFlank, 80,
            { cigarOp(59, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 6, 6 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(5, 5), composition->locus.motifs.at(2));
}

TEST(MotifComposition, FlankBasesAlignedOntoTheLocusAreCutAtTheFlankAnchor)
{
    // A right flank that starts out resembling the repeat. Reads from a CAG x 6 allele, aligned without the
    // deletion, put 12 bases of that flank on the reference repeat sequence. With the reference flanks known, the
    // right flank's anchor marks where the repeat sequence ends.
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
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 6, 6 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(48, 8), composition->locus.motifs.at(1));

    // Without the flanks, the flank bases are split into repeat units and reported as new motifs.
    composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 6, 6 }), false);
    ASSERT_TRUE(composition);
    EXPECT_GT(composition->motifs.size(), 1u);
}

TEST(MotifComposition, RepeatSequenceContainingTheFlankAnchorIsNotCutThere)
{
    // HG002's chr22:16261470-16261537 (TGATTCCATT x 6.7 in the reference). Its 47bp allele contains the right
    // flank's anchor, AATGATGATTCC, 21 bases in, ahead of a TGATTCCATT and a CGATTCCATT unit. Reads of that allele,
    // aligned with the deletion, end their repeat sequence exactly where the flank begins; a cut at the inner anchor
    // match would drop those two units.
    const string leftFlank = "ATGATGATTCCACTCGATTCCATATGATAA";
    const string referenceRepeatSequence = "TGATTCCATTTGATTCCACTCGATGATTCCATTTGATTCCATTCAATGGTTCCATTCGATTCCATTC";
    const string rightFlank = "AATGATGATTCCATTCGAGTTCATTGATTA";
    const string shortAllele = "TGATTCCATTTGATTCCATTCAATGATGATTCCATTCGATTCCATTC";
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.locusStart = 16261470;
    locus.locusEnd = 16261537;
    locus.catalogMotif = "TGATTCCATT";
    locus.referenceRepeatSequence = referenceRepeatSequence;
    locus.leftFlankSequence = leftFlank;
    locus.rightFlankSequence = rightFlank;
    locus.meanFragmentLength = 300;
    ReadSet reads;
    for (int index = 0; index != 8; ++index)
    {
        // 30 flank bases, the allele's first 25 bases, the 20bp deletion, its last 22 bases, 20 flank bases.
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "short" + std::to_string(index), MateNumber::kFirstMate, leftFlank + shortAllele + rightFlank.substr(0, 20),
            locus.locusStart - 30, { cigarOp(55, BAM_CMATCH), cigarOp(20, BAM_CDEL), cigarOp(42, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(10, { 5, 5 }), false);
    ASSERT_TRUE(composition);
    // Three TGATTCCATT units per read and the CGATTCCATT unit. Cut at the inner anchor match, a read kept a single
    // TGATTCCATT and no CGATTCCATT at all. The first unit also tests the frame: the reference repeat sequence's best
    // tiling is at offset 3 (its third and fourth units follow a 3-base shift), and a split that started there would
    // skip the exact unit at the tract's start.
    EXPECT_EQ(vector<string>({ "TGATTCCATT", "CGATTCCATT" }), composition->motifs);
    EXPECT_EQ(std::make_pair(24, 8), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(8, 8), composition->locus.motifs.at(2));
}

TEST(MotifComposition, FlankAnchorFoundOutsideTheTractDoesNotMoveItsEdge)
{
    // HG002's chr22:15864961-15864977 (AG x 8). The left flank ends in TATA repeats, so its anchor CATATGTATATA
    // sits 9 bases before the locus. A read with 4 of those TATA bases deleted still matches the anchor exactly,
    // 12 bases before its repeat sequence; placing the flank's edge 9 bases past the anchor would cut two AG units
    // that the alignment already put in the locus.
    const string leftFlank = "ACACACACACATATGTATATATATATATAT";
    const string rightFlank = "GAAAGTTTGAATTTACCTATATTAAAAGAT";
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.locusStart = 15864961;
    locus.locusEnd = 15864977;
    locus.catalogMotif = "AG";
    locus.referenceRepeatSequence = repeatMotif("AG", 8);
    locus.leftFlankSequence = leftFlank;
    locus.rightFlankSequence = rightFlank;
    locus.meanFragmentLength = 300;
    ReadSet reads;
    for (int index = 0; index != 8; ++index)
    {
        // The flank's first 15 bases, the 4-base deletion, its last 11 bases, the repeat, 20 bases of right flank.
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "deletion" + std::to_string(index), MateNumber::kFirstMate,
            leftFlank.substr(0, 15) + leftFlank.substr(19) + repeatMotif("AG", 8) + rightFlank.substr(0, 20),
            locus.locusStart - 30, { cigarOp(15, BAM_CMATCH), cigarOp(4, BAM_CDEL), cigarOp(47, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(2, { 8, 8 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "AG" }), composition->motifs);
    EXPECT_EQ(std::make_pair(64, 8), composition->locus.motifs.at(1));
}

TEST(MotifComposition, RepeatUnitsNextToARepeatLikeFlankAreNotCutOff)
{
    // A GATC locus (GATC GATC GATC T) whose right flank is ATCT repeats. Every 12-mer of that flank also occurs one
    // period closer to the locus, across the junction with the reference repeat sequence, so the flank gives no
    // anchor. Otherwise ATCTATCTATCT would be the anchor: reads with GATC GATC TATC T come within one mismatch of it
    // inside their repeat sequence, and reads with 8 extra bases at E (GATC GATC GATC TATC TATC T) match it exactly
    // inside their third GATC. All of their GATC repeat units are kept.
    const string referenceRepeatSequence = "GATCGATCGATCT";
    const string rightFlank = repeatMotif("ATCT", 7) + "AT";
    MotifCompositionLocus locus;
    locus.contigIndex = 0;
    locus.locusStart = 100;
    locus.locusEnd = 113;
    locus.catalogMotif = "GATC";
    locus.referenceRepeatSequence = referenceRepeatSequence;
    locus.leftFlankSequence = kLeftFlank;
    locus.rightFlankSequence = rightFlank;
    locus.meanFragmentLength = 300;
    ReadSet reads;
    for (int index = 0; index != 12; ++index)
    {
        const string repeatSequence = index < 6 ? referenceRepeatSequence : "GATCGATCTATCT";
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "read" + std::to_string(index), MateNumber::kFirstMate,
            kLeftFlank + repeatSequence + rightFlank.substr(0, 20), 80, { cigarOp(53, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    for (int index = 0; index != 6; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "insertion" + std::to_string(index), MateNumber::kFirstMate,
            kLeftFlank + "GATCGATCGATCTATCTATCT" + rightFlank.substr(0, 20), 80,
            { cigarOp(33, BAM_CMATCH), cigarOp(8, BAM_CINS), cigarOp(20, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(4, { 3, 3 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "GATC", "TATC" }), composition->motifs);
    EXPECT_EQ(std::make_pair(6 * 3 + 6 * 2 + 6 * 3, 18), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(6 + 6 * 2, 12), composition->locus.motifs.at(2));
}

TEST(MotifComposition, SoftClippedRepeatSequenceExtendsAFlankingRead)
{
    ReadSet reads;
    // Aligned through the left flank and 15 bases of repeat, then 15 more repeat bases soft-clipped.
    FullReadPair pair;
    pair.firstMate = makeMappedRead(
        "clipped", MateNumber::kFirstMate, kLeftFlank + repeatMotif("CAG", 10), 80,
        { cigarOp(35, BAM_CMATCH), cigarOp(15, BAM_CSOFT_CLIP) });
    reads.pairs.push_back(std::move(pair));
    const auto composition = computeMotifComposition(
        makeCagLocus(), kRegionExtensionLength, reads.pointers(), boost::none, false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(std::make_pair(10, 1), composition->locus.motifs.at(1));
}

TEST(MotifComposition, NoReadsGiveAnEmptyComposition)
{
    const auto composition
        = computeMotifComposition(makeCagLocus(), kRegionExtensionLength, { }, boost::none, false);
    ASSERT_TRUE(composition);
    EXPECT_TRUE(composition->motifs.empty());
    EXPECT_TRUE(composition->locus.motifs.empty());
}

TEST(MotifCompositionCatalogKnownMotifs, UnusableEntriesAreSkipped)
{
    EXPECT_EQ(
        vector<string>({ "CAA", "CAG", "CCG" }),
        validateKnownMotifs({ "caa", "CAG", "CAGG", "CNG", "CAA", "CcG", "" }, 3, "locus"));
    EXPECT_TRUE(validateKnownMotifs({}, 3, "locus").empty());
    // AAC and ACA are CAA written from other starting bases, so they repeat it.
    EXPECT_EQ(vector<string>({ "AAC", "GCA" }), validateKnownMotifs({ "AAC", "GCA", "ACA", "CAA" }, 3, "locus"));
}

TEST(MotifComposition, CatalogKnownMotifInTwoReadPairsIsCounted)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt0", kInterruptedRepeatSequence));
    MotifCompositionLocus locus = makeCagLocus();
    // CAT is never seen, so it gets no ID.
    locus.knownMotifs = { "CAG", "CAA", "CAT" };

    // One read pair is not enough, even for a listed motif: it could be a sequencing error.
    auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    // Two are: without the catalog list, two CAA read pairs among ten CAG ones still fail the error test.
    reads.pairs.push_back(makeSpanningRead("alt1", kInterruptedRepeatSequence));
    locus.knownMotifs.clear();
    composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    locus.knownMotifs = { "CAG", "CAA", "CAT" };
    composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(118, 12), composition->locus.motifs.at(1));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifs.at(2));
    EXPECT_EQ(std::make_pair(90 + 2 * 7, 12), composition->locus.motifPairs.at({ 1, 1 }));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifPairs.at({ 1, 2 }));
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifPairs.at({ 2, 1 }));

    // CAA is not in the reference repeat sequence, so it counts as a motif not in the reference repeat sequence.
    EXPECT_TRUE(computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}

TEST(MotifComposition, CatalogKnownMotifIsMatchedInTheRotationTheReadsShow)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    reads.pairs.push_back(makeSpanningRead("alt0", kInterruptedRepeatSequence));
    reads.pairs.push_back(makeSpanningRead("alt1", kInterruptedRepeatSequence));
    MotifCompositionLocus locus = makeCagLocus();
    // ACA is CAA written from another starting base; the reads are cut in line with CAG, so they show CAA.
    locus.knownMotifs = { "ACA" };
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifs.at(2));
}

TEST(MotifComposition, CatalogKnownMotifNextToAnInsertionIsCounted)
{
    // Two read pairs carry CAA as the third repeat unit, followed by a 2 bp insertion, so CAA touches a motif on one
    // side only.
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    for (int index = 0; index != 2; ++index)
    {
        FullReadPair pair;
        pair.firstMate = makeMappedRead(
            "alt" + std::to_string(index), MateNumber::kFirstMate,
            kLeftFlank + "CAGCAGCAA" + "GG" + repeatMotif("CAG", 7) + kRightFlank, 80,
            { cigarOp(29, BAM_CMATCH), cigarOp(2, BAM_CINS), cigarOp(41, BAM_CMATCH) });
        reads.pairs.push_back(std::move(pair));
    }
    MotifCompositionLocus locus = makeCagLocus();

    // With CAA listed, the list keeps it as a motif candidate, and it is trusted like any listed motif seen in two read
    // pairs. Otherwise the anchoring rule would remove it before it could be trusted, since the in-frame test never
    // keeps a motif candidate below 10 bp.
    locus.knownMotifs = { "CAG", "CAA" };
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG", "CAA" }), composition->motifs);
    EXPECT_EQ(std::make_pair(2, 2), composition->locus.motifs.at(2));
}

TEST(MotifComposition, UnseenCatalogKnownMotifIsNotCounted)
{
    ReadSet reads;
    for (int index = 0; index != 10; ++index)
    {
        reads.pairs.push_back(makeSpanningRead("ref" + std::to_string(index), repeatMotif("CAG", 10)));
    }
    MotifCompositionLocus locus = makeCagLocus();
    locus.knownMotifs = { "CAA" };
    const auto composition = computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), false);
    ASSERT_TRUE(composition);
    EXPECT_EQ(vector<string>({ "CAG" }), composition->motifs);

    EXPECT_FALSE(computeMotifComposition(
        locus, kRegionExtensionLength, reads.pointers(), RepeatGenotype(3, { 10, 10 }), true));
}
