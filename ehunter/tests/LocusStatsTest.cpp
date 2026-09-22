//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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

#include "core/LocusStats.hh"

#include "gtest/gtest.h"
#include <htslib/sam.h>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphcore/GraphBuilders.hh"

using namespace ehunter;
using graphtools::decodeGraphAlignment;
using graphtools::GraphAlignment;

TEST(LocusStatsCalculator, NoDataGiven_StatsNotCalculated)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    GraphAlignment alignmentStartingInsideRepeat = decodeGraphAlignment(0, "1[3M]", &graph);
    GraphAlignment alignmentStartingOnRightFlank = decodeGraphAlignment(0, "2[4M]", &graph);

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 0, 0, 0.0), statsCalculator.estimate(Sex::kFemale));
}

TEST(LocusStatsCalculator, TypicalReadLengths_StatsCalculated)
{
    graphtools::Graph graph = graphtools::makeStrGraph("TAATG", "CCG", "CCTTATTA");

    LocusStatsCalculator statsCalculator(ChromType::kAutosome, graph);

    GraphAlignment alignmentStartingOnLeftFlank = decodeGraphAlignment(3, "0[2M]1[2M]", &graph);
    GraphAlignment alignmentStartingInsideRepeat = decodeGraphAlignment(0, "1[3M]", &graph);
    GraphAlignment alignmentStartingOnRightFlank = decodeGraphAlignment(0, "2[3M]", &graph);

    for (int index = 0; index != 29; ++index)
    {
        statsCalculator.recordReadLen(alignmentStartingOnLeftFlank);
        statsCalculator.recordReadLen(alignmentStartingInsideRepeat);
        statsCalculator.recordReadLen(alignmentStartingOnRightFlank);
    }
    statsCalculator.recordReadLen(alignmentStartingOnRightFlank);
    statsCalculator.recordReadLen(alignmentStartingOnRightFlank);

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 3, 0, 18), statsCalculator.estimate(Sex::kFemale));
}

// -- LocusStatsCalculatorFromReadAlignments (the optimized-streaming fast path) ------------------------
//
// Its depth is deliberately the same quantity LocusStatsCalculator reports: coverage of the two reference
// flanks, from reads whose alignment STARTS in a flank. The tests below pin that definition, since the
// emitted Coverage field is only comparable across genotyping modes while it holds.

namespace
{

const int kTestContigIndex = 0;
const int kTestExtensionLength = 100;
// Repeat region [1000, 1030) with 100bp flanks: left flank [900, 1000), right flank [1030, 1130).
const GenomicRegion kTestRepeatRegion(kTestContigIndex, 1000, 1030);

FullRead makeFullRead(const std::string& fragmentId, MateNumber mateNumber, int readLength, int32_t position)
{
    LinearAlignmentStats alignmentStats;
    alignmentStats.chromId = kTestContigIndex;
    alignmentStats.pos = position;
    alignmentStats.isMapped = true;
    return FullRead(Read(ReadId(fragmentId, mateNumber), std::string(readLength, 'A'), false), alignmentStats);
}

FullReadPair makeReadPair(const std::string& fragmentId, int readLength, int32_t readPos, int32_t matePos)
{
    return FullReadPair(
        makeFullRead(fragmentId, MateNumber::kFirstMate, readLength, readPos),
        makeFullRead(fragmentId, MateNumber::kSecondMate, readLength, matePos));
}

} // namespace

TEST(LocusStatsCalculatorFromReadAlignments, NoReadsGiven_StatsNotCalculated)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 0, 0, 0.0), statsCalculator.estimate(Sex::kFemale));
}

TEST(LocusStatsCalculatorFromReadAlignments, ReadsStartingOnFlanks_DepthMatchesFlankCoverage)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    // 10 pairs, one mate starting on each flank, all 50bp. Both mates count, so readCount is 20 and the
    // number of start positions is 100 + 100 - 50, giving 50 * 20 / 150.
    for (int index = 0; index != 10; ++index)
    {
        statsCalculator.inspect(makeReadPair("frag" + std::to_string(index), 50, 950, 1050));
    }

    const LocusStats stats = statsCalculator.estimate(Sex::kFemale);
    EXPECT_EQ(50, stats.meanReadLength());
    EXPECT_NEAR(50.0 * 20 / 150, stats.depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, ReadsStartingInsideRepeatOrOutsideFlanks_NotCounted)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    // Starting inside the repeat region, and starting beyond either flank: none of these are the linear
    // stand-in for "the alignment path's first node is a flank node".
    statsCalculator.inspect(makeReadPair("insideRepeat", 50, 1010, 1020));
    statsCalculator.inspect(makeReadPair("beforeLeftFlank", 50, 899, 800));
    statsCalculator.inspect(makeReadPair("afterRightFlank", 50, 1130, 1200));

    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 0, 0, 0.0), statsCalculator.estimate(Sex::kFemale));
}

TEST(LocusStatsCalculatorFromReadAlignments, FlankBoundaries_AreHalfOpen)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    // The first base of each flank is inside it; the first base past each flank is not. The read starting
    // at 1129 is inside the right flank but runs to 1179, past the locus window's end at 1130, so it is
    // not counted (see ReadsExtendingPastTheLocusWindow_NotCounted below).
    statsCalculator.inspect(makeReadPair("leftEdges", 50, 900, 999));
    statsCalculator.inspect(makeReadPair("rightEdges", 50, 1030, 1129));

    const LocusStats stats = statsCalculator.estimate(Sex::kFemale);
    EXPECT_NEAR(50.0 * 3 / 150, stats.depth(), 1e-9);

    // The repeat's own first base (1000) belongs to the repeat, not to the left flank. Without this the
    // left boundary could be widened to `<=` and every other case here would still pass.
    LocusStatsCalculatorFromReadAlignments repeatStartCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    repeatStartCalculator.inspect(makeReadPair("repeatFirstBase", 50, 1000, 1000));
    ASSERT_EQ(LocusStats(AlleleCount::kTwo, 0, 0, 0.0), repeatStartCalculator.estimate(Sex::kFemale));

    // ...while the base immediately before it is the last base of the left flank.
    LocusStatsCalculatorFromReadAlignments lastLeftFlankBaseCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    lastLeftFlankBaseCalculator.inspect(makeReadPair("lastLeftFlankBase", 50, 999, 999));
    EXPECT_NEAR(50.0 * 2 / 150, lastLeftFlankBaseCalculator.estimate(Sex::kFemale).depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, LeadingSoftClip_ReadPlacedAtItsFirstBase)
{
    // A read from an expanded allele: its first 30 bases are repeat sequence the aligner soft-clipped, so
    // BAM POS is 1040 (in the right flank) but the read's first base sits at 1010, inside the repeat. Placed
    // by POS it would count as right-flank-anchored, and such reads only exist for alleles longer than the
    // reference.
    const auto cigarOp = [](int length, int op) { return (static_cast<uint32_t>(length) << BAM_CIGAR_SHIFT) | op; };

    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    FullReadPair clippedPair = makeReadPair("clippedIntoRepeat", 50, 950, 1040);
    clippedPair.secondMate->s.cigar = { cigarOp(30, BAM_CSOFT_CLIP), cigarOp(20, BAM_CMATCH) };
    statsCalculator.inspect(clippedPair);
    EXPECT_NEAR(50.0 * 1 / 150, statsCalculator.estimate(Sex::kFemale).depth(), 1e-9);

    // A hard clip before the soft clip does not change that, and a read whose soft-clipped bases still
    // start in the flank stays counted.
    LocusStatsCalculatorFromReadAlignments flankClipCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    FullReadPair flankClippedPair = makeReadPair("clippedInFlank", 50, 950, 1060);
    flankClippedPair.secondMate->s.cigar
        = { cigarOp(5, BAM_CHARD_CLIP), cigarOp(10, BAM_CSOFT_CLIP), cigarOp(40, BAM_CMATCH) };
    flankClipCalculator.inspect(flankClippedPair);
    EXPECT_NEAR(50.0 * 2 / 150, flankClipCalculator.estimate(Sex::kFemale).depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, ReadsExtendingPastTheLocusWindow_NotCounted)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    // The fast path's cache keeps a whole pair when EITHER mate is contained in the locus window, so a
    // mate can start inside the right flank yet run past the window's end. The full genotyper never
    // counts such a read, and the `- meanReadLength` term in the denominator assumes it is not counted,
    // so it must be excluded here too. 1081 + 50 = 1131 is one base past the window end of 1130.
    statsCalculator.inspect(makeReadPair("straddlesWindowEnd", 50, 950, 1081));

    const LocusStats stats = statsCalculator.estimate(Sex::kFemale);
    EXPECT_EQ(0, stats.meanFragLength());
    EXPECT_NEAR(50.0 * 1 / 150, stats.depth(), 1e-9);

    // The last start position that still fits inside the window is counted.
    LocusStatsCalculatorFromReadAlignments lastFittingCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    lastFittingCalculator.inspect(makeReadPair("lastFitting", 50, 950, 1080));
    EXPECT_NEAR(50.0 * 2 / 150, lastFittingCalculator.estimate(Sex::kFemale).depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, UnmappedMateOrWrongContig_NotCounted)
{
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);

    FullReadPair readPair = makeReadPair("halfMapped", 50, 950, 1050);
    readPair.secondMate->s.isMapped = false;
    statsCalculator.inspect(readPair);

    FullReadPair otherContigPair = makeReadPair("otherContig", 50, 950, 1050);
    otherContigPair.firstMate->s.chromId = kTestContigIndex + 1;
    otherContigPair.secondMate->s.chromId = kTestContigIndex + 1;
    statsCalculator.inspect(otherContigPair);

    // Only the mapped mate of the first pair contributes, and a pair with one countable mate records no
    // fragment length.
    const LocusStats stats = statsCalculator.estimate(Sex::kFemale);
    EXPECT_EQ(0, stats.meanFragLength());
    EXPECT_NEAR(50.0 * 1 / 150, stats.depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, FragmentLengthRecorded_OnlyForMatesOnTheSameFlank)
{
    LocusStatsCalculatorFromReadAlignments sameFlankCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    // Both mates on the left flank: fragment spans 950 to 980 + 50.
    sameFlankCalculator.inspect(makeReadPair("sameFlank", 50, 950, 980));
    EXPECT_EQ(80, sameFlankCalculator.estimate(Sex::kFemale).meanFragLength());

    // One mate per flank: a pair straddling the repeat has an apparent span that depends on the allele
    // size, so LocusStatsCalculator excludes it and so does this one.
    LocusStatsCalculatorFromReadAlignments oppositeFlanksCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength);
    oppositeFlanksCalculator.inspect(makeReadPair("oppositeFlanks", 50, 950, 1050));
    EXPECT_EQ(0, oppositeFlanksCalculator.estimate(Sex::kFemale).meanFragLength());
}

TEST(LocusStatsCalculatorFromReadAlignments, LeftFlankClampedAtContigStart)
{
    // Repeat at [30, 60) with a 100bp extension: the left flank is clamped to [0, 30), so the denominator
    // is 30 + 100 - 50 rather than 100 + 100 - 50.
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, GenomicRegion(kTestContigIndex, 30, 60), kTestExtensionLength);

    statsCalculator.inspect(makeReadPair("clamped", 50, 10, 70));

    EXPECT_NEAR(50.0 * 2 / 80, statsCalculator.estimate(Sex::kFemale).depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, RightFlankClampedAtContigEnd)
{
    // Repeat at [1000, 1030) with a 100bp extension on a contig only 1080 bases long: the right flank is
    // clamped to [1030, 1080), so the denominator is 100 + 50 - 50 rather than 100 + 100 - 50. Without the
    // clamp the extra 50 nonexistent start positions would push the reported coverage down.
    LocusStatsCalculatorFromReadAlignments statsCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength, 1080);

    statsCalculator.inspect(makeReadPair("nearContigEnd", 50, 950, 1030));

    EXPECT_NEAR(50.0 * 2 / 100, statsCalculator.estimate(Sex::kFemale).depth(), 1e-9);

    // A read that would fit in the unclamped flank but runs off the contig is not counted.
    LocusStatsCalculatorFromReadAlignments pastContigEndCalculator(
        ChromType::kAutosome, kTestRepeatRegion, kTestExtensionLength, 1080);
    pastContigEndCalculator.inspect(makeReadPair("pastContigEnd", 50, 950, 1050));
    EXPECT_NEAR(50.0 * 1 / 100, pastContigEndCalculator.estimate(Sex::kFemale).depth(), 1e-9);
}

TEST(LocusStatsCalculatorFromReadAlignments, ReadLongerThanBothFlanks_DepthIsZeroNotNegative)
{
    // Flanks of 20bp each with 50bp reads leave no valid start positions. LocusStatsCalculator would
    // divide by a negative number here; this calculator reports no coverage instead.
    LocusStatsCalculatorFromReadAlignments statsCalculator(ChromType::kAutosome, kTestRepeatRegion, 20);

    statsCalculator.inspect(makeReadPair("longReads", 50, 985, 1035));

    const LocusStats stats = statsCalculator.estimate(Sex::kFemale);
    EXPECT_EQ(50, stats.meanReadLength());
    EXPECT_DOUBLE_EQ(0.0, stats.depth());
}
