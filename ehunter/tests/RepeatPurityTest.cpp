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

#include "core/RepeatPurity.hh"

#include "gmock/gmock.h"

using namespace ehunter;

TEST(ComputeRepeatSequencePurity, PerfectRepeat_AllBasesMatch)
{
    const auto purity = computeRepeatSequencePurity("CAGCAGCAGCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 12);
    EXPECT_EQ(purity.totalBases, 12);
}

TEST(ComputeRepeatSequencePurity, PartialTrailingMotif_CountedAgainstPhase)
{
    // "CAGCA" = perfect "CAG" plus the first two bases of the next motif copy; all 5 match.
    const auto purity = computeRepeatSequencePurity("CAGCA", "CAG");
    EXPECT_EQ(purity.matchedBases, 5);
    EXPECT_EQ(purity.totalBases, 5);
}

TEST(ComputeRepeatSequencePurity, SingleSubstitution_OneMismatch)
{
    // Position 4 (second copy, offset 1) is T instead of A.
    const auto purity = computeRepeatSequencePurity("CAGCTGCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 8);
    EXPECT_EQ(purity.totalBases, 9);
}

TEST(ComputeRepeatSequencePurity, IndelShiftsDownstreamPhase)
{
    // A single inserted base ("CAGXCAGCAG") throws off the phase for everything after it, so the
    // flat, gap-free comparison undercounts matches past the insertion. Documents the heuristic's limit.
    const auto purity = computeRepeatSequencePurity("CAGTCAGCAG", "CAG");
    EXPECT_EQ(purity.totalBases, 10);
    EXPECT_LT(purity.matchedBases, 9);
}

TEST(ComputeRepeatSequencePurity, LowercaseSoftMaskedReference_Matches)
{
    const auto purity = computeRepeatSequencePurity("cagCAGcag", "CAG");
    EXPECT_EQ(purity.matchedBases, 9);
    EXPECT_EQ(purity.totalBases, 9);
}

TEST(ComputeRepeatSequencePurity, NBaseCountsAsMismatch)
{
    const auto purity = computeRepeatSequencePurity("CANCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 5);
    EXPECT_EQ(purity.totalBases, 6);
}

TEST(ComputeRepeatSequencePurity, EmptyMotif_ReturnsZeroTotal)
{
    const auto purity = computeRepeatSequencePurity("CAGCAG", "");
    EXPECT_EQ(purity.matchedBases, 0);
    EXPECT_EQ(purity.totalBases, 0);
}

TEST(ComputeRepeatSequencePurity, EmptySequence_ReturnsZeroTotal)
{
    const auto purity = computeRepeatSequencePurity("", "CAG");
    EXPECT_EQ(purity.matchedBases, 0);
    EXPECT_EQ(purity.totalBases, 0);
}
