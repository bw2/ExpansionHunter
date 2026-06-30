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

TEST(MotifTilingPurity, PerfectRepeat_AllBasesMatch)
{
    const auto purity = motifTilingPurity("CAGCAGCAGCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 12);
    EXPECT_EQ(purity.totalBases, 12);
}

TEST(MotifTilingPurity, PartialTrailingMotif_CountedAgainstPhase)
{
    // "CAGCA" = perfect "CAG" plus the first two bases of the next motif copy; all 5 match.
    const auto purity = motifTilingPurity("CAGCA", "CAG");
    EXPECT_EQ(purity.matchedBases, 5);
    EXPECT_EQ(purity.totalBases, 5);
}

TEST(MotifTilingPurity, SingleSubstitution_OneMismatch)
{
    // Position 4 (second copy, offset 1) is T instead of A.
    const auto purity = motifTilingPurity("CAGCTGCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 8);
    EXPECT_EQ(purity.totalBases, 9);
}

TEST(MotifTilingPurity, IndelShiftsDownstreamPhase)
{
    // A single inserted base ("CAGXCAGCAG") throws off the phase for everything after it, so the
    // flat tiling comparison undercounts matches past the insertion. Documents the heuristic's limit.
    const auto purity = motifTilingPurity("CAGTCAGCAG", "CAG");
    EXPECT_EQ(purity.totalBases, 10);
    EXPECT_LT(purity.matchedBases, 9);
}

TEST(MotifTilingPurity, LowercaseSoftMaskedReference_Matches)
{
    const auto purity = motifTilingPurity("cagCAGcag", "CAG");
    EXPECT_EQ(purity.matchedBases, 9);
    EXPECT_EQ(purity.totalBases, 9);
}

TEST(MotifTilingPurity, NBaseCountsAsMismatch)
{
    const auto purity = motifTilingPurity("CANCAG", "CAG");
    EXPECT_EQ(purity.matchedBases, 5);
    EXPECT_EQ(purity.totalBases, 6);
}

TEST(MotifTilingPurity, EmptyMotif_ReturnsZeroTotal)
{
    const auto purity = motifTilingPurity("CAGCAG", "");
    EXPECT_EQ(purity.matchedBases, 0);
    EXPECT_EQ(purity.totalBases, 0);
}

TEST(MotifTilingPurity, EmptySequence_ReturnsZeroTotal)
{
    const auto purity = motifTilingPurity("", "CAG");
    EXPECT_EQ(purity.matchedBases, 0);
    EXPECT_EQ(purity.totalBases, 0);
}
