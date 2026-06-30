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

#include "sample/HtsLowMemStreamingHelpers.hh"

#include <string>
#include <vector>

#include <htslib/sam.h>

#include "gtest/gtest.h"

#include "core/Read.hh"

using namespace ehunter;

namespace
{

// Encode one CIGAR operation into htslib's packed uint32_t representation.
uint32_t cigarOp(int length, int op)
{
    return (static_cast<uint32_t>(length) << BAM_CIGAR_SHIFT) | static_cast<uint32_t>(op);
}

// Build a minimal FullRead carrying only the fields processRead() reads:
// the query sequence, the 0-based alignment start (pos), and the CIGAR.
FullRead makeRead(const std::string& sequence, int64_t pos, const std::vector<uint32_t>& cigar)
{
    Read read(ReadId("frag", MateNumber::kFirstMate), sequence, false);
    LinearAlignmentStats stats;
    stats.chromId = 0;
    stats.pos = static_cast<int32_t>(pos);
    stats.isMapped = true;
    stats.cigar = cigar;
    return FullRead(std::move(read), std::move(stats));
}

} // namespace

// Worked example: a fully-spanning read, single 30M op, no soft clips.
// Read aligns to ref [90, 120); locus is the 0-based half-open interval [100, 105) (5 bp).
// The 5 repeat bases sit at read-seq indices 10..14, flanked by non-repetitive bases.
// Confirms repeat_sequence_size == 5 (the pre-fix code reported 6 — off-by-one).
TEST(ProcessRead, SpanningReadCountsExactRepeatBases)
{
    // Whole read is 'T's: the M branch counts overlap by position, not by content,
    // and 'T' never matches the "CAG" motif so no flank/soft-clip extension occurs.
    const std::string sequence(30, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(30, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/105, "CAG");

    EXPECT_TRUE(result.overlaps_repeats);
    EXPECT_EQ(5, result.repeat_sequence_size_in_base_pairs);
    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_FALSE(result.soft_clipped_bases_contain_repetitive_sequence);
}

// Soft-clip case: a flanking read whose aligned portion enters the repeat and whose
// soft-clipped tail is pure repeat. Read aligns to ref [90, 110) with CIGAR 20M9S;
// locus is [100, 130). The 20M contributes 10 bp of in-repeat sequence (ref 100..109),
// and the 9 bp soft clip "CAGCAGCAG" adds 3 motif copies (9 bp) via flank extension.
// Total = 19 bp; the soft clip is flagged as repetitive; the read does not span the
// right edge so is_spanning_read is false.
TEST(ProcessRead, SoftClippedRepeatBasesAreCounted)
{
    const std::string sequence = std::string(20, 'T') + "CAGCAGCAG";
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(20, BAM_CMATCH), cigarOp(9, BAM_CSOFT_CLIP)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.overlaps_repeats);
    EXPECT_EQ(19, result.repeat_sequence_size_in_base_pairs);
    EXPECT_FALSE(result.is_spanning_read);
    EXPECT_TRUE(result.soft_clipped_bases_contain_repetitive_sequence);
}

// Insertion exactly at the locus boundary: locus is [100, 130) (half-open, so the
// repeat's last base is ref 129 and ref 130 is the first base past it). A 6 bp insertion
// sits at ref 130 — just outside the locus but within the ±motif_size padding the code
// applies for insertions adjacent to the repeat. CIGAR is 40M6I10M aligned at ref 90.
// The 40M contributes 30 bp (ref 100..129), the boundary insertion adds 6 bp (length is
// a multiple of the 3 bp motif), and the trailing 10M lies entirely past the locus.
// Total = 36 bp.
TEST(ProcessRead, InsertionAtLocusBoundaryIsCounted)
{
    const std::string sequence(56, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(40, BAM_CMATCH), cigarOp(6, BAM_CINS), cigarOp(10, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.overlaps_repeats);
    EXPECT_EQ(36, result.repeat_sequence_size_in_base_pairs);
    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_FALSE(result.soft_clipped_bases_contain_repetitive_sequence);
}

// Companion to the boundary test: an insertion one base beyond the ±motif_size padding
// must NOT be counted. With locus [100, 130) and a 3 bp motif the padding extends the
// closed right bound to ref 133, so a 6 bp insertion at ref 134 is excluded. CIGAR is
// 44M6I6M aligned at ref 90; only the 44M's 30 in-repeat bases (ref 100..129) count.
// Total = 30 bp (not 36). This pins the right boundary and guards against re-widening it.
TEST(ProcessRead, InsertionJustBeyondPaddingIsNotCounted)
{
    const std::string sequence(56, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(44, BAM_CMATCH), cigarOp(6, BAM_CINS), cigarOp(6, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.overlaps_repeats);
    EXPECT_EQ(30, result.repeat_sequence_size_in_base_pairs);
    EXPECT_FALSE(result.soft_clipped_bases_contain_repetitive_sequence);
}

// ReadRepeatPurity: a spanning read whose in-repeat tract is a perfect "CAG" tiling.
// Read [90,120) with 30M; locus [100, 106) (6 bp). The 6 in-repeat read bases (indices
// 10..15) are "CAGCAG" -> all 6 match the motif, so matched == total == 6.
TEST(ProcessRead, PureSpanningReadRepeatPurityIsOne)
{
    const std::string sequence = std::string(10, 'T') + "CAGCAG" + std::string(14, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(30, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/106, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(6, result.repeat_read_bases);
    EXPECT_EQ(6, result.matched_bases_within_repeat);
}

// ReadRepeatPurity: same layout but the in-repeat tract "CAGCTG" carries one substitution
// (position 4 is T instead of A), so 5 of 6 bases match the motif tiling.
TEST(ProcessRead, ImpureSpanningReadRepeatPurityCountsMismatch)
{
    const std::string sequence = std::string(10, 'T') + "CAGCTG" + std::string(14, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(30, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/106, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(6, result.repeat_read_bases);
    EXPECT_EQ(5, result.matched_bases_within_repeat);
}
