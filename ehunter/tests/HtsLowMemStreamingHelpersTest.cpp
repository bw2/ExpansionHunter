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

// ReadRepeatPurity: a spanning read whose in-repeat tract is a perfect consecutive "CAG" repeat.
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
// (position 4 is T instead of A), so 5 of 6 bases match the motif.
TEST(ProcessRead, ImpureSpanningReadRepeatPurityCountsMismatch)
{
    const std::string sequence = std::string(10, 'T') + "CAGCTG" + std::string(14, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(30, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/106, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(6, result.repeat_read_bases);
    EXPECT_EQ(5, result.matched_bases_within_repeat);
}

// Consensus tract: a clean spanning read with a single M op. Read [90, 140) with 50M; locus
// [100, 130). Both locus edges fall inside the M op, so the tract is read-seq [10, 40) -- the
// 30 bases the read carries between the locus edges -- and it equals the vote size.
TEST(ProcessRead, SpanningReadExportsRepeatTract)
{
    const std::string sequence = std::string(10, 'T') + std::string(10, 'C') + std::string(30, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(50, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(30, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(10, result.repeat_tract_start);
    EXPECT_EQ(30, result.repeat_tract_length);
    EXPECT_EQ(sequence.substr(10, 30), sequence.substr(result.repeat_tract_start, result.repeat_tract_length));
}

// A whole-motif insertion sitting exactly on the RIGHT locus edge is an expansion the aligner pushed
// just outside the repeat, so the tract must swallow it -- otherwise an expanded allele would report
// the reference-length sequence. Same read as InsertionAtLocusBoundaryIsCounted (40M6I10M at ref 90,
// locus [100, 130), 3 bp motif): the vote counts the 6 bp insertion via the padding rule, and the
// tract must agree at 36 bp rather than stopping at the 30 aligned bases.
TEST(ProcessRead, RepeatTractSwallowsAWholeMotifInsertionAtTheRightLocusEdge)
{
    const std::string sequence(56, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(40, BAM_CMATCH), cigarOp(6, BAM_CINS), cigarOp(10, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(36, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(10, result.repeat_tract_start);
    EXPECT_EQ(36, result.repeat_tract_length);
}

// Mirror of the above for the LEFT locus edge: 10M6I40M at ref 90 with locus [100, 130) puts a 6 bp
// whole-motif insertion at ref 100. The vote counts it (30 aligned bases + 6 inserted = 36), so the
// tract must start at the first inserted base rather than after it.
TEST(ProcessRead, RepeatTractSwallowsAWholeMotifInsertionAtTheLeftLocusEdge)
{
    const std::string sequence(56, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(10, BAM_CMATCH), cigarOp(6, BAM_CINS), cigarOp(40, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(36, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(10, result.repeat_tract_start);
    EXPECT_EQ(36, result.repeat_tract_length);
}

// A non-whole-motif insertion at a locus edge stays OUT of the tract, because the vote does not count
// it either (the BAM_CINS branch requires a whole number of motifs). 40M1I10M at ref 90 with locus
// [100, 130) and a 3 bp motif: vote 30, tract 30.
TEST(ProcessRead, RepeatTractIgnoresANonWholeMotifInsertionAtTheLocusEdge)
{
    const std::string sequence(51, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(40, BAM_CMATCH), cigarOp(1, BAM_CINS), cigarOp(10, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(30, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(30, result.repeat_tract_length);
}

// A locus edge landing inside a deletion still yields a usable tract: the read has no base at the
// deleted reference positions, and a deletion consumes no read bases, so the read position at the D
// op is already the first base past the deletion. CIGAR 8M6D42M at ref 90 puts the 6 bp deletion over
// ref [98, 104), straddling the locus start (ref 100). The tract runs from read index 8 to the locus
// end, 26 bases -- exactly the vote size, since the vote also credits a deletion with nothing.
TEST(ProcessRead, RepeatTractSurvivesADeletionStraddlingTheLocusStart)
{
    const std::string sequence(50, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(8, BAM_CMATCH), cigarOp(6, BAM_CDEL), cigarOp(42, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(26, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(8, result.repeat_tract_start);
    EXPECT_EQ(26, result.repeat_tract_length);
}

// Mirror of the above for the locus END. CIGAR 38M6D12M at ref 90 puts the deletion over ref
// [128, 134), straddling the locus end (ref 130). The tract stops at the last read base before the
// deletion (read index 38), so it spans read [10, 38) = 28 bases, matching the vote.
TEST(ProcessRead, RepeatTractSurvivesADeletionStraddlingTheLocusEnd)
{
    const std::string sequence(50, 'T');
    FullRead read = makeRead(
        sequence, /*pos=*/90, {cigarOp(38, BAM_CMATCH), cigarOp(6, BAM_CDEL), cigarOp(12, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(28, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(10, result.repeat_tract_start);
    EXPECT_EQ(28, result.repeat_tract_length);
}

// Regression: a deletion that ENDS exactly at the locus start must not steal the tract's start edge.
// CIGAR 93M1D55M at ref 11347586 with locus [11347680, 11347687) (motif "C"): the 1 bp deletion covers
// ref [11347679, 11347680), i.e. entirely before the half-open locus, and the following 55M is what
// actually carries the repeat. Modelled on read HISEQ1:18:H8VC6ADXX:1:2101:2771:5050 in
// profile/cache_test/test.bam. A closed-interval start test would let the deletion claim the edge and
// drop this read from the consensus even though it votes for the genotype.
TEST(ProcessRead, DeletionEndingAtLocusStartDoesNotSuppressTheTract)
{
    const std::string sequence = std::string(93, 'A') + "CCCCCCC" + std::string(48, 'A');
    FullRead read = makeRead(
        sequence, /*pos=*/11347586,
        {cigarOp(93, BAM_CMATCH), cigarOp(1, BAM_CDEL), cigarOp(55, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/11347680, /*locus_end=*/11347687, "C");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(93, result.repeat_tract_start);
    EXPECT_EQ(7, result.repeat_tract_length);
    EXPECT_EQ("CCCCCCC", sequence.substr(result.repeat_tract_start, result.repeat_tract_length));
}

// Regression: a locus whose reference length is not a whole multiple of the motif still yields a tract,
// and that tract is the full locus length (11 bp), NOT the truncated whole-motif length (10 bp). The
// consensus lookup must therefore match buckets to alleles by unit count -- 11 / 2 == 5 units -- rather
// than by alleleUnits * motifSize, which would look for a 10 bp bucket that can never exist here.
// Modelled on catalog locus 1-479076-479087-CA (chr1:479076-479087, 11 bp, motif "CA").
TEST(ProcessRead, RepeatTractKeepsTheRemainderWhenLocusIsNotAWholeNumberOfMotifs)
{
    const std::string sequence = std::string(10, 'T') + "CACACACACAC" + std::string(9, 'T');
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(30, BAM_CMATCH)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/111, "CA");

    EXPECT_TRUE(result.is_spanning_read);
    EXPECT_EQ(11, result.repeat_sequence_size_in_base_pairs);
    EXPECT_EQ(10, result.repeat_tract_start);
    EXPECT_EQ(11, result.repeat_tract_length);
    // 11 / 2 == 5 units, so a bucket lookup on 5 * 2 == 10 bp would miss this read entirely.
    EXPECT_NE(result.repeat_tract_length, (result.repeat_sequence_size_in_base_pairs / 2) * 2);
}

// A non-spanning read that never reaches the right locus edge has no end anchor, so no tract.
// CIGAR 20M9S at ref 90 with locus [100, 130): the aligned portion ends at ref 110.
TEST(ProcessRead, RepeatTractIsUnusableWhenLocusEndIsNotReached)
{
    const std::string sequence = std::string(20, 'T') + "CAGCAGCAG";
    FullRead read = makeRead(sequence, /*pos=*/90, {cigarOp(20, BAM_CMATCH), cigarOp(9, BAM_CSOFT_CLIP)});

    FastReadAnalysisResult result = processRead(read, /*locus_start=*/100, /*locus_end=*/130, "CAG");

    EXPECT_FALSE(result.is_spanning_read);
    EXPECT_EQ(-1, result.repeat_tract_length);
}
