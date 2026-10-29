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

#include <algorithm>
#include <memory>
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

namespace
{
std::vector<std::shared_ptr<FullReadPair>> makeReadPairs(const std::vector<std::pair<std::string, std::string>>& sequences)
{
    std::vector<std::shared_ptr<FullReadPair>> readPairs;
    for (const auto& [firstMateSequence, secondMateSequence] : sequences) {
        auto readPair = std::make_shared<FullReadPair>();
        readPair->firstMate = makeRead(firstMateSequence, /*pos=*/100, {cigarOp(firstMateSequence.size(), BAM_CMATCH)});
        if (!secondMateSequence.empty()) {
            readPair->secondMate = makeRead(secondMateSequence, /*pos=*/400, {cigarOp(secondMateSequence.size(), BAM_CMATCH)});
        }
        readPairs.push_back(readPair);
    }
    return readPairs;
}

std::string repeat(const std::string& motif, int count)
{
    std::string sequence;
    for (int i = 0; i < count; ++i) {
        sequence += motif;
    }
    return sequence;
}
} // namespace

// Two in-repeat reads anchored by flanking mates are the signature of an allele longer than a read, so the
// fast path must hand the locus to the full genotyper. The second one is on the other strand (CTG for a CAG
// locus) and starts mid-motif, which the purity test also accepts.
TEST(HasEnoughInRepeatReadsForFullGenotyping, TwoInRepeatReadsAreEnough)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    const auto readPairs = makeReadPairs({
        {flank, repeat("CAG", 50)},
        {flank, "TG" + repeat("CTG", 49) + "C"},
        {flank, flank}});
    EXPECT_TRUE(hasEnoughInRepeatReadsForFullGenotyping(readPairs, "CAG"));
}

// On two-color sequencers, reads whose signal fails become runs of G, which look like in-repeat reads at C and G
// homopolymers (G is the reverse complement of C), so homopolymers never go to the full genotyper this way.
TEST(HasEnoughInRepeatReadsForFullGenotyping, HomopolymersNeverQualify)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    const auto readPairs = makeReadPairs({{flank, repeat("G", 150)}, {flank, repeat("G", 150)}});
    EXPECT_FALSE(hasEnoughInRepeatReadsForFullGenotyping(readPairs, "G"));
    EXPECT_FALSE(hasEnoughInRepeatReadsForFullGenotyping(readPairs, "C"));
}

TEST(HasEnoughInRepeatReadsForFullGenotyping, OneInRepeatReadIsNotEnough)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    const auto readPairs = makeReadPairs({{flank, repeat("CAG", 50)}, {flank, flank}, {flank, ""}});
    EXPECT_FALSE(hasEnoughInRepeatReadsForFullGenotyping(readPairs, "CAG"));
}

// A short locus inside a reference stretch of the motif longer than a read: such loci get in-repeat reads
// from the reference itself, so the fast path keeps them.
TEST(HasReadLengthStretchOfInRepeatSequence, LongReferenceRepeatIsFound)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    EXPECT_TRUE(hasReadLengthStretchOfInRepeatSequence(flank + repeat("ACC", 70) + flank, "CCA", 150));
}

// An imperfect VNTR like chr22:20,162,000-20,164,800 (ACC with an ATCACT every 21 bases): reads from it pass the
// in-repeat read test, so it must count as a stretch of repeat too.
TEST(HasReadLengthStretchOfInRepeatSequence, ImperfectReferenceRepeatIsFound)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    EXPECT_TRUE(hasReadLengthStretchOfInRepeatSequence(flank + repeat("ACCACCACCACCACCATCACT", 12) + flank, "ACC", 150));
}

TEST(HasReadLengthStretchOfInRepeatSequence, ShortReferenceRepeatIsNotFound)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    EXPECT_FALSE(hasReadLengthStretchOfInRepeatSequence(flank + flank + repeat("CAG", 10) + flank + flank, "CAG", 150));
    EXPECT_FALSE(hasReadLengthStretchOfInRepeatSequence(repeat("CAG", 10), "CAG", 150));
}

// Reads that are only partly repeat (a flank plus a normal-length allele) don't count.
TEST(HasEnoughInRepeatReadsForFullGenotyping, ReadsWithFlanksAreNotInRepeatReads)
{
    const std::string flank = "ACGTTGCAAGTCCATGGTACCTTAGACGTACGATCGGATCCTAGCTAGGTACCATGGTTAACG";
    const auto readPairs = makeReadPairs({
        {flank + repeat("CAG", 20), repeat("CAG", 20) + flank},
        {flank + repeat("CAG", 25), repeat("CAG", 25) + flank}});
    EXPECT_FALSE(hasEnoughInRepeatReadsForFullGenotyping(readPairs, "CAG"));
}

// Reads are sized in bp, so with a 12 bp motif reads of 24, 30 and 35 bp all count toward 2 repeat units. Their counts
// must add up rather than the last size's count replacing the others (CSTB once reported 1 spanning read instead of 31).
TEST(CountReadsByRepeatSizeInUnits, CountsOfSizesInTheSameUnitAddUp)
{
    const CountTable readCountsByRepeatSizeInUnits = countReadsByRepeatSizeInUnits({{24, 28}, {30, 2}, {35, 1}, {36, 4}}, 12);
    EXPECT_EQ(31, readCountsByRepeatSizeInUnits.countOf(2));
    EXPECT_EQ(4, readCountsByRepeatSizeInUnits.countOf(3));
    EXPECT_EQ(std::vector<int32_t>({2, 3}), readCountsByRepeatSizeInUnits.getElementsWithNonzeroCounts());
}

namespace
{
// A read pair (read 1 only) of fragment `fragmentId`, with read 1 flagged as a duplicate when isFlaggedAsDuplicate.
std::shared_ptr<FullReadPair> makeReadPairForLocusCache(const std::string& fragmentId, bool isFlaggedAsDuplicate)
{
    LinearAlignmentStats stats;
    stats.isMapped = true;
    stats.isFlaggedAsDuplicate = isFlaggedAsDuplicate;
    auto readPair = std::make_shared<FullReadPair>();
    readPair->firstMate = FullRead(Read(ReadId(fragmentId, MateNumber::kFirstMate), "ACGT", false), stats);
    return readPair;
}

// Offers a read pair the way the optimized-streaming loop does: admitFragment, then offerReadPair if admitted.
void offerToLocusCache(LocusCache& locusCache, const std::shared_ptr<FullReadPair>& readPair)
{
    const bool isFlaggedPair = isFlaggedAsDuplicateOrFailingQualityChecks(readPair->firstMate->s);
    if (locusCache.admitFragment(readPair->firstMate->r.fragmentId(), isFlaggedPair, locusCache.maxReadPairs()))
    {
        locusCache.offerReadPair(readPair, isFlaggedPair, locusCache.maxReadPairs());
    }
}

int countFlaggedReadPairsInLocusCache(const LocusCache& locusCache)
{
    int count = 0;
    for (const auto& readPair : locusCache.readPairs)
    {
        count += readPair->firstMate->s.isFlaggedAsDuplicate ? 1 : 0;
    }
    return count;
}

std::vector<std::string> sortedFragmentIdsInLocusCache(const LocusCache& locusCache)
{
    std::vector<std::string> fragmentIds;
    for (const auto& readPair : locusCache.readPairs)
    {
        fragmentIds.push_back(readPair->firstMate->r.fragmentId());
    }
    std::sort(fragmentIds.begin(), fragmentIds.end());
    return fragmentIds;
}

// Depth 10 over a 50 bp window with 50 bp reads is 10 * (50 + 50) / 50 = 20 reads, so a cap of 10 read pairs.
LocusCache makeLocusCacheWithCapOf10ReadPairs() { return LocusCache(0, GenomicRegion(0, 0, 50), "LOCUS", 10, 50); }
} // namespace

// Past the cap, read pairs without a duplicate or QC-fail flag take the slots of flagged pairs, so a pile of flagged
// pairs (like the ~5,600 duplicate-flagged pairs seen at ATXN7 in one sample) can't crowd out the other reads.
TEST(LocusCache, UnflaggedReadPairsDisplaceFlaggedReadPairsPastTheCap)
{
    LocusCache locusCache = makeLocusCacheWithCapOf10ReadPairs();
    ASSERT_EQ(10u, locusCache.maxReadPairs());
    int unflaggedCount = 0;
    for (int i = 0; i < 100; ++i)
    {
        const bool isFlaggedAsDuplicate = i % 13 != 0;
        unflaggedCount += isFlaggedAsDuplicate ? 0 : 1;
        offerToLocusCache(locusCache, makeReadPairForLocusCache("fragment" + std::to_string(i), isFlaggedAsDuplicate));
    }
    ASSERT_EQ(8, unflaggedCount);
    EXPECT_EQ(10u, locusCache.readPairs.size());
    EXPECT_EQ(2, countFlaggedReadPairsInLocusCache(locusCache));
    EXPECT_TRUE(locusCache.reservoirSampled());
    EXPECT_EQ(100u, locusCache.offeredCount());
}

// With more unflagged read pairs than slots, no flagged pair is kept, and each kept pair is a different fragment.
TEST(LocusCache, FlaggedReadPairsAreDroppedWhenThereAreEnoughUnflaggedReadPairs)
{
    LocusCache locusCache = makeLocusCacheWithCapOf10ReadPairs();
    for (int i = 0; i < 100; ++i)
    {
        offerToLocusCache(locusCache, makeReadPairForLocusCache("fragment" + std::to_string(i), i % 3 != 0));
    }
    EXPECT_EQ(0, countFlaggedReadPairsInLocusCache(locusCache));
    std::vector<std::string> fragmentIds = sortedFragmentIdsInLocusCache(locusCache);
    EXPECT_EQ(10u, fragmentIds.size());
    EXPECT_EQ(fragmentIds.end(), std::adjacent_find(fragmentIds.begin(), fragmentIds.end()));
}

// Without flagged read pairs, and with each fragment offered once, the sample is plain reservoir sampling (Algorithm R)
// with the locus's random number generator, the same as before the flags were considered.
TEST(LocusCache, WithoutFlaggedReadPairsTheSampleIsPlainReservoirSampling)
{
    LocusCache locusCache = makeLocusCacheWithCapOf10ReadPairs();
    uint64_t rngState = fnv1a64("LOCUS");
    auto nextRandom = [&rngState](uint64_t bound) {
        rngState += 0x9E3779B97F4A7C15ULL;
        uint64_t z = rngState;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        return z % bound;
    };
    std::vector<std::string> expectedFragmentIds;
    for (int i = 0; i < 200; ++i)
    {
        const std::string fragmentId = "fragment" + std::to_string(i);
        offerToLocusCache(locusCache, makeReadPairForLocusCache(fragmentId, false));
        if (expectedFragmentIds.size() < 10)
        {
            expectedFragmentIds.push_back(fragmentId);
        }
        else
        {
            const uint64_t j = nextRandom(i + 1);
            if (j < 10)
            {
                expectedFragmentIds[j] = fragmentId;
            }
        }
    }
    std::sort(expectedFragmentIds.begin(), expectedFragmentIds.end());
    EXPECT_EQ(expectedFragmentIds, sortedFragmentIdsInLocusCache(locusCache));
}

// A fragment offered a second time past the cap (a far-away pair whose two ends are both in the locus window) is not
// admitted again when its first offer took a flagged pair's slot, so it can't fill two slots.
TEST(LocusCache, AFragmentThatDisplacedAFlaggedReadPairIsNotAdmittedAgain)
{
    LocusCache locusCache = makeLocusCacheWithCapOf10ReadPairs();
    for (int i = 0; i < 10; ++i)
    {
        offerToLocusCache(locusCache, makeReadPairForLocusCache("flagged" + std::to_string(i), true));
    }
    offerToLocusCache(locusCache, makeReadPairForLocusCache("unflagged", false));
    EXPECT_EQ(9, countFlaggedReadPairsInLocusCache(locusCache));
    EXPECT_FALSE(locusCache.admitFragment("unflagged", false, locusCache.maxReadPairs()));
}

// A fragment first offered below the cap is not admitted again once the cap is reached (the other end of a far-away
// pair arriving later), so it isn't counted twice and doesn't take part in the sampling a second time.
TEST(LocusCache, AFragmentOfferedBelowTheCapIsNotAdmittedAgainPastIt)
{
    LocusCache locusCache = makeLocusCacheWithCapOf10ReadPairs();
    for (int i = 0; i < 20; ++i)
    {
        offerToLocusCache(locusCache, makeReadPairForLocusCache("fragment" + std::to_string(i), false));
    }
    EXPECT_FALSE(locusCache.admitFragment("fragment0", false, locusCache.maxReadPairs()));
    EXPECT_TRUE(locusCache.admitFragment("fragment20", false, locusCache.maxReadPairs()));
    EXPECT_EQ(20u, locusCache.offeredCount());
}
