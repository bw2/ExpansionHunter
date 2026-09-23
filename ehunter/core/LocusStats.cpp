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

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include <htslib/sam.h>

using graphtools::GraphAlignment;

namespace ehunter
{

namespace
{

// Reference position of the read's first base. BAM POS is the first ALIGNED base, so without this a read
// whose leading bases the aligner soft-clipped (typically repeat sequence from an allele longer than the
// reference) would look like it starts in the right flank. The graph aligner aligns the whole read, so it
// starts that read in the repeat and does not count it.
int64_t readStartIncludingSoftClip(const FullRead& read)
{
    int64_t start = read.s.pos;
    for (const uint32_t operation : read.s.cigar)
    {
        if (bam_cigar_op(operation) == BAM_CHARD_CLIP)
        {
            continue;
        }
        if (bam_cigar_op(operation) == BAM_CSOFT_CLIP)
        {
            start -= bam_cigar_oplen(operation);
        }
        break;
    }
    return start;
}

}

bool LocusStats::operator==(const LocusStats& other) const
{
    return alleleCount_ == other.alleleCount_ && meanReadLen_ == other.meanReadLen_
        && meanFragLen_ == other.meanFragLen_ && depth_ == other.depth_;
}

std::ostream& operator<<(std::ostream& out, const LocusStats& stats)
{
    out << "LocusStats(meanReadLength=" << stats.meanReadLength() << ", depth=" << stats.depth() << ")";
    return out;
}

LocusStatsCalculator::LocusStatsCalculator(ChromType chromType, const graphtools::Graph& graph)
    : chromType_(chromType)
{
    // As elsewhere in the program, assuming that the fist and last node are flanks
    leftFlankId_ = 0;
    rightFlankId_ = graph.numNodes() - 1;

    leftFlankLength_ = graph.nodeSeq(leftFlankId_).length();
    rightFlankLength_ = graph.nodeSeq(rightFlankId_).length();
}

void LocusStatsCalculator::inspect(const GraphAlignment& readAlign, const GraphAlignment& mateAlign)
{
    recordReadLen(readAlign);
    recordReadLen(mateAlign);
    recordFragLen(readAlign, mateAlign);
}

void LocusStatsCalculator::inspectRead(const GraphAlignment& readAlign) { recordReadLen(readAlign); }

AlleleCount determineExpectedAlleleCount(ChromType chromType, Sex sex)
{
    switch (chromType)
    {
    case ChromType::kY:
        return AlleleCount::kOne; // Assume that chrY always has copy number one
    case ChromType::kX:
        return (sex == Sex::kFemale ? AlleleCount::kTwo : AlleleCount::kOne);
    case ChromType::kAutosome:
        return AlleleCount::kTwo;
    }

    return AlleleCount::kTwo; // To remove spurious control reaches end of non-void function warning
}

LocusStats LocusStatsCalculator::estimate(Sex sampleSex)
{
    const int readCount = boost::accumulators::count(readLengthAccumulator_);
    AlleleCount alleleCount = determineExpectedAlleleCount(chromType_, sampleSex);

    if (readCount == 0)
    {
        return { alleleCount, 0, 0, 0.0 };
    }

    const int meanReadLength = boost::accumulators::mean(readLengthAccumulator_);
    const int numberOfStartPositions = leftFlankLength_ + rightFlankLength_ - meanReadLength;
    const double depth = meanReadLength * (static_cast<double>(readCount) / numberOfStartPositions);

    int meanFragLen = 0;
    const int fragCount = boost::accumulators::count(fragLengthAccumulator_);
    if (fragCount != 0)
    {
        meanFragLen = boost::accumulators::mean(fragLengthAccumulator_);
    }

    return { alleleCount, meanReadLength, meanFragLen, depth };
}

void LocusStatsCalculator::recordReadLen(const GraphAlignment& readAlign)
{
    const graphtools::NodeId firstNode = readAlign.path().getNodeIdByIndex(0);
    if (firstNode == leftFlankId_ || firstNode == rightFlankId_)
    {
        readLengthAccumulator_(readAlign.queryLength());
    }
}

void LocusStatsCalculator::recordFragLen(const GraphAlignment& readAlign, const GraphAlignment& mateAlign)
{
    const auto readStartNode = readAlign.path().getNodeIdByIndex(0);
    const auto mateStartNode = mateAlign.path().getNodeIdByIndex(0);
    const bool matesStartOnLeftFlank = readStartNode == leftFlankId_ && mateStartNode == leftFlankId_;
    const bool matesStartOnRightFlank = readStartNode == rightFlankId_ && mateStartNode == rightFlankId_;

    if (!matesStartOnLeftFlank && !matesStartOnRightFlank)
    {
        return;
    }

    const int readStart = readAlign.path().startPosition();
    const int readEnd = readStart + static_cast<int>(readAlign.queryLength());

    const int mateStart = mateAlign.path().startPosition();
    const int mateEnd = mateStart + static_cast<int>(mateAlign.queryLength());

    if (readEnd < mateEnd)
    {
        fragLengthAccumulator_(mateEnd - readStart);
    }
    else if (mateEnd < readEnd)
    {
        fragLengthAccumulator_(readEnd - mateStart);
    }
}


LocusStatsCalculatorFromReadAlignments::LocusStatsCalculatorFromReadAlignments(
    ChromType chromType, const GenomicRegion& locusRegion, int regionExtensionLength, int64_t contigLength)
    : chromType_(chromType)
    , locusRegion_(locusRegion)
    // The flanks the full genotyper's graph is built from: regionExtensionLength bases of reference on
    // either side of the repeat (io/LocusSpecDecoding.cpp addFlankingRegions). Clamped to the contig at
    // both ends, since a coordinate outside the contig is not a real reference interval and would inflate
    // the denominator below with start positions no read can occupy. Loci close enough to a contig edge
    // for this to bite are exactly the ones the full genotyper cannot handle at all: it reads the whole
    // unclamped flank, and FastaReference::getSequence throws when that runs off the contig.
    , leftFlankStart_(std::max(int64_t(0), locusRegion.start() - regionExtensionLength))
    , rightFlankEnd_(
          contigLength > 0 ? std::min(contigLength, locusRegion.end() + regionExtensionLength)
                           : locusRegion.end() + regionExtensionLength)
{
}

LocusStatsCalculatorFromReadAlignments::Flank
LocusStatsCalculatorFromReadAlignments::flankAnchoringRead(const FullRead& read) const
{
    if (!read.s.isMapped || read.s.chromId != locusRegion_.contigIndex())
    {
        return Flank::kNone;
    }

    const int64_t alignmentStart = readStartIncludingSoftClip(read);
    const int64_t alignmentEnd = alignmentStart + static_cast<int64_t>(read.r.sequence().size());

    // The read must fit inside the locus window, not merely start in a flank. estimate() below counts on
    // that: its denominator drops one read length precisely because a read starting near the far edge of
    // the right flank runs past the window and is dropped. Seeking and streaming enforce exactly this rule
    // before a read reaches the full genotyper's stats calculator (AnalyzerFinder::query admits only
    // contained reads). Low-mem streaming is looser on both sides -- its cache keeps a whole pair when
    // EITHER mate is contained, and genotypeLocusFull's containment test only picks single-ended vs paired
    // routing for NEARBY pairs, so a far-apart pair passes both mates through regardless. Applying the
    // strict rule here keeps the numerator consistent with the denominator. It matches seeking except for
    // reads with a leading soft clip: seeking tests containment from POS, while this tests it from the read's
    // first base (readStartIncludingSoftClip), so such a read near either edge of the window can be kept by
    // one and dropped by the other.
    if (alignmentStart < leftFlankStart_ || alignmentEnd > rightFlankEnd_)
    {
        return Flank::kNone;
    }

    if (alignmentStart < locusRegion_.start())
    {
        return Flank::kLeft;
    }
    if (locusRegion_.end() <= alignmentStart)
    {
        return Flank::kRight;
    }
    return Flank::kNone;
}

void LocusStatsCalculatorFromReadAlignments::inspect(const FullReadPair& readPair)
{
    // Account for each present, mapped mate independently so that single-ended entries (e.g. reads whose
    // mate is unmapped, used by the optimized-streaming fast path for genotyping) still contribute to the
    // read-length and coverage stats. The fragment-length contribution is only computed when both mates
    // are present and start in the same flank.
    const Flank readFlank = readPair.firstMate ? flankAnchoringRead(*readPair.firstMate) : Flank::kNone;
    const Flank mateFlank = readPair.secondMate ? flankAnchoringRead(*readPair.secondMate) : Flank::kNone;

    if (readFlank != Flank::kNone)
    {
        recordReadLen(readPair.firstMate->r);
    }
    if (mateFlank != Flank::kNone)
    {
        recordReadLen(readPair.secondMate->r);
    }

    // Mirrors LocusStatsCalculator::recordFragLen, which only accepts pairs whose two alignments start on
    // the SAME flank node. That restriction is what keeps the estimate a fragment length rather than the
    // span of a pair straddling the repeat, whose apparent length depends on the allele size.
    if (readFlank != Flank::kNone && readFlank == mateFlank)
    {
        recordFragLen(*readPair.firstMate, *readPair.secondMate);
    }

    // -- RepeatCoverage (disabled) --------------------------------------------------------------------
    // What this function used to accumulate, kept for reference: the read bases overlapping the repeat
    // region itself, which estimate() below divided by the repeat region's length. See the header for why
    // it was replaced by the flank-based figure the full genotyper reports.
    //
    // static int64_t readBasesOverlapInterval(const Read& r, const LinearAlignmentStats& a,
    //                                         const GenomicRegion& region) {
    //     if (a.chromId != region.contigIndex()) {
    //         return 0;
    //     }
    //     const int64_t readStart = a.pos;
    //     const int64_t readEnd = a.pos + r.sequence().size();
    //     return std::max(int64_t(0), std::min(readEnd, region.end()) - std::max(readStart, region.start()));
    // }
    //
    // if (readPair.firstMate && readPair.firstMate->s.isMapped) {
    //     const int64_t readBasesOverlapLocus
    //         = readBasesOverlapInterval(readPair.firstMate->r, readPair.firstMate->s, locusRegion_);
    //     if (readBasesOverlapLocus > 0) {
    //         basesOverlappingLocus_ += readBasesOverlapLocus;
    //     }
    // }
}

void LocusStatsCalculatorFromReadAlignments::recordReadLen(const Read& read)
{
    readLengthAccumulator_(read.sequence().size());
}

void LocusStatsCalculatorFromReadAlignments::recordFragLen(const FullRead& read, const FullRead& mate)
{
    const int64_t readStart = readStartIncludingSoftClip(read);
    const int64_t readEnd = readStart + static_cast<int64_t>(read.r.sequence().size());

    const int64_t mateStart = readStartIncludingSoftClip(mate);
    const int64_t mateEnd = mateStart + static_cast<int64_t>(mate.r.sequence().size());

    if (readEnd < mateEnd)
    {
        fragLengthAccumulator_(mateEnd - readStart);
    }
    else if (mateEnd < readEnd)
    {
        fragLengthAccumulator_(readEnd - mateStart);
    }
}

LocusStats LocusStatsCalculatorFromReadAlignments::estimate(Sex sampleSex)
{
    const int readCount = boost::accumulators::count(readLengthAccumulator_);
    AlleleCount alleleCount = determineExpectedAlleleCount(chromType_, sampleSex);

    if (readCount == 0)
    {
        return { alleleCount, 0, 0, 0.0 };
    }

    const int meanReadLength = boost::accumulators::mean(readLengthAccumulator_);
    const int fragCount = boost::accumulators::count(fragLengthAccumulator_);
    const int meanFragLen = (fragCount != 0) ? boost::accumulators::mean(fragLengthAccumulator_) : 0;

    // Same estimator as LocusStatsCalculator::estimate: the reads counted above each occupy one start
    // position in the flanks, and the number of start positions a countable read can occupy is the two
    // flank lengths less one read length (a read starting within meanReadLength of the far edge of the
    // right flank extends past the locus window and is never collected).
    //
    // Unlike the graph-based calculator this guards the denominator: --region-extension-length small
    // relative to the read length (or a locus at the very start of a contig) would otherwise divide by
    // zero or by a negative number.
    const int64_t leftFlankLength = locusRegion_.start() - leftFlankStart_;
    const int64_t rightFlankLength = rightFlankEnd_ - locusRegion_.end();
    const int64_t numberOfStartPositions = leftFlankLength + rightFlankLength - meanReadLength;
    const double depth = (numberOfStartPositions > 0)
        ? meanReadLength * (static_cast<double>(readCount) / numberOfStartPositions)
        : 0.0;

    return { alleleCount, meanReadLength, meanFragLen, depth };

    // -- RepeatCoverage (disabled) --------------------------------------------------------------------
    // The depth this function used to return, kept for reference: total read bases overlapping the repeat
    // region divided by that region's length, guarded against zero-length loci (e.g. a pure-insertion
    // variant whose ReferenceRegion has start == end).
    //
    // const int locusSize = locusRegion_.end() - locusRegion_.start();
    // const double repeatCoverage
    //     = (locusSize > 0) ? (basesOverlappingLocus_ / static_cast<double>(locusSize)) : 0.0;
}

}
