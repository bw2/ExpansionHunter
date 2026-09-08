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

#pragma once

#include <iostream>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Read.hh"
#include "core/Reference.hh"

namespace ehunter
{

class LocusStats
{
public:
    LocusStats(
        AlleleCount alleleCount = AlleleCount::kOne, int meanReadLen = 0, int meanFragLen = 0, double depth = 0)
        : alleleCount_(alleleCount)
        , meanReadLen_(meanReadLen)
        , meanFragLen_(meanFragLen)
        , depth_(depth)
    {
    }

    AlleleCount alleleCount() const { return alleleCount_; }
    int meanReadLength() const { return meanReadLen_; }
    int meanFragLength() const { return meanFragLen_; }
    double depth() const { return depth_; }
    void setDepth(double depth) { depth_ = depth; }

    bool operator==(const LocusStats& other) const;

private:
    AlleleCount alleleCount_;
    int meanReadLen_;
    int meanFragLen_;
    double depth_;
};

std::ostream& operator<<(std::ostream& out, const LocusStats& stats);

// Computes read and coverage statistics for each locus from reads aligning to the flanks
class LocusStatsCalculator
{
public:
    LocusStatsCalculator(ChromType chromType, const graphtools::Graph& graph);

    void inspect(const graphtools::GraphAlignment& readAlign, const graphtools::GraphAlignment& mateAlign);
    void inspectRead(const graphtools::GraphAlignment& readAlign);

    LocusStats estimate(Sex sampleSex);
    void recordReadLen(const graphtools::GraphAlignment& readAlign);

private:
    using AccumulatorStats
        = boost::accumulators::features<boost::accumulators::tag::count, boost::accumulators::tag::mean>;
    using Accumulator = boost::accumulators::accumulator_set<int, AccumulatorStats>;

    void recordFragLen(const graphtools::GraphAlignment& readAlign, const graphtools::GraphAlignment& mateAlign);

    ChromType chromType_;
    Accumulator readLengthAccumulator_;
    Accumulator fragLengthAccumulator_;
    graphtools::NodeId leftFlankId_;
    graphtools::NodeId rightFlankId_;
    int leftFlankLength_;
    int rightFlankLength_;
};


AlleleCount determineExpectedAlleleCount(ChromType chromType, Sex sex);


// Computes the same KIND of read and coverage statistics as LocusStatsCalculator, but from the reads'
// linear (BAM/CRAM) alignments instead of their graph alignments, for the optimized-streaming fast path
// which performs no read-to-graph realignment.
//
// What is mirrored from LocusStatsCalculator::estimate, so that the emitted Coverage field measures the
// same quantity in both modes: the flank window, the rule that only reads anchored in a flank are
// counted, and the number-of-start-positions denominator. `s.pos` inside a flank interval is the linear
// stand-in for the graph rule "the alignment path's first node is a flank node".
//
// What is NOT mirrored, so the two will not agree read-for-read: LocusStatsCalculator only ever sees
// reads the graph aligner accepted (LocusAnalyzer::processOntargetMates skips a mate whose alignment came
// back empty), while this calculator has no realignment to gate on and counts every mapped in-window read
// anchored in a flank. It also runs before the fast path's own supplementary/secondary/mapQ screen, which
// is deliberate: the full genotyper applies no such screen to its stats either.
class LocusStatsCalculatorFromReadAlignments
{
public:
    // contigLength bounds the right flank so it cannot run off the end of the contig; pass 0 when it is
    // not known, which leaves the right flank unbounded.
    LocusStatsCalculatorFromReadAlignments(
        ChromType chromType, const GenomicRegion& locusRegion, int regionExtensionLength,
        int64_t contigLength = 0);

    void inspect(const FullReadPair& readPair);
    LocusStats estimate(Sex sampleSex);

private:
    using AccumulatorStats
        = boost::accumulators::features<boost::accumulators::tag::count, boost::accumulators::tag::mean>;
    using Accumulator = boost::accumulators::accumulator_set<int, AccumulatorStats>;

    // Which flank a read's linear alignment starts in, if any.
    enum class Flank
    {
        kNone,
        kLeft,
        kRight
    };

    // The flank a read is anchored in: its alignment starts inside that flank AND the whole read fits
    // inside the locus window. Reads that fail either test are not counted.
    Flank flankAnchoringRead(const FullRead& read) const;
    void recordReadLen(const Read& read);
    void recordFragLen(const FullRead& read, const FullRead& mate);

    ChromType chromType_;
    Accumulator readLengthAccumulator_;
    Accumulator fragLengthAccumulator_;
    GenomicRegion locusRegion_;
    // Reference flanks of the locus, matching the graph's first and last node under full genotyping
    // (io/LocusSpecDecoding.cpp addFlankingRegions). Half-open, clamped to the start of the contig.
    int64_t leftFlankStart_;
    int64_t rightFlankEnd_;

    // -- RepeatCoverage (disabled) --------------------------------------------------------------------
    // The fast path originally reported depth over the REPEAT REGION itself: the read bases overlapping
    // locusRegion_ divided by that region's length. That is a different quantity from the flank coverage
    // the full genotyper reports under the same Coverage field, on a denominator of tens of bp rather than
    // ~2 * regionExtensionLength, so the two branches disagreed systematically. Kept commented out, under
    // the name it would carry if it were ever emitted alongside the flank-based Coverage.
    // unsigned int basesOverlappingLocus_;
};

}
