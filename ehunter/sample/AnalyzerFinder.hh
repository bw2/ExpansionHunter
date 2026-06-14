//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <memory>
#include <unordered_map>
#include <vector>

#include "thirdparty/intervaltree/IntervalTree.h"

#include "core/Read.hh"
#include "locus/LocusAnalyzer.hh"
#include "sample/GenomeMask.hh"

namespace ehunter
{

// Specifies which mates should be processed with a given locus analyzer
enum class AnalyzerInputType
{
    kReadOnly,
    kMateOnly,
    kBothReads
};

// Stores information needed to properly pass reads to the analyzer
struct AnalyzerBundle
{
    AnalyzerBundle(locus::RegionType regionType, const size_t initLocusAnalyzerIndex)
        : regionType(regionType)
        , inputType(AnalyzerInputType::kBothReads)
        , locusIndex(initLocusAnalyzerIndex)
    {
    }

    locus::RegionType regionType;
    AnalyzerInputType inputType;
    size_t locusIndex;
};

// Two mates are classified as "nearby" (rather than far-apart) when they map to the same contig and
// their start positions are within this many bases. Shared by AnalyzerFinder (seeking/streaming) and the
// low-mem-streaming / optimized-streaming full-genotyping path so that the two read-pair routing
// implementations cannot drift apart.
constexpr int kMaxMateDistance = 1000;

// True if the two mates are nearby (same contig and start positions within kMaxMateDistance bases).
bool areMatesNearby(int32_t readContigId, int64_t readPosition, int32_t mateContigId, int64_t matePosition);

// True if a read alignment that spans [alignmentStart, alignmentEnd) is fully contained within the
// read-extraction interval [intervalStart, intervalEnd) (start positions and ends are half-open).
bool isAlignmentContainedInInterval(
    int64_t intervalStart, int64_t intervalEnd, int64_t alignmentStart, int64_t alignmentEnd);

void processAnalyzerBundleReadPair(
    locus::LocusAnalyzer& locusAnalyzer, locus::RegionType regionType, AnalyzerInputType inputType, Read& read,
    Read& mate, graphtools::AlignerSelector& alignerSelector);

// Enables retrieval of appropriate locus analyzers by genomic coordinates of read alignments
class AnalyzerFinder
{
public:
    AnalyzerFinder(std::vector<std::unique_ptr<locus::LocusAnalyzer>>& locusAnalyzers);

    AnalyzerFinder(const LocusDescriptionCatalog& locusDescriptions);

    // Retrieves analyzers appropriate for the given read pair
    std::vector<AnalyzerBundle> query(
        int32_t readContigId, int64_t readStart, int64_t readEnd, int32_t mateContigId, int64_t mateStart,
        int64_t mateEnd) const;

    // Retrieves analyzers appropriate for the given read
    std::vector<AnalyzerBundle> query(int32_t readContigId, int64_t readStart, int64_t readEnd) const;

private:
    using AnalyzerIntervalTree = IntervalTree<std::size_t, AnalyzerBundle>;
    using AnalyzerIntervalTrees = std::unordered_map<int32_t, AnalyzerIntervalTree>;

    AnalyzerIntervalTrees intervalTrees_;
};

}
