//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#include "locus/LocusAligner.hh"

#include "alignment/AlignmentFilters.hh"
#include "alignment/OperationsOnAlignments.hh"

namespace ehunter
{
namespace locus
{

LocusAligner::LocusAligner(
    std::string locusId, GraphPtr graph, const HeuristicParameters& params, BamletWriterPtr writer,
    AlignmentBufferPtr buffer)
    : locusId_(std::move(locusId))
    , aligner_(graph, params.kmerLenForAlignment(), params.paddingLength(), params.seedAffixTrimLength())
    , orientationPredictor_(graph, params.orientationPredictorKmerLen(), params.orientationPredictorMinKmerCount())
    , writer_(std::move(writer))
    , alignmentBuffer_(std::move(buffer))
{
}

LocusAligner::AlignedPair LocusAligner::align(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector)
{
    auto readAlign = align(read, alignerSelector);
    auto mateAlign = mate ? align(*mate, alignerSelector) : boost::none;

    int numMatchingBases = static_cast<int>(static_cast<double>(read.sequence().length()) / 7.5);
    numMatchingBases = std::max(numMatchingBases, 10);
    LinearAlignmentParameters parameters;
    const int kMinNonRepeatAlignmentScore = numMatchingBases * parameters.matchScore;

    if (!checkIfLocallyPlacedReadPair(readAlign, mateAlign, kMinNonRepeatAlignmentScore))
    {
        return { boost::none, boost::none };
    }

    if (readAlign && mateAlign)
    {
        // Optionally buffer fragments for specialized caller extensions (RFC1, REViewer):
        if (alignmentBuffer_)
        {
            AlignedFragment fragment {
                read.fragmentId(),
                read.sequence(),
                *readAlign,
                !read.isReversed(),
                mate->sequence(),
                *mateAlign,
                !mate->isReversed()
            };
            alignmentBuffer_->tryInsertFragment(std::move(fragment));
        }

        // Output realigned reads to bam:
        writer_->write(
            locusId_, read.fragmentId(), read.sequence(), read.isFirstMate(), read.isReversed(), mate->isReversed(),
            *readAlign);
        writer_->write(
            locusId_, mate->fragmentId(), mate->sequence(), mate->isFirstMate(), mate->isReversed(), read.isReversed(),
            *mateAlign);
    }

    return { readAlign, mateAlign };
}

LocusAligner::OptionalAlign LocusAligner::align(Read& read, graphtools::AlignerSelector& alignerSelector) const
{
    // Key on the original (pre-orientation) sequence. Look up by const reference first so a cache hit (or a
    // miss on an all-unique locus) costs no string copy; only copy when we need an owned key to insert,
    // before reverseComplement() below mutates read.sequence().
    auto cached = alignCache_.find(read.sequence());
    if (cached != alignCache_.end())
    {
        // Replay align()'s in-place side effect on the read before returning the memoized result.
        if (cached->second.orientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
        {
            read.reverseComplement();
        }
        return cached->second.align;
    }

    std::string sequence = read.sequence();
    OrientationPrediction predictedOrientation = orientationPredictor_.predict(sequence);

    OptionalAlign result;
    if (predictedOrientation == OrientationPrediction::kAlignsInReverseComplementOrientation)
    {
        read.reverseComplement();
    }
    if (predictedOrientation != OrientationPrediction::kDoesNotAlign)
    {
        auto readAligns = aligner_.align(read.sequence(), alignerSelector);
        if (!readAligns.empty())
        {
            result = computeCanonicalAlignment(readAligns);
        }
    }

    if (alignCache_.size() < kMaxAlignCacheEntries)
    {
        alignCache_.emplace(std::move(sequence), CachedAlign { predictedOrientation, result });
    }

    return result;
}

}
}
