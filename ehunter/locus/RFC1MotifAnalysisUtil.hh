//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Chris Saunders <csaunders@illumina.com>
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

/// \file
///
/// \brief Items used by the RFC1 motif analyzer. They are separated here to enable unit testing without adding them
/// to the interface.

#pragma once

#include <numeric>
#include <string>
#include <vector>

#include "boost/optional.hpp"

#include "graphalign/GraphAlignment.hh"
#include "locus/AlignmentBuffer.hh"

namespace ehunter
{

/// \brief Single-read alignment data for RFC1 motif analysis
///
/// This structure preserves the original RFC1 interface which operates on individual reads.
/// It is populated from the new paired-fragment AlignmentBuffer via extractRFC1ReadAlignments().
///
struct RFC1ReadAlignment
{
    std::string read;
    bool isReversed;
    graphtools::GraphAlignment readAlignment;
};

/// \brief Extract read alignments from paired-fragment buffer for RFC1 analysis
///
/// RFC1 motif analysis operates on individual reads, not paired fragments.
/// This function converts the paired-fragment AlignmentBuffer to a vector of
/// single-read alignments suitable for RFC1 analysis.
///
/// Only the read (not mate) from each fragment is extracted, matching the original
/// RFC1 behavior where only the read side was buffered.
///
/// \param alignmentBuffer The paired-fragment alignment buffer
/// \return Vector of read alignments that overlap the repeat region
///
std::vector<RFC1ReadAlignment> extractRFC1ReadAlignments(const locus::AlignmentBuffer& alignmentBuffer);

/// \brief Division with divide by zero guard
///
/// This should consistently do the sane thing for all integral and floating point types, except for floating point
/// types wider than double..
///
template <typename A, typename B> double safeFrac(const A a, const B b)
{
    const auto bd(static_cast<double>(b));
    return (((bd <= 0.) && (bd >= 0)) ? 0. : (a / bd));
}

/// \brief Return mean of the elements in an iterator range
///
template <typename T> double mean(T begin, T end)
{
    return (std::accumulate(begin, end, 0.0) / std::distance(begin, end));
}

/// \brief Return the lexicographical minimum rotation from all rotations of \p str
///
std::string getMinRotation(std::string str);

/// \breif Determine the range of bases in a read which are usable for repeat motif extraction.
///
/// This routine will trim off the 3' end of the read at a defined distance before the first low-quality base in the
/// read
///
/// \param[in] binaryQuals Quality vector for the read, reduced to 2 {low, high} quality states
///
/// \param[in] isReversed True if the read is aligned in reverse orientation
///
/// \return 2-tuple of start and end positions representing a zero-indexed, closed interval of usable base positions in
/// read coordinates. none is returned when no usable bases are found.
///
boost::optional<std::pair<unsigned, unsigned>>
findUsableReadBaseRange(std::vector<uint8_t> binaryQuals, bool isReversed);

}
