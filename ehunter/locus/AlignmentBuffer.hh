//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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

#pragma once

#include <map>
#include <string>

#include "graphalign/GraphAlignment.hh"

namespace ehunter
{
namespace locus
{

/// \brief A paired-end fragment with both read and mate alignments
///
/// Stores both reads of a fragment along with their graph alignments.
/// Used by REViewer for SVG visualization and by RFC1 motif analysis.
///
struct AlignedFragment
{
    std::string fragmentId;

    // Read 1
    std::string readBases;
    graphtools::GraphAlignment readAlignment;
    bool readIsForwardStrand;

    // Read 2 (mate)
    std::string mateBases;
    graphtools::GraphAlignment mateAlignment;
    bool mateIsForwardStrand;
};

/// \brief Buffer for paired-end fragment alignments at a single locus
///
/// Stores fragments (read pairs) that have both mates successfully aligned to the locus graph.
/// This buffer supports:
/// - REViewer SVG visualization (needs paired data for fragment length and phasing)
/// - RFC1 motif analysis (uses read data only, via adapter)
///
/// Fragments are stored in a map keyed by fragment ID to enable efficient lookup and
/// to ensure each fragment is stored only once.
///
class AlignmentBuffer
{
public:
    using buffer_t = std::map<std::string, AlignedFragment>;

    /// Insert a fragment into the buffer if at least one mate overlaps the repeat region
    ///
    /// \param fragment The aligned fragment to insert
    /// \return true if fragment was inserted, false if it was filtered out
    bool tryInsertFragment(AlignedFragment fragment);

    const buffer_t& getBuffer() const { return bufData_; }

    size_t size() const { return bufData_.size(); }

    bool empty() const { return bufData_.empty(); }

private:
    buffer_t bufData_;
};

}
}
