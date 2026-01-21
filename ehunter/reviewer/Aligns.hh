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
// Adapted from REViewer (Copyright 2020 Illumina, Inc., GPL-3.0 license)
//

#pragma once

#include <cassert>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Path.hh"

#include "core/Read.hh"

namespace ehunter
{
namespace reviewer
{

using GraphPath = graphtools::Path;
using GraphPaths = std::vector<GraphPath>;

using GraphAlign = graphtools::GraphAlignment;
using GraphAlignPtr = std::shared_ptr<GraphAlign>;

/// A read with its graph alignment, wrapping the core ehunter::Read
struct ReadWithAlign
{
    ReadWithAlign(ehunter::Read read, GraphAlign align)
        : read(std::move(read))
        , align(std::move(align))
    {
    }

    /// Convenience accessor for the sequence
    const std::string& bases() const { return read.sequence(); }

    ehunter::Read read;
    GraphAlign align;
};

/// A paired-end fragment containing both reads (using ReadWithAlign wrapper)
struct Frag
{
    Frag(ReadWithAlign read, ReadWithAlign mate)
        : read(std::move(read))
        , mate(std::move(mate))
    {
    }

    ReadWithAlign read;
    ReadWithAlign mate;
};

/// Map of fragments indexed by fragment ID
using FragById = std::map<std::string, Frag>;

/// Assignment of fragments to alignment indices for visualization
struct FragAssignment
{
    FragAssignment(std::vector<std::string> fragIds, std::vector<int> alignIndexByFrag)
        : fragIds(std::move(fragIds))
        , alignIndexByFrag(std::move(alignIndexByFrag))
    {
    }

    std::vector<std::string> fragIds;
    std::vector<int> alignIndexByFrag;
};

/// Alignment of a single read projected onto a haplotype path
struct ReadPathAlign
{
    ReadPathAlign(const graphtools::Path& hapPath, int pathIndex, int startIndexOnPath, GraphAlignPtr align);

    int pathIndex;       // Which haplotype (0 or 1)
    int startIndexOnPath;
    int begin;           // Position on path
    int end;             // Position on path
    GraphAlignPtr align;
};

/// Alignment of a fragment (both reads) to a haplotype path
struct FragPathAlign
{
    FragPathAlign(ReadPathAlign readAlign, ReadPathAlign mateAlign)
        : readAlign(std::move(readAlign))
        , mateAlign(std::move(mateAlign))
    {
        assert(this->readAlign.pathIndex == this->mateAlign.pathIndex);
    }

    ReadPathAlign readAlign;
    ReadPathAlign mateAlign;
};

/// Map of fragment path alignments indexed by fragment ID
using FragPathAlignsById = std::map<std::string, std::vector<FragPathAlign>>;

}  // namespace reviewer
}  // namespace ehunter
