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

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Path.hh"

#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

namespace ehunter
{
namespace reviewer
{

/// Pair path alignment for a single fragment
struct PairPathAlign
{
    PairPathAlign() = default;
    PairPathAlign(std::vector<ReadPathAlign> readAligns, std::vector<ReadPathAlign> mateAligns)
        : readAligns(std::move(readAligns))
        , mateAligns(std::move(mateAligns))
    {
    }

    std::vector<ReadPathAlign> readAligns;
    std::vector<ReadPathAlign> mateAligns;
};

/// Calculate alignment score
int score(const graphtools::GraphAlignment& alignment, int matchScore = 5, int mismatchScore = -4, int gapScore = -8);

/// Map of pair path alignments indexed by fragment ID
using PairPathAlignById = std::map<std::string, PairPathAlign>;

/// Project fragments onto genotype paths
PairPathAlignById project(const std::vector<graphtools::Path>& genotypePaths, const FragById& fragById);

}  // namespace reviewer
}  // namespace ehunter
