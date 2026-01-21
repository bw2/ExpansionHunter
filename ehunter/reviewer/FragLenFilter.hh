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
#include <vector>

#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"
#include "reviewer/Projection.hh"

namespace ehunter
{
namespace reviewer
{

/// Calculate mean fragment length from flanking reads
/// @param fragById Map of fragments by ID
/// @return Mean fragment length
int getMeanFragLen(const FragById& fragById);

/// Filter fragment alignments by fragment length consistency
/// @param meanFragLen Mean fragment length
/// @param paths Diplotype paths
/// @param pairPathAlignById Pair path alignments
/// @return Filtered fragment path alignments
FragPathAlignsById resolveByFragLen(
    int meanFragLen,
    const Diplotype& paths,
    const PairPathAlignById& pairPathAlignById);

}  // namespace reviewer
}  // namespace ehunter
