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

#include <utility>
#include <vector>

#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

namespace ehunter
{
namespace reviewer
{

using ScoredDiplotype = std::pair<Diplotype, int>;
using ScoredDiplotypes = std::vector<ScoredDiplotype>;

/// Score diplotypes by read alignment support
/// @param fragById Map of fragments by ID
/// @param diplotypes Candidate diplotypes to score
/// @return Scored diplotypes sorted by score (highest first)
ScoredDiplotypes scoreDiplotypes(const FragById& fragById, const std::vector<Diplotype>& diplotypes);

}  // namespace reviewer
}  // namespace ehunter
