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
#include <string>
#include <vector>

#include "graphcore/Path.hh"

#include "locus/LocusSpecification.hh"
#include "locus/LocusFindings.hh"

namespace ehunter
{
namespace reviewer
{

using Diplotype = std::vector<graphtools::Path>;

std::ostream& operator<<(std::ostream& out, const Diplotype& diplotype);

/// Computes all possible diplotype paths at the given locus
/// @param meanFragLen Mean fragment length
/// @param locusSpec Locus specification
/// @param findings Locus findings containing genotypes
/// @return Vector of all possible diplotype paths
///
/// Assumption:
/// All haplotype paths start at the first base of left flank
/// (node 0) and end at the last base of the right flank
/// (last node)
///
std::vector<Diplotype> getCandidateDiplotypes(
    int meanFragLen,
    const LocusSpecification& locusSpec,
    const LocusFindings& findings);

}  // namespace reviewer
}  // namespace ehunter
