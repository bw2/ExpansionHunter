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

#include <optional>
#include <utility>
#include <vector>

#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

namespace ehunter
{
namespace reviewer
{

using ScoredDiplotype = std::pair<Diplotype, int>;
using ScoredDiplotypes = std::vector<ScoredDiplotype>;

/// How strongly the fragments favor the top-scoring candidate diplotype over the next-best one
struct DiplotypeChoiceSupport
{
    int numberOfCandidateDiplotypes = 0;
    // Fragments whose best pair alignment score is higher on the top diplotype than on the next-best one
    int fragmentsFavoringTopDiplotype = 0;
    // Fragments whose best pair alignment score is higher on the next-best diplotype than on the top one
    int fragmentsFavoringNextBestDiplotype = 0;
    // True when the top two diplotypes have equal scores, so which one ranks first is arbitrary
    bool topDiplotypeIsTiedWithNextBest = false;
};

/// Score diplotypes by read alignment support
/// @param fragById Map of fragments by ID
/// @param diplotypes Candidate diplotypes to score
/// @return Scored diplotypes sorted by score (highest first); diplotypes with equal scores keep their input order
ScoredDiplotypes scoreDiplotypes(const FragById& fragById, const std::vector<Diplotype>& diplotypes);

/// Count the fragments that favor the top-scoring diplotype over the next-best one, and vice versa
/// @param fragById Map of fragments by ID
/// @param scoredDiplotypes Output of scoreDiplotypes (sorted highest score first)
/// @return Support summary; the fragment counts stay 0 when there are fewer than two candidates
DiplotypeChoiceSupport
compareTopTwoDiplotypes(const FragById& fragById, const ScoredDiplotypes& scoredDiplotypes);

/// Describe which allele of each repeat variant the chosen diplotype places on each haplotype
/// @param locusSpec Locus specification
/// @param findings Locus findings with the genotype of each repeat variant
/// @param chosenDiplotype Top-scoring diplotype
/// @param support Output of compareTopTwoDiplotypes for the same candidates
/// @return Phasing summary, or empty if a variant has no genotype or if its two alleles were capped to the
///         same path length (see capLengths in GenotypePaths.cpp), which hides which haplotype carries which
std::optional<RepeatAllelePhasing> summarizeRepeatAllelePhasing(
    const LocusSpecification& locusSpec, const LocusFindings& findings, const Diplotype& chosenDiplotype,
    const DiplotypeChoiceSupport& support);

}  // namespace reviewer
}  // namespace ehunter
