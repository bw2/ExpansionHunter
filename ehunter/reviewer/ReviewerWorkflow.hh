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
#include <string>

#include "locus/AlignmentBuffer.hh"
#include "locus/LocusSpecification.hh"
#include "locus/LocusFindings.hh"
#include "reviewer/Aligns.hh"

namespace ehunter
{
namespace reviewer
{

/// Context containing intermediate data from the reviewer workflow.
/// This allows reuse of computed data for both SVG generation and metrics computation.
struct ReviewerContext
{
    GraphPaths paths;                    // Best diplotype (haplotype paths)
    FragById fragById;                   // Fragments indexed by ID
    FragPathAlignsById fragPathAlignsById;  // Fragment path alignments
    FragAssignment fragAssignment;       // Assignment of fragments to haplotypes

    ReviewerContext()
        : fragAssignment({}, {})
    {
    }

    ReviewerContext(
        GraphPaths paths_,
        FragById fragById_,
        FragPathAlignsById fragPathAlignsById_,
        FragAssignment fragAssignment_)
        : paths(std::move(paths_))
        , fragById(std::move(fragById_))
        , fragPathAlignsById(std::move(fragPathAlignsById_))
        , fragAssignment(std::move(fragAssignment_))
    {
    }
};

/// Check if read visualization should be generated for this locus based on its plot policy
/// @param locusSpec Locus specification with plot policy and conditions
/// @param findings Genotype results to check against conditions
/// @return true if visualization should be generated
bool shouldPlotReadVisualization(
    const LocusSpecification& locusSpec,
    const LocusFindings& findings);

/// Run REViewer workflow and return context for further processing.
/// This computes all intermediate data structures but does not generate SVG.
/// @param locusSpec Locus specification
/// @param alignmentBuffer Paired fragment alignments
/// @param findings Genotype results
/// @param meanFragLen Mean fragment length
/// @return ReviewerContext if successful, empty optional if workflow fails
std::optional<ReviewerContext> runReviewerWorkflow(
    const LocusSpecification& locusSpec,
    const locus::AlignmentBuffer& alignmentBuffer,
    const LocusFindings& findings,
    int meanFragLen);

/// Run REViewer workflow and generate SVG
/// @param locusSpec Locus specification
/// @param alignmentBuffer Paired fragment alignments
/// @param findings Genotype results
/// @param outputPrefix Output file prefix (SVG will be <prefix>.<locusId>.svg)
/// @param meanFragLen Mean fragment length
void runReviewer(
    const LocusSpecification& locusSpec,
    const locus::AlignmentBuffer& alignmentBuffer,
    const LocusFindings& findings,
    const std::string& outputPrefix,
    int meanFragLen);

}  // namespace reviewer
}  // namespace ehunter
