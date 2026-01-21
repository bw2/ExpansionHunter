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

#include "reviewer/ReviewerWorkflow.hh"

#include <cstdlib>
#include <optional>
#include <stdexcept>

#include "spdlog/spdlog.h"

#include "reviewer/Aligns.hh"
#include "reviewer/FragLenFilter.hh"
#include "reviewer/GenerateSvg.hh"
#include "reviewer/GenotypePaths.hh"
#include "reviewer/LanePlot.hh"
#include "reviewer/Origin.hh"
#include "reviewer/Phasing.hh"
#include "reviewer/Projection.hh"

namespace ehunter
{
namespace reviewer
{

/// Evaluate a single plot condition against the genotype
static bool evaluateCondition(const PlotReadVisualization& cond, const LocusFindings& findings)
{
    // Get the first repeat variant findings
    int shortAllele = 0;
    int longAllele = 0;

    for (const auto& [variantId, variantFindings] : findings.findingsForEachVariant)
    {
        const RepeatFindings* repeatFindings = dynamic_cast<const RepeatFindings*>(variantFindings.get());
        if (repeatFindings && repeatFindings->optionalGenotype())
        {
            shortAllele = repeatFindings->optionalGenotype()->shortAlleleSizeInUnits();
            longAllele = repeatFindings->optionalGenotype()->longAlleleSizeInUnits();
            break;
        }
    }

    int alleleValue = (cond.appliedTo == PlotThresholdAppliedTo::kShortAllele) ? shortAllele : longAllele;

    switch (cond.op)
    {
    case PlotThresholdComparisonOp::kLessThan:
        return alleleValue < cond.threshold;
    case PlotThresholdComparisonOp::kLessThanOrEqual:
        return alleleValue <= cond.threshold;
    case PlotThresholdComparisonOp::kEqual:
        return alleleValue == cond.threshold;
    case PlotThresholdComparisonOp::kNotEqual:
        return alleleValue != cond.threshold;
    case PlotThresholdComparisonOp::kGreaterThanOrEqual:
        return alleleValue >= cond.threshold;
    case PlotThresholdComparisonOp::kGreaterThan:
        return alleleValue > cond.threshold;
    }
    return false;
}

/// Check if any plot condition is met
static bool anyConditionMet(const std::vector<PlotReadVisualization>& conditions, const LocusFindings& findings)
{
    for (const auto& cond : conditions)
    {
        if (evaluateCondition(cond, findings))
        {
            return true;
        }
    }
    return false;
}

bool shouldPlotReadVisualization(
    const LocusSpecification& locusSpec,
    const LocusFindings& findings)
{
    switch (locusSpec.plotPolicy())
    {
    case PlotPolicy::kNone:
        return false;
    case PlotPolicy::kAll:
        return true;
    case PlotPolicy::kConditional:
        if (locusSpec.hasPlotConditions())
        {
            return anyConditionMet(locusSpec.plotConditions(), findings);
        }
        return false;
    }
    return false;
}

/// Convert ExpansionHunter AlignmentBuffer to FragById format using ehunter::Read
static FragById convertToFragById(const locus::AlignmentBuffer& buffer)
{
    FragById fragById;
    for (const auto& [fragId, alignedFrag] : buffer.getBuffer())
    {
        // Create ReadId for the read (mate 1)
        ReadId readId(fragId, MateNumber::kFirstMate);
        // isReversed is the opposite of isForwardStrand
        bool readIsReversed = !alignedFrag.readIsForwardStrand;
        ehunter::Read ehRead(readId, alignedFrag.readBases, readIsReversed);
        ReadWithAlign read(std::move(ehRead), alignedFrag.readAlignment);

        // Create ReadId for the mate (mate 2)
        ReadId mateId(fragId, MateNumber::kSecondMate);
        bool mateIsReversed = !alignedFrag.mateIsForwardStrand;
        ehunter::Read ehMate(mateId, alignedFrag.mateBases, mateIsReversed);
        ReadWithAlign mate(std::move(ehMate), alignedFrag.mateAlignment);

        fragById.emplace(fragId, Frag(std::move(read), std::move(mate)));
    }
    return fragById;
}

std::optional<ReviewerContext> runReviewerWorkflow(
    const LocusSpecification& locusSpec,
    const locus::AlignmentBuffer& alignmentBuffer,
    const LocusFindings& findings,
    int meanFragLen)
{
    const std::string& locusId = locusSpec.locusId();

    spdlog::debug("REViewer workflow: Processing locus {}", locusId);

    // For reproducibility - same seed as standalone REViewer
    srand(14345);

    // Convert alignment buffer to REViewer format
    auto fragById = convertToFragById(alignmentBuffer);
    if (fragById.empty())
    {
        spdlog::warn("REViewer workflow: No fragments available for locus {}", locusId);
        return std::nullopt;
    }
    spdlog::debug("REViewer workflow: Converted {} fragments", fragById.size());

    // Use provided mean fragment length, or calculate from flanking reads
    int effectiveMeanFragLen = meanFragLen;
    if (effectiveMeanFragLen <= 0)
    {
        try
        {
            effectiveMeanFragLen = getMeanFragLen(fragById);
        }
        catch (const std::exception& e)
        {
            spdlog::warn("REViewer workflow: Unable to determine fragment length for {}: {}", locusId, e.what());
            effectiveMeanFragLen = 350;  // Default fallback
        }
    }
    spdlog::debug("REViewer workflow: Using mean fragment length {}", effectiveMeanFragLen);

    // Generate candidate diplotypes from genotype findings (not VCF)
    auto candidateDiplotypes = getCandidateDiplotypes(effectiveMeanFragLen, locusSpec, findings);
    if (candidateDiplotypes.empty())
    {
        spdlog::warn("REViewer workflow: No candidate diplotypes generated for {}", locusId);
        return std::nullopt;
    }
    spdlog::debug("REViewer workflow: Generated {} candidate diplotypes", candidateDiplotypes.size());

    // Score diplotypes by alignment support
    auto scoredDiplotypes = scoreDiplotypes(fragById, candidateDiplotypes);
    auto topDiplotype = scoredDiplotypes.front().first;
    spdlog::debug("REViewer workflow: Top diplotype has {} haplotype paths", topDiplotype.size());

    // Project reads onto best diplotype
    auto pairPathAlignById = project(topDiplotype, fragById);
    spdlog::debug("REViewer workflow: Projected {} read pairs", pairPathAlignById.size());

    // Filter by fragment length
    auto fragPathAlignsById = resolveByFragLen(effectiveMeanFragLen, topDiplotype, pairPathAlignById);
    spdlog::debug("REViewer workflow: Resolved {} fragment alignments", fragPathAlignsById.size());

    // Assign fragment origins
    auto fragAssignment = getBestFragAssignment(topDiplotype, fragPathAlignsById);
    spdlog::debug("REViewer workflow: Assigned {} fragments", fragAssignment.fragIds.size());

    return ReviewerContext(
        std::move(topDiplotype),
        std::move(fragById),
        std::move(fragPathAlignsById),
        std::move(fragAssignment));
}

void runReviewer(
    const LocusSpecification& locusSpec,
    const locus::AlignmentBuffer& alignmentBuffer,
    const LocusFindings& findings,
    const std::string& outputPrefix,
    int meanFragLen)
{
    const std::string& locusId = locusSpec.locusId();

    spdlog::info("REViewer: Processing locus {}", locusId);

    // Run the workflow to get context
    auto contextOpt = runReviewerWorkflow(locusSpec, alignmentBuffer, findings, meanFragLen);
    if (!contextOpt)
    {
        spdlog::warn("REViewer: Workflow failed for locus {}", locusId);
        return;
    }

    const auto& context = *contextOpt;

    // Generate visualization blueprint
    auto lanePlots = generateBlueprint(context.paths, context.fragById, context.fragAssignment, context.fragPathAlignsById);

    // Generate SVG output
    const std::string svgPath = outputPrefix + "." + locusId + ".svg";
    generateSvg(lanePlots, svgPath);

    spdlog::info("REViewer: Generated SVG for locus {} at {}", locusId, svgPath);
}

}  // namespace reviewer
}  // namespace ehunter
