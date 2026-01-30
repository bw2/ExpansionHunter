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
//

#include "locus/LocusAnalyzer.hh"
#include "locus/AlleleQualityMetrics.hh"
#include "locus/LocusAligner.hh"
#include "locus/RFC1MotifAnalysis.hh"
#include "locus/RepeatAnalyzer.hh"
#include "locus/SmallVariantAnalyzer.hh"
#include "reviewer/GenerateSvg.hh"
#include "reviewer/LanePlot.hh"
#include "reviewer/Metrics.hh"
#include "reviewer/ReviewerWorkflow.hh"

#include "spdlog/spdlog.h"

using boost::optional;
using graphtools::GraphAlignment;
using graphtools::NodeId;
using std::string;
using std::vector;

namespace ehunter
{
namespace locus
{

LocusAnalyzer::LocusAnalyzer(LocusSpecification locusSpec, const HeuristicParameters& params, BamletWriterPtr writer,
                             bool enableAlleleQualityMetrics)
    : locusSpec_(std::move(locusSpec))
    , enableAlleleQualityMetrics_(enableAlleleQualityMetrics)
    // Create buffer if requiresAlignmentBuffer() (RFC1, plot-all policy, or has plot conditions)
    // or if allele quality metrics are enabled
    , alignmentBuffer_((locusSpec_.requiresAlignmentBuffer() || enableAlleleQualityMetrics)
                       ? std::make_shared<locus::AlignmentBuffer>() : nullptr)
    , aligner_(locusSpec_.locusId(), &locusSpec_.regionGraph(), params, std::move(writer), alignmentBuffer_)
    , statsCalc_(locusSpec_.typeOfChromLocusLocatedOn(), locusSpec_.regionGraph())
{
    if (aligner_.bamletWriter())
    {
        aligner_.bamletWriter()->initLocusSpec(locusSpec_);
    }

    for (const auto& variantSpec : locusSpec_.variantSpecs())
    {
        if (variantSpec.classification().type == VariantType::kRepeat)
        {
            const auto& graph = locusSpec_.regionGraph();
            const int repeatNodeId = static_cast<int>(variantSpec.nodes().front());
            const auto& motif = graph.nodeSeq(repeatNodeId);

            if (variantSpec.classification().subtype == VariantSubtype::kRareRepeat)
            {
                if (irrPairFinder())
                {
                    const string message
                        = "Region " + locusSpec_.locusId() + " must not have more than one rare repeat";
                    throw std::logic_error(message);
                }
                addIrrPairFinder(motif);
            }

            addRepeatAnalyzer(variantSpec.id(), repeatNodeId);
        }
        else if (variantSpec.classification().type == VariantType::kSmallVariant)
        {
            addSmallVariantAnalyzer(
                variantSpec.id(), variantSpec.classification().subtype, variantSpec.nodes(),
                variantSpec.optionalRefNode());
        }
        else
        {
            std::stringstream encoding;
            encoding << variantSpec.classification().type << "/" << variantSpec.classification().subtype;
            throw std::logic_error("Missing logic to create an analyzer for " + encoding.str());
        }
    }
}

void LocusAnalyzer::processMates(
    Read& read, Read* mate, RegionType regionType, graphtools::AlignerSelector& alignerSelector)
{
    if (regionType == RegionType::kTarget)
    {
        processOntargetMates(read, mate, alignerSelector);
    }
    else if (mate)
    {
        processOfftargetMates(read, *mate);
    }
}

void LocusAnalyzer::processOntargetMates(Read& read, Read* mate, graphtools::AlignerSelector& alignerSelector)
{
    auto alignedPair = aligner_.align(read, mate, alignerSelector);

    const bool neitherMateAligned = !alignedPair.first && !alignedPair.second;
    const bool bothMatesAligned = alignedPair.first && alignedPair.second;

    if (irrPairFinder_ && neitherMateAligned && mate)
    {
        processOfftargetMates(read, *mate);
        return;
    }

    if (bothMatesAligned)
    {
        statsCalc_.inspect(*alignedPair.first, *alignedPair.second);
        runVariantAnalysis(read, *alignedPair.first, *mate, *alignedPair.second);
    }
    else
    {
        if (alignedPair.first)
        {
            statsCalc_.inspectRead(*alignedPair.first);
        }
        if (alignedPair.second)
        {
            statsCalc_.inspectRead(*alignedPair.second);
        }
    }
}

void LocusAnalyzer::processOfftargetMates(const Read& read, const Read& mate)
{
    if (!irrPairFinder_)
    {
        const string message = "Locus " + locusSpec_.locusId() + " is not supposed to have offtarget read pairs";
        throw std::logic_error(message);
    }

    if (irrPairFinder_->check(read.sequence(), mate.sequence()))
    {
        int numAnalyzersFound = 0;
        for (auto& variantAnalyzer : variantAnalyzers_)
        {
            auto repeatAnalyzer = dynamic_cast<RepeatAnalyzer*>(variantAnalyzer.get());
            if (repeatAnalyzer != nullptr && repeatAnalyzer->repeatUnit() == irrPairFinder_->targetMotif())
            {
                numAnalyzersFound++;
                repeatAnalyzer->addInrepeatReadPair();
            }
        }

        if (numAnalyzersFound != 1)
        {
            const string message = "Locus " + locusSpec_.locusId() + " must have exactly one rare motif";
            throw std::logic_error(message);
        }
    }
}

LocusFindings LocusAnalyzer::analyze(
    Sex sampleSex, boost::optional<double> genomeWideDepth,
    const std::string& outputPrefix)
{
    LocusFindings locusFindings(statsCalc_.estimate(sampleSex));
    if (genomeWideDepth && locusSpec_.requiresGenomeWideDepth())
    {
        locusFindings.stats.setDepth(*genomeWideDepth);
    }

    for (auto& variantAnalyzer : variantAnalyzers_)
    {
        std::unique_ptr<VariantFindings> variantFindingsPtr = variantAnalyzer->analyze(locusFindings.stats);
        const string& variantId = variantAnalyzer->variantId();
        locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));
    }

    // Run RFC1 caller if required for this locus:
    if (locusSpec().useRFC1MotifAnalysis())
    {
        assert(alignmentBuffer_);
        runRFC1MotifAnalysis(*alignmentBuffer_, locusFindings);
    }

    // Run REViewer workflow once if alignment buffer exists (for metrics and/or SVG)
    std::optional<reviewer::ReviewerContext> reviewerContext;
    const bool needsReviewer = alignmentBuffer_ &&
        (reviewer::shouldPlotReadVisualization(locusSpec_, locusFindings) || (enableAlleleQualityMetrics_ && alignmentBuffer_->getBuffer().size() > 0));

    if (needsReviewer)
    {
        try
        {
            reviewerContext = reviewer::runReviewerWorkflow(
                locusSpec_, *alignmentBuffer_, locusFindings, locusFindings.stats.meanFragLength());
        }
        catch (const std::exception& e)
        {
            spdlog::debug("Failed to run reviewer workflow for locus {}: {}", locusSpec_.locusId(), e.what());
        }
    }

    // Compute and attach allele quality metrics if enabled and we have reviewer context
    if (enableAlleleQualityMetrics_ && reviewerContext)
    {
        try
        {
            auto metricsByVariant = reviewer::getMetrics(
                locusSpec_, reviewerContext->paths, reviewerContext->fragById,
                reviewerContext->fragAssignment, reviewerContext->fragPathAlignsById);

            for (const auto& metrics : metricsByVariant)
            {
                auto it = locusFindings.findingsForEachVariant.find(metrics.variantId);
                if (it != locusFindings.findingsForEachVariant.end())
                {
                    RepeatFindings* repeatFindings = dynamic_cast<RepeatFindings*>(it->second.get());
                    if (repeatFindings && repeatFindings->optionalGenotype())
                    {
                        const auto& genotype = *repeatFindings->optionalGenotype();
                        RepeatAlleleQualityMetrics alleleQualityMetrics;
                        alleleQualityMetrics.variantId = metrics.variantId;
                        alleleQualityMetrics.hasMetrics = true;

                        // Determine if homozygous (single allele metrics) or heterozygous (two)
                        bool isHomozygous = (genotype.numAlleles() == 1) ||
                            (genotype.shortAlleleSizeInUnits() == genotype.longAlleleSizeInUnits());

                        if (isHomozygous)
                        {
                            // Homozygous/hemizygous: produce one AlleleMetrics with combined depth
                            AlleleMetrics allele;
                            allele.alleleNumber = 1;
                            allele.alleleSize = genotype.shortAlleleSizeInUnits();
                            // Sum depths from all paths for homozygous case
                            double totalDepth = 0.0;
                            for (double d : metrics.alleleDepth)
                            {
                                totalDepth += d;
                            }
                            allele.depth = totalDepth;
                            // For homozygous, average QD across paths (weighted by depth if available)
                            double totalQualitySum = 0.0;
                            for (size_t i = 0; i < metrics.qd.size() && i < metrics.alleleDepth.size(); ++i)
                            {
                                totalQualitySum += metrics.qd[i] * metrics.alleleDepth[i];
                            }
                            allele.qd = (totalDepth > 0.0) ? (totalQualitySum / totalDepth) : 0.0;
                            // For homozygous, average meanInsertedBasesWithinRepeats across paths (weighted by depth)
                            double totalMeanInsertedBasesWeighted = 0.0;
                            for (size_t i = 0; i < metrics.meanInsertedBasesWithinRepeats.size() && i < metrics.alleleDepth.size(); ++i)
                            {
                                totalMeanInsertedBasesWeighted += metrics.meanInsertedBasesWithinRepeats[i] * metrics.alleleDepth[i];
                            }
                            allele.meanInsertedBasesWithinRepeats = (totalDepth > 0.0) ? (totalMeanInsertedBasesWeighted / totalDepth) : 0.0;
                            // For homozygous, average meanDeletedBasesWithinRepeats across paths (weighted by depth)
                            double totalMeanDeletedBasesWeighted = 0.0;
                            for (size_t i = 0; i < metrics.meanDeletedBasesWithinRepeats.size() && i < metrics.alleleDepth.size(); ++i)
                            {
                                totalMeanDeletedBasesWeighted += metrics.meanDeletedBasesWithinRepeats[i] * metrics.alleleDepth[i];
                            }
                            allele.meanDeletedBasesWithinRepeats = (totalDepth > 0.0) ? (totalMeanDeletedBasesWeighted / totalDepth) : 0.0;
                            // For homozygous, average the already-normalized flank depths across haplotypes
                            // Note: metrics.leftFlankNormalizedDepth already contains depth/alleleDepth per haplotype
                            double totalLeftFlankNormalized = 0.0;
                            for (size_t i = 0; i < metrics.leftFlankNormalizedDepth.size(); ++i)
                            {
                                totalLeftFlankNormalized += metrics.leftFlankNormalizedDepth[i];
                            }
                            allele.leftFlankNormalizedDepth = (!metrics.leftFlankNormalizedDepth.empty()) ? (totalLeftFlankNormalized / metrics.leftFlankNormalizedDepth.size()) : 0.0;
                            // For homozygous, average the already-normalized right flank depths across haplotypes
                            double totalRightFlankNormalized = 0.0;
                            for (size_t i = 0; i < metrics.rightFlankNormalizedDepth.size(); ++i)
                            {
                                totalRightFlankNormalized += metrics.rightFlankNormalizedDepth[i];
                            }
                            allele.rightFlankNormalizedDepth = (!metrics.rightFlankNormalizedDepth.empty()) ? (totalRightFlankNormalized / metrics.rightFlankNormalizedDepth.size()) : 0.0;
                            // Sum highQualityUnambiguousReads from all paths for homozygous case
                            int totalHighQualUnambiguous = 0;
                            for (int count : metrics.highQualityUnambiguousReads)
                            {
                                totalHighQualUnambiguous += count;
                            }
                            allele.highQualityUnambiguousReads = totalHighQualUnambiguous;
                            // For homozygous, combine raw strand counts and recompute the binomial p-value
                            // (averaging phred-scaled p-values is statistically invalid)
                            int totalForwardReads = 0;
                            int totalReverseReads = 0;
                            for (size_t i = 0; i < metrics.forwardStrandReads.size(); ++i)
                            {
                                totalForwardReads += metrics.forwardStrandReads[i];
                            }
                            for (size_t i = 0; i < metrics.reverseStrandReads.size(); ++i)
                            {
                                totalReverseReads += metrics.reverseStrandReads[i];
                            }
                            allele.strandBiasBinomialPhred = reviewer::computeStrandBiasBinomialPhred(totalForwardReads, totalReverseReads);
                            // Compute confidenceIntervalDividedByAlleleSize for homozygous case
                            const auto ci = genotype.shortAlleleSizeInUnitsCi();
                            if (allele.alleleSize > 0)
                            {
                                allele.confidenceIntervalDividedByAlleleSize = static_cast<double>(ci.end() - ci.start()) / allele.alleleSize;
                            }
                            alleleQualityMetrics.alleles.push_back(allele);
                        }
                        else
                        {
                            // Heterozygous: produce two AlleleMetrics
                            // metrics.genotype[i] is repeat units for haplotype i
                            // metrics.alleleDepth[i] is depth for haplotype i
                            for (size_t i = 0; i < metrics.genotype.size() && i < metrics.alleleDepth.size(); ++i)
                            {
                                AlleleMetrics allele;
                                allele.alleleNumber = static_cast<int>(i) + 1;
                                allele.alleleSize = metrics.genotype[i];
                                allele.depth = metrics.alleleDepth[i];
                                // Copy QD for this haplotype
                                if (i < metrics.qd.size())
                                {
                                    allele.qd = metrics.qd[i];
                                }
                                // Copy meanInsertedBasesWithinRepeats for this haplotype
                                if (i < metrics.meanInsertedBasesWithinRepeats.size())
                                {
                                    allele.meanInsertedBasesWithinRepeats = metrics.meanInsertedBasesWithinRepeats[i];
                                }
                                // Copy meanDeletedBasesWithinRepeats for this haplotype
                                if (i < metrics.meanDeletedBasesWithinRepeats.size())
                                {
                                    allele.meanDeletedBasesWithinRepeats = metrics.meanDeletedBasesWithinRepeats[i];
                                }
                                // Copy already-normalized left flank depth for this haplotype
                                if (i < metrics.leftFlankNormalizedDepth.size())
                                {
                                    allele.leftFlankNormalizedDepth = metrics.leftFlankNormalizedDepth[i];
                                }
                                // Copy already-normalized right flank depth for this haplotype
                                if (i < metrics.rightFlankNormalizedDepth.size())
                                {
                                    allele.rightFlankNormalizedDepth = metrics.rightFlankNormalizedDepth[i];
                                }
                                // Copy highQualityUnambiguousReads for this haplotype
                                if (i < metrics.highQualityUnambiguousReads.size())
                                {
                                    allele.highQualityUnambiguousReads = metrics.highQualityUnambiguousReads[i];
                                }
                                // Copy strandBiasBinomialPhred for this haplotype
                                if (i < metrics.strandBiasBinomialPhred.size())
                                {
                                    allele.strandBiasBinomialPhred = metrics.strandBiasBinomialPhred[i];
                                }
                                // Compute confidenceIntervalDividedByAlleleSize for heterozygous case
                                // Match CI by comparing actual allele size to genotype short/long allele sizes
                                // (don't rely on index ordering which depends on reviewer's diplotype path ordering)
                                const bool isShortAllele = (allele.alleleSize == genotype.shortAlleleSizeInUnits());
                                const auto ci = isShortAllele ? genotype.shortAlleleSizeInUnitsCi() : genotype.longAlleleSizeInUnitsCi();
                                if (allele.alleleSize > 0)
                                {
                                    allele.confidenceIntervalDividedByAlleleSize = static_cast<double>(ci.end() - ci.start()) / allele.alleleSize;
                                }
                                alleleQualityMetrics.alleles.push_back(allele);
                            }
                        }

                        repeatFindings->setAlleleQualityMetrics(alleleQualityMetrics);
                        repeatFindings->setCountsOfHighQualityUnambiguousReads(metrics.countsOfHighQualityUnambiguousReads);
                    }
                }
            }
        }
        catch (const std::exception& e)
        {
            spdlog::error("Failed to compute allele quality metrics for locus {}: {}", locusSpec_.locusId(), e.what());
        }
    }

    // Propagate consensus sequences and statistics from ReviewerContext to RepeatFindings
    // NOTE: Consensus building currently only supports loci with a single repeat variant.
    // For loci with multiple repeat variants, the consensus would be a mix of all variants,
    // which is not meaningful to attach to individual variant findings.
    if (reviewerContext && !reviewerContext->consensusResult.alleleConsensuses.empty())
    {
        // Count repeat variants and find the single repeat variant ID if there's only one
        int repeatVariantCount = 0;
        std::string singleRepeatVariantId;
        for (const auto& variantSpec : locusSpec_.variantSpecs())
        {
            if (variantSpec.classification().type == VariantType::kRepeat)
            {
                repeatVariantCount++;
                singleRepeatVariantId = variantSpec.id();
            }
        }

        if (repeatVariantCount > 1)
        {
            spdlog::warn("Skipping consensus for locus {} - multiple repeat variants not supported",
                         locusSpec_.locusId());
        }
        else if (repeatVariantCount == 1)
        {
            // Extract consensus strings and read support from the consensus result
            std::vector<std::string> consensusStrings;
            std::vector<std::string> consensusReadSupport;
            for (const auto& alleleConsensus : reviewerContext->consensusResult.alleleConsensuses)
            {
                consensusStrings.push_back(alleleConsensus.toString());
                consensusReadSupport.push_back(alleleConsensus.toReadSupportString());
            }

            // Attach consensus sequences and read support only to the single repeat variant
            auto it = locusFindings.findingsForEachVariant.find(singleRepeatVariantId);
            if (it != locusFindings.findingsForEachVariant.end())
            {
                RepeatFindings* repeatFindings = dynamic_cast<RepeatFindings*>(it->second.get());
                if (repeatFindings)
                {
                    repeatFindings->setConsensusSequences(consensusStrings);
                    repeatFindings->setConsensusReadSupport(consensusReadSupport);
                }
            }

            spdlog::debug("Consensus sequences for locus {}: {} alleles, {} total anchors",
                          locusSpec_.locusId(),
                          consensusStrings.size(),
                          reviewerContext->consensusResult.totalAnchors);
        }
    }

    // Generate SVG if visualization is requested (controlled separately from quality metrics)
    if (reviewerContext && reviewer::shouldPlotReadVisualization(locusSpec_, locusFindings))
    {
        try
        {
            // Use the context we already have to generate the SVG
            auto lanePlots = reviewer::generateBlueprint(
                reviewerContext->paths, reviewerContext->fragById,
                reviewerContext->fragAssignment, reviewerContext->fragPathAlignsById);
            const std::string svgPath = outputPrefix + "." + locusSpec_.locusId() + ".svg";
            reviewer::generateSvg(lanePlots, svgPath);
            spdlog::info("REViewer: Generated SVG for locus {} at {}", locusSpec_.locusId(), svgPath);
        }
        catch (const std::exception& e)
        {
            spdlog::error("Failed to generate image for locus {}: {}", locusSpec_.locusId(), e.what());
        }
    }

    return locusFindings;
}

void LocusAnalyzer::addIrrPairFinder(std::string motif) { irrPairFinder_ = IrrPairFinder(std::move(motif)); }

void LocusAnalyzer::addRepeatAnalyzer(std::string variantId, graphtools::NodeId nodeId)
{
    variantAnalyzers_.emplace_back(std::make_unique<RepeatAnalyzer>(
        std::move(variantId), locusSpec_.regionGraph(), nodeId, locusSpec_.genotyperParameters()));
}

void LocusAnalyzer::addSmallVariantAnalyzer(
    string variantId, VariantSubtype subtype, vector<NodeId> nodes, optional<NodeId> refNode)
{
    variantAnalyzers_.emplace_back(std::make_unique<SmallVariantAnalyzer>(
        std::move(variantId), subtype, locusSpec_.regionGraph(), std::move(nodes), refNode,
        locusSpec_.genotyperParameters()));
}

void LocusAnalyzer::runVariantAnalysis(
    const Read& read, const LocusAnalyzer::Align& readAlign, const Read& mate, const LocusAnalyzer::Align& mateAlign)
{
    for (auto& analyzer : variantAnalyzers_)
    {
        analyzer->processMates(read, readAlign, mate, mateAlign);
    }
}

}
}
