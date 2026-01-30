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

#include "reviewer/Metrics.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "graphalign/Operation.hh"
#include "htslib/kfunc.h"

namespace ehunter
{
namespace reviewer
{

using graphtools::NodeId;
using graphtools::OperationType;
using std::map;
using std::string;
using std::vector;

// Per-haplotype accumulators for single-pass metrics computation
struct HaplotypeAccumulators
{
    // For flank depth (20bp windows)
    int leftFlankCoveredBases = 0;
    int rightFlankCoveredBases = 0;

    // Alignments for depth calculation
    vector<GraphAlignPtr> alignments;
};

// Per-variant-per-haplotype accumulators (only for reads overlapping the repeat)
struct VariantHaplotypeAccumulators
{
    int numReadsOverlappingRepeat = 0;
    int totalInsertedBases = 0;
    int totalDeletedBases = 0;
    int highQualUnambiguousCount = 0;
    double qualitySum = 0.0;        // For QD calculation
    int forwardStrandCount = 0;     // For strand bias
    int reverseStrandCount = 0;     // For strand bias
};

// All accumulated data from single pass
struct AllAccumulators
{
    // Per-haplotype data
    vector<HaplotypeAccumulators> perHaplotype;

    // Per-variant-per-haplotype data (variant ID -> per-haplotype accumulators)
    map<string, vector<VariantHaplotypeAccumulators>> perVariantHaplotype;
};

double computeReadQuality(const graphtools::GraphAlignment& graphAlign)
{
    uint32_t matchedBases = 0;
    uint32_t totalBases = graphAlign.queryLength();

    if (totalBases == 0)
    {
        return 0.0;
    }

    for (const auto& nodeAlign : graphAlign.alignments())
    {
        for (const auto& op : nodeAlign)
        {
            if (op.type() == OperationType::kMatch)
            {
                matchedBases += op.queryLength();
            }
        }
    }

    return static_cast<double>(matchedBases) / static_cast<double>(totalBases);
}

// Results for repeat insertion metrics for a single alignment
struct RepeatInsertionResult
{
    int insertedBases = 0;
};

// Count insertions in the repeat node for a single graph alignment
// Returns {hasInsertion, insertedBases} for all occurrences of repeatNodeId in the alignment
static RepeatInsertionResult countRepeatNodeInsertions(
    const graphtools::GraphAlignment& graphAlign,
    NodeId repeatNodeId)
{
    RepeatInsertionResult result;

    for (size_t nodeIndex = 0; nodeIndex < graphAlign.path().numNodes(); ++nodeIndex)
    {
        if (graphAlign.getNodeIdByIndex(nodeIndex) != repeatNodeId)
        {
            continue;
        }

        // Get the alignment for this node occurrence
        const auto& nodeAlign = graphAlign.alignments()[nodeIndex];

        // Count inserted bases in this node
        for (const auto& op : nodeAlign)
        {
            if (op.type() == OperationType::kInsertionToRef)
            {
                result.insertedBases += static_cast<int>(op.queryLength());
            }
        }
    }

    return result;
}

// Compute QD (quality by depth) for all variants and haplotypes
// QD = sum(readQuality) / alleleDepth. Returns 0.0 if alleleDepth is zero.
static map<string, vector<double>> computeQD(
    const LocusSpecification& locusSpec,
    const map<string, vector<double>>& qualitySums,
    const map<string, vector<double>>& genotypeDepths)
{
    map<string, vector<double>> qd;

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        const string& variantId = variantSpec.id();
        const auto& sums = qualitySums.at(variantId);
        const auto& depths = genotypeDepths.at(variantId);

        vector<double> variantQD;
        for (size_t hapIndex = 0; hapIndex < sums.size(); ++hapIndex)
        {
            if (depths[hapIndex] > 0.0)
            {
                variantQD.push_back(sums[hapIndex] / depths[hapIndex]);
            }
            else
            {
                variantQD.push_back(0.0);
            }
        }
        qd[variantId] = variantQD;
    }

    return qd;
}

// Results for repeat deletion metrics for a single alignment
struct RepeatDeletionResult
{
    int deletedBases = 0;
};

// Count deletions in the repeat node for a single graph alignment
// Returns {hasDeletion, deletedBases} for all occurrences of repeatNodeId in the alignment
static RepeatDeletionResult countRepeatNodeDeletions(
    const graphtools::GraphAlignment& graphAlign,
    NodeId repeatNodeId)
{
    RepeatDeletionResult result;

    for (size_t nodeIndex = 0; nodeIndex < graphAlign.path().numNodes(); ++nodeIndex)
    {
        if (graphAlign.getNodeIdByIndex(nodeIndex) != repeatNodeId)
        {
            continue;
        }

        // Get the alignment for this node occurrence
        const auto& nodeAlign = graphAlign.alignments()[nodeIndex];

        // Use cached count - O(1)
        int deletedInNode = static_cast<int>(nodeAlign.numDeleted());
        result.deletedBases += deletedInNode;
    }

    return result;
}

// Compute binomial test p-value and convert to Phred scale
// Tests whether the forward/reverse read counts deviate significantly from 50/50
// Two-tailed binomial test: given n total reads and k forward reads,
// compute P(X <= min(k, n-k)) * 2 under Binomial(n, 0.5)
//
// Uses htslib's kf_betai(a, b, x) which computes the regularized incomplete beta function I_x(a,b).
// The cumulative binomial P(X <= k) for Binomial(n, p) can be computed as kf_betai(n-k, k+1, 1-p).
//
// Returns -10 * log10(p_value) (phred-scaled). Returns 0.0 if no reads.
double computeStrandBiasBinomialPhred(int forwardReads, int reverseReads)
{
    int n = forwardReads + reverseReads;

    // Edge case: no reads, no bias measurable
    if (n == 0)
    {
        return 0.0;
    }

    // For a two-tailed test with p=0.5, we compute P(X <= min(k, n-k)) * 2
    // where k is the number of forward reads
    int k = std::min(forwardReads, reverseReads);

    // P(X <= k) for Binomial(n, 0.5) = kf_betai(n-k, k+1, 1-0.5) = kf_betai(n-k, k+1, 0.5)
    double cumulativeProb = kf_betai(n - k, k + 1, 0.5);

    // Two-tailed p-value (capped at 1.0)
    double pValue = std::min(1.0, 2.0 * cumulativeProb);

    // Cap at 0 if p >= 1 (no bias)
    if (pValue >= 1.0)
    {
        return 0.0;
    }

    // Avoid log(0) - cap at very high value for very small p-values
    constexpr double MIN_P = 1e-300;
    double p = std::max(pValue, MIN_P);
    return -10.0 * std::log10(p);
}

// Compute the number of reference bases covered by the alignment within [rangeStart, rangeEnd)
// Returns overlap length - mismatches count as coverage
static int computeCoveredBasesInRange(const graphtools::Alignment& align, int rangeStart, int rangeEnd)
{
    const int alignEnd = static_cast<int>(align.referenceLength());
    const int overlapStart = std::max(0, rangeStart);
    const int overlapEnd = std::min(alignEnd, rangeEnd);
    return std::max(0, overlapEnd - overlapStart);
}

// Accumulate covered bases in the 20bp windows adjacent to the repeat for a single graph alignment
// leftWindowStart/End define [flankLen-20, flankLen) region in node-local coordinates
// rightWindowStart/End define [0, 20) region in node-local coordinates
static void accumulateFlankDepth(
    const graphtools::GraphAlignment& graphAlign,
    NodeId leftFlankNode,
    NodeId rightFlankNode,
    int leftWindowStart,
    int leftWindowEnd,
    int rightWindowStart,
    int rightWindowEnd,
    int& leftCoveredBases,
    int& rightCoveredBases)
{
    for (size_t nodeIndex = 0; nodeIndex < graphAlign.path().numNodes(); ++nodeIndex)
    {
        const auto nodeId = graphAlign.getNodeIdByIndex(nodeIndex);
        const auto& nodeAlign = graphAlign.alignments()[nodeIndex];

        if (nodeId == leftFlankNode)
        {
            // Compute overlap between alignment and left window [leftWindowStart, leftWindowEnd)
            const int alignStart = static_cast<int>(nodeAlign.referenceStart());
            const int alignEnd = alignStart + static_cast<int>(nodeAlign.referenceLength());

            const int overlapStart = std::max(alignStart, leftWindowStart);
            const int overlapEnd = std::min(alignEnd, leftWindowEnd);

            if (overlapEnd > overlapStart)
            {
                // Count covered bases within the overlap region
                // Convert to alignment-local coordinates: [overlapStart - alignStart, overlapEnd - alignStart)
                leftCoveredBases += computeCoveredBasesInRange(
                    nodeAlign, overlapStart - alignStart, overlapEnd - alignStart);
            }
        }
        else if (nodeId == rightFlankNode)
        {
            // Compute overlap between alignment and right window [rightWindowStart, rightWindowEnd)
            const int alignStart = static_cast<int>(nodeAlign.referenceStart());
            const int alignEnd = alignStart + static_cast<int>(nodeAlign.referenceLength());

            const int overlapStart = std::max(alignStart, rightWindowStart);
            const int overlapEnd = std::min(alignEnd, rightWindowEnd);

            if (overlapEnd > overlapStart)
            {
                rightCoveredBases += computeCoveredBasesInRange(
                    nodeAlign, overlapStart - alignStart, overlapEnd - alignStart);
            }
        }
    }
}

// Compute match rate for an entire GraphAlignment (across all nodes)
// matchRate = matchedBases / referenceLength
// Returns 0.0 if referenceLength is 0 (treat as low quality)
static double computeGraphAlignmentMatchRate(const graphtools::GraphAlignment& graphAlign)
{
    size_t matchedBases = 0;
    size_t totalRefLength = 0;

    for (const auto& nodeAlign : graphAlign.alignments())
    {
        matchedBases += nodeAlign.numMatched();
        totalRefLength += nodeAlign.referenceLength();
    }

    if (totalRefLength == 0)
    {
        return 0.0;  // Treat as low quality if no reference alignment
    }
    return static_cast<double>(matchedBases) / static_cast<double>(totalRefLength);
}

// Check if the genotype is homozygous or hemizygous (all paths are identical)
// For such genotypes, all reads are unambiguous by definition since there's only
// one distinct allele they could support
static bool isHomozygousOrHemizygous(const GraphPaths& paths)
{
    if (paths.size() <= 1)
    {
        return true;  // Hemizygous (single allele)
    }

    // Check if all paths are identical (homozygous)
    const auto& firstPath = paths.front();
    for (size_t i = 1; i < paths.size(); ++i)
    {
        if (!(paths[i] == firstPath))
        {
            return false;  // Heterozygous (different alleles)
        }
    }
    return true;  // All paths identical (homozygous)
}

// Check if a fragment is unambiguous (only consistent with one haplotype)
// A fragment is unambiguous when:
// 1. All of its possible alignments have the same pathIndex
// 2. At least one alignment actually exists (read or mate) to confirm support
static bool isUnambiguousFragment(const vector<FragPathAlign>& fragPathAligns)
{
    if (fragPathAligns.empty())
    {
        return false;
    }

    const int firstPathIndex = fragPathAligns.front().readAlign.pathIndex;
    bool hasAnyAlignment = false;

    for (const auto& fragPathAlign : fragPathAligns)
    {
        if (fragPathAlign.readAlign.pathIndex != firstPathIndex)
        {
            return false;  // Fragment aligns to multiple paths
        }
        // Check if this alignment has at least one actual alignment
        if (fragPathAlign.readAlign.align || fragPathAlign.mateAlign.align)
        {
            hasAnyAlignment = true;
        }
    }

    // Only consider unambiguous if at least one alignment exists
    return hasAnyAlignment;
}

using Genotype = vector<int>;

static map<string, Genotype> getGenotypes(const LocusSpecification& locusSpec, const GraphPaths& hapPaths)
{
    map<string, Genotype> genotypes;
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.nodes().size() == 1);
        const auto strNode = variantSpec.nodes().front();

        Genotype genotype;
        for (const auto& hapPath : hapPaths)
        {
            const int strLen = static_cast<int>(std::count(hapPath.nodeIds().begin(), hapPath.nodeIds().end(), strNode));
            genotype.push_back(strLen);
        }
        genotypes[variantSpec.id()] = genotype;
    }

    return genotypes;
}

static int getTotalMatchesToNode(NodeId targetNode, const vector<GraphAlignPtr>& hapAligns)
{
    int numMatches = 0;
    for (const auto& alignPtr : hapAligns)
    {
        for (size_t nodeIndex = 0; nodeIndex != alignPtr->path().numNodes(); ++nodeIndex)
        {
            auto node = alignPtr->getNodeIdByIndex(nodeIndex);
            if (targetNode != node)
            {
                continue;
            }

            numMatches += static_cast<int>(alignPtr->alignments()[nodeIndex].numMatched());
        }
    }

    return numMatches;
}

static map<string, double>
getAlleleDepths(const LocusSpecification& locusSpec, const GraphPath& hapPath, const vector<GraphAlignPtr>& hapAligns)
{
    map<string, double> alleleDepths;

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.nodes().size() == 1);
        const auto strNode = variantSpec.nodes().front();
        int strLen = static_cast<int>(std::count(hapPath.nodeIds().begin(), hapPath.nodeIds().end(), strNode));
        strLen *= static_cast<int>(locusSpec.regionGraph().nodeSeq(strNode).length());
        const int numMatches = getTotalMatchesToNode(strNode, hapAligns);
        const double depth = strLen > 0 ? numMatches / static_cast<double>(strLen) : 0.0;
        alleleDepths[variantSpec.id()] = depth;
    }

    return alleleDepths;
}

// Single-pass function to accumulate all metrics data
static AllAccumulators accumulateAllMetrics(
    const LocusSpecification& locusSpec,
    const GraphPaths& paths,
    const FragById& fragById,
    const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById)
{
    AllAccumulators accum;
    const double matchRateThreshold = 0.9;

    // Get graph and flank info
    const auto& graph = locusSpec.regionGraph();
    const NodeId leftFlankNode = 0;
    const NodeId rightFlankNode = graph.numNodes() - 1;
    const int leftFlankLen = static_cast<int>(graph.nodeSeq(leftFlankNode).length());
    const int rightFlankLen = static_cast<int>(graph.nodeSeq(rightFlankNode).length());

    // Define 20 bp windows adjacent to repeat
    const int leftWindowStart = std::max(0, leftFlankLen - 20);
    const int leftWindowEnd = leftFlankLen;
    const int rightWindowStart = 0;
    const int rightWindowEnd = std::min(20, rightFlankLen);

    // Build repeat node map for quick lookup
    map<string, NodeId> repeatNodeByVariant;
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        assert(variantSpec.nodes().size() == 1);
        repeatNodeByVariant[variantSpec.id()] = variantSpec.nodes().front();
    }

    // Initialize accumulators
    accum.perHaplotype.resize(paths.size());
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        accum.perVariantHaplotype[variantSpec.id()].resize(paths.size());
    }

    // Single pass over all fragments
    for (size_t fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
    {
        const auto& fragId = fragAssignment.fragIds[fragIndex];
        const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
        const auto& fragPathAlign = fragPathAlignsById.at(fragId)[alignIndex];
        const size_t hapIndex = static_cast<size_t>(fragPathAlign.readAlign.pathIndex);

        auto& hapAccum = accum.perHaplotype[hapIndex];

        // Check if fragment is unambiguous (for high-quality unambiguous count)
        // For homozygous/hemizygous (identical paths), all reads are unambiguous by definition
        // since there's only one distinct allele they could support
        const bool isUnambiguous = isHomozygousOrHemizygous(paths) || isUnambiguousFragment(fragPathAlignsById.at(fragId));

        // Process read alignment
        if (fragPathAlign.readAlign.align)
        {
            const auto& readAlign = *fragPathAlign.readAlign.align;

            // Depth: collect alignment (includes all reads for depth calculation)
            hapAccum.alignments.push_back(fragPathAlign.readAlign.align);

            // Flank depth (includes all reads - intentionally about flanks)
            accumulateFlankDepth(
                readAlign, leftFlankNode, rightFlankNode,
                leftWindowStart, leftWindowEnd, rightWindowStart, rightWindowEnd,
                hapAccum.leftFlankCoveredBases, hapAccum.rightFlankCoveredBases);

            // Per-variant metrics (only for reads overlapping the repeat)
            for (const auto& variantSpec : locusSpec.variantSpecs())
            {
                const auto repeatNodeId = repeatNodeByVariant.at(variantSpec.id());
                if (readAlign.overlapsNode(repeatNodeId))
                {
                    auto& varHapAccum = accum.perVariantHaplotype[variantSpec.id()][hapIndex];
                    varHapAccum.numReadsOverlappingRepeat++;
                    varHapAccum.qualitySum += computeReadQuality(readAlign);
                    varHapAccum.totalInsertedBases += countRepeatNodeInsertions(readAlign, repeatNodeId).insertedBases;
                    varHapAccum.totalDeletedBases += countRepeatNodeDeletions(readAlign, repeatNodeId).deletedBases;
                    // High-quality unambiguous count
                    if (isUnambiguous && computeGraphAlignmentMatchRate(readAlign) >= matchRateThreshold)
                    {
                        varHapAccum.highQualUnambiguousCount++;
                    }
                }
            }
        }

        // Process mate alignment
        if (fragPathAlign.mateAlign.align)
        {
            const auto& mateAlign = *fragPathAlign.mateAlign.align;

            // Depth: collect alignment (includes all reads for depth calculation)
            hapAccum.alignments.push_back(fragPathAlign.mateAlign.align);

            // Flank depth (includes all reads - intentionally about flanks)
            accumulateFlankDepth(
                mateAlign, leftFlankNode, rightFlankNode,
                leftWindowStart, leftWindowEnd, rightWindowStart, rightWindowEnd,
                hapAccum.leftFlankCoveredBases, hapAccum.rightFlankCoveredBases);

            // Per-variant metrics (only for reads overlapping the repeat)
            for (const auto& variantSpec : locusSpec.variantSpecs())
            {
                const auto repeatNodeId = repeatNodeByVariant.at(variantSpec.id());
                if (mateAlign.overlapsNode(repeatNodeId))
                {
                    auto& varHapAccum = accum.perVariantHaplotype[variantSpec.id()][hapIndex];
                    varHapAccum.numReadsOverlappingRepeat++;
                    varHapAccum.qualitySum += computeReadQuality(mateAlign);
                    varHapAccum.totalInsertedBases += countRepeatNodeInsertions(mateAlign, repeatNodeId).insertedBases;
                    varHapAccum.totalDeletedBases += countRepeatNodeDeletions(mateAlign, repeatNodeId).deletedBases;
                    // High-quality unambiguous count
                    if (isUnambiguous && computeGraphAlignmentMatchRate(mateAlign) >= matchRateThreshold)
                    {
                        varHapAccum.highQualUnambiguousCount++;
                    }
                }
            }
        }

        // Strand bias: count only the read strand (not mate) to avoid balanced cancellation
        // Only count if the read overlaps a repeat
        if (fragPathAlign.readAlign.align)
        {
            const auto& readAlign = *fragPathAlign.readAlign.align;
            const auto& frag = fragById.at(fragId);
            for (const auto& variantSpec : locusSpec.variantSpecs())
            {
                const auto repeatNodeId = repeatNodeByVariant.at(variantSpec.id());
                if (readAlign.overlapsNode(repeatNodeId))
                {
                    auto& varHapAccum = accum.perVariantHaplotype[variantSpec.id()][hapIndex];
                    if (frag.read.read.isReversed())
                    {
                        varHapAccum.reverseStrandCount++;
                    }
                    else
                    {
                        varHapAccum.forwardStrandCount++;
                    }
                }
            }
        }
    }

    return accum;
}

MetricsByVariant getMetrics(
    const LocusSpecification& locusSpec, const GraphPaths& paths, const FragById& fragById,
    const FragAssignment& fragAssignment, const FragPathAlignsById& fragPathAlignsById)
{
    MetricsByVariant metricsByVariant;
    const auto genotypes = getGenotypes(locusSpec, paths);

    // Single-pass accumulation of all metrics data
    const auto accum = accumulateAllMetrics(locusSpec, paths, fragById, fragAssignment, fragPathAlignsById);

    // Get graph info for flank window sizes
    const auto& graph = locusSpec.regionGraph();
    const NodeId leftFlankNode = 0;
    const NodeId rightFlankNode = graph.numNodes() - 1;
    const int leftFlankLen = static_cast<int>(graph.nodeSeq(leftFlankNode).length());
    const int rightFlankLen = static_cast<int>(graph.nodeSeq(rightFlankNode).length());
    const int leftWindowSize = std::min(20, leftFlankLen);
    const int rightWindowSize = std::min(20, rightFlankLen);

    // Compute per-haplotype depth and derived metrics
    map<string, vector<double>> genotypeDepths;
    map<string, vector<double>> qualitySums;

    // Initialize quality sums per variant
    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        qualitySums[variantSpec.id()].resize(paths.size(), 0.0);
    }

    for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
    {
        // Compute allele depths using collected alignments
        auto alleleDepths = getAlleleDepths(locusSpec, paths[hapIndex], accum.perHaplotype[hapIndex].alignments);
        for (const auto& variantAndDepth : alleleDepths)
        {
            genotypeDepths[variantAndDepth.first].push_back(variantAndDepth.second);
        }

        // Store quality sums (per-variant, only from reads overlapping the repeat)
        for (const auto& variantSpec : locusSpec.variantSpecs())
        {
            qualitySums[variantSpec.id()][hapIndex] = accum.perVariantHaplotype.at(variantSpec.id())[hapIndex].qualitySum;
        }
    }

    // Compute QD from quality sums and depths
    const auto qdValues = computeQD(locusSpec, qualitySums, genotypeDepths);

    for (const auto& variantSpec : locusSpec.variantSpecs())
    {
        const string& variantId = variantSpec.id();
        Metrics metrics;
        metrics.variantId = variantId;
        metrics.genotype = genotypes.at(variantId);
        metrics.alleleDepth = genotypeDepths.at(variantId);
        metrics.qd = qdValues.at(variantId);

        // Compute mean insertion/deletion from accumulated data
        const auto& varHapAccums = accum.perVariantHaplotype.at(variantId);
        metrics.meanInsertedBasesWithinRepeats.resize(paths.size());
        metrics.meanDeletedBasesWithinRepeats.resize(paths.size());
        for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
        {
            const auto& varHapAccum = varHapAccums[hapIndex];
            metrics.meanInsertedBasesWithinRepeats[hapIndex] = (varHapAccum.numReadsOverlappingRepeat > 0)
                ? static_cast<double>(varHapAccum.totalInsertedBases) / varHapAccum.numReadsOverlappingRepeat
                : 0.0;
            metrics.meanDeletedBasesWithinRepeats[hapIndex] = (varHapAccum.numReadsOverlappingRepeat > 0)
                ? static_cast<double>(varHapAccum.totalDeletedBases) / varHapAccum.numReadsOverlappingRepeat
                : 0.0;
        }

        // Strand bias (only from reads overlapping the repeat)
        metrics.strandBiasBinomialPhred.resize(paths.size());
        metrics.forwardStrandReads.resize(paths.size());
        metrics.reverseStrandReads.resize(paths.size());
        for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
        {
            const auto& varHapAccum = varHapAccums[hapIndex];
            metrics.strandBiasBinomialPhred[hapIndex] = computeStrandBiasBinomialPhred(
                varHapAccum.forwardStrandCount, varHapAccum.reverseStrandCount);
            metrics.forwardStrandReads[hapIndex] = varHapAccum.forwardStrandCount;
            metrics.reverseStrandReads[hapIndex] = varHapAccum.reverseStrandCount;
        }

        // Normalize flank depths by allele depth
        const auto& alleleDepths = genotypeDepths.at(variantId);
        metrics.leftFlankNormalizedDepth.resize(paths.size());
        metrics.rightFlankNormalizedDepth.resize(paths.size());

        for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
        {
            const auto& hapAccum = accum.perHaplotype[hapIndex];
            // Compute mean flank depth
            double leftFlankMeanDepth = (leftWindowSize > 0)
                ? static_cast<double>(hapAccum.leftFlankCoveredBases) / leftWindowSize
                : 0.0;
            double rightFlankMeanDepth = (rightWindowSize > 0)
                ? static_cast<double>(hapAccum.rightFlankCoveredBases) / rightWindowSize
                : 0.0;

            if (alleleDepths[hapIndex] > 0.0)
            {
                metrics.leftFlankNormalizedDepth[hapIndex] = leftFlankMeanDepth / alleleDepths[hapIndex];
                metrics.rightFlankNormalizedDepth[hapIndex] = rightFlankMeanDepth / alleleDepths[hapIndex];
            }
            else
            {
                // If allele depth is zero, normalized depth is undefined; use 0.0
                metrics.leftFlankNormalizedDepth[hapIndex] = 0.0;
                metrics.rightFlankNormalizedDepth[hapIndex] = 0.0;
            }
        }

        // High-quality unambiguous reads (only those overlapping the repeat)
        metrics.highQualityUnambiguousReads.resize(paths.size());
        for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
        {
            metrics.highQualityUnambiguousReads[hapIndex] = varHapAccums[hapIndex].highQualUnambiguousCount;
        }

        // Build CountTable of high-quality unambiguous reads by allele size
        // Uses genotype to map haplotype index to allele size (in repeat units)
        for (size_t hapIndex = 0; hapIndex < paths.size(); ++hapIndex)
        {
            int alleleSize = metrics.genotype[hapIndex];
            int count = varHapAccums[hapIndex].highQualUnambiguousCount;
            if (count > 0)
            {
                metrics.countsOfHighQualityUnambiguousReads.incrementCountOf(alleleSize, count);
            }
        }

        metricsByVariant.push_back(metrics);
    }

    return metricsByVariant;
}

}  // namespace reviewer
}  // namespace ehunter
