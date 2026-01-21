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

#include <string>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "locus/LocusSpecification.hh"
#include "reviewer/Aligns.hh"

namespace ehunter
{
namespace reviewer
{

/// Compute read quality score for a GraphAlignment
/// Returns a value in [0.0, 1.0] where 1.0 means perfect alignment (all matches)
/// Formula: readQuality = matchedBases / totalBases
/// matchedBases = bases with OperationType::kMatch
/// totalBases = GraphAlignment::queryLength()
/// Returns 0.0 if totalBases is zero
double computeReadQuality(const graphtools::GraphAlignment& graphAlign);

/// Compute strand bias binomial test p-value (phred-scaled)
/// Tests whether forward/reverse counts deviate significantly from 50/50
/// Returns -10 * log10(p_value). Returns 0.0 if no reads.
double computeStrandBiasBinomialPhred(int forwardReads, int reverseReads);

/// Metrics for a single variant
struct Metrics
{
    std::string variantId = "NA";
    std::vector<int> genotype;
    std::vector<double> alleleDepth;
    // QD: sum of read quality divided by allele depth, per haplotype
    // If alleleDepth is zero, QD is 0.0
    std::vector<double> qd;
    // Mean inserted bases in the repeat region per read, per haplotype
    // If readCount is zero, meanInsertedBasesWithinRepeats is 0.0
    std::vector<double> meanInsertedBasesWithinRepeats;
    // Mean deleted bases in the repeat node per read, per haplotype
    // If no reads overlap the repeat node, meanDeletedBasesWithinRepeats is 0.0
    std::vector<double> meanDeletedBasesWithinRepeats;
    // Strand bias: Binomial test p-value (phred-scaled) testing whether the forward/reverse
    // strand distribution deviates significantly from 50/50. Value = -10 * log10(p).
    // Higher values indicate greater strand bias. Returns 0.0 if no reads.
    std::vector<double> strandBiasBinomialPhred;
    // Raw forward strand read counts per haplotype (for proper statistical combination in homozygous case)
    std::vector<int> forwardStrandReads;
    // Raw reverse strand read counts per haplotype (for proper statistical combination in homozygous case)
    std::vector<int> reverseStrandReads;
    // Left flank depth (20 bp window) normalized by allele depth, per haplotype
    std::vector<double> leftFlankNormalizedDepth;
    // Right flank depth (20 bp window) normalized by allele depth, per haplotype
    std::vector<double> rightFlankNormalizedDepth;
    // Count of high-quality unambiguous reads per haplotype
    // High-quality: matchRate >= 0.9 (same threshold as overhang filtering)
    // Unambiguous: fragment is consistent with only one haplotype (singlePath == true)
    std::vector<int> highQualityUnambiguousReads;
};

using MetricsByVariant = std::vector<Metrics>;

/// Calculate metrics for all variants in a locus
/// @param locusSpec Locus specification
/// @param paths Haplotype paths
/// @param fragById Map of fragments by ID
/// @param fragAssignment Fragment origin assignment
/// @param fragPathAlignsById Fragment path alignments
/// @return Metrics for each variant
MetricsByVariant getMetrics(
    const LocusSpecification& locusSpec,
    const GraphPaths& paths,
    const FragById& fragById,
    const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById);

}  // namespace reviewer
}  // namespace ehunter
