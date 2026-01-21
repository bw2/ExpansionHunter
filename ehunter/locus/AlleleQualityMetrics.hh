//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
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

#pragma once

#include <array>
#include <string>
#include <vector>

namespace ehunter
{

// Per-allele quality metrics computed from aligned reads
struct AlleleMetrics
{
    int alleleNumber = 0;  // 1-based: 1 for short/only allele, 2 for long allele in het calls
    int alleleSize = 0;

    // DP: allele-specific depth (count of reads supporting this allele)
    double depth = 0.0;

    // QD: read quality / depth (average read quality for allele-supporting reads)
    double qd = 0.0;

    // Mean inserted bases in the repeat region per read for allele-supporting reads
    double meanInsertedBasesWithinRepeats = 0.0;

    // Mean deleted bases in the repeat node per read
    double meanDeletedBasesWithinRepeats = 0.0;

    // Binomial test p-value (phred-scaled) for strand bias (deviation from 50/50)
    double strandBiasBinomialPhred = 0.0;

    // Left flank depth (20 bp window) normalized by allele depth
    double leftFlankNormalizedDepth = 0.0;

    // Right flank depth (20 bp window) normalized by allele depth
    double rightFlankNormalizedDepth = 0.0;

    // Count of high-quality reads that unambiguously support this allele
    int highQualityUnambiguousReads = 0;

    // Confidence interval width divided by allele size
    // Formula: (CI_upper - CI_lower) / AlleleSize
    double confidenceIntervalDividedByAlleleSize = 0.0;
};

// Container for all allele quality metrics for a repeat variant
// Homozygous/hemizygous calls produce one AlleleMetrics; heterozygous calls produce two
struct RepeatAlleleQualityMetrics
{
    std::string variantId;
    std::vector<AlleleMetrics> alleles; // size 1 for hom/hemi, size 2 for het
    bool hasMetrics = false;
};

}
