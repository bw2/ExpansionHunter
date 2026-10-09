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

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/LocusStats.hh"
#include "locus/VariantFindings.hh"

namespace ehunter
{

// For a locus with more than one repeat variant (e.g. HTT's CAG and CCG repeats): which allele of each repeat
// sits on the same haplotype as which allele of the others. The pairing is the one the REViewer step picks
// (reviewer/Phasing.cpp), since each variant is otherwise genotyped independently with its alleles sorted by size.
struct RepeatAllelePhasing
{
    // One entry per haplotype, mapping variant id to that haplotype's repeat size in repeat units. Haplotypes are
    // ordered by their repeat sizes, compared in catalog variant order.
    std::vector<std::map<std::string, int>> repeatSizeByVariantIdOnEachHaplotype;
    // Distinct ways of pairing the genotyped alleles onto haplotypes: 1 when at most one repeat is heterozygous,
    // and the fields below are only meaningful when it is more than 1
    int numberOfPossiblePairings = 0;
    // Fragments whose best alignment scores higher on the chosen pairing's haplotypes than on the next-best
    // pairing's haplotypes, and vice versa
    int fragmentsSupportingChosenPairingOverNextBest = 0;
    int fragmentsSupportingNextBestPairingOverChosen = 0;
    // True when the chosen and next-best pairings got equal REViewer scores, so the choice between them is arbitrary
    bool chosenPairingIsTiedWithNextBest = false;
};

// Container with per-locus analysis results
struct LocusFindings
{
    explicit LocusFindings(LocusStats stats = {})
        : stats(stats)
    {
    }
    LocusStats stats;
    // VariantFindings is an abstract class from which findings for all variant types are derived
    std::unordered_map<std::string, std::unique_ptr<VariantFindings>> findingsForEachVariant;
    // Set only for loci with more than one repeat variant, when the REViewer step ran successfully
    std::optional<RepeatAllelePhasing> repeatAllelePhasing;
};

using SampleFindings = std::vector<LocusFindings>;

}
