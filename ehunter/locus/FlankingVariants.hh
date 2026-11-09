//
// ExpansionHunter
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

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include "core/GenomicRegion.hh"

namespace ehunter
{

// Small variants found in the two flanks of a repeat locus by --output-flanking-variants, with the reads of each
// repeat allele that show each variant's REF and ALT allele. The variants are phased with the repeat alleles by
// the reads themselves: a read only counts toward an allele when its read pair's alignment score in the REViewer
// step is best on that allele's path alone (see reviewer/FlankingVariantsByAllele.hh).

// A SNV or indel in VCF form: the 0-based position of the first REF base (for an indel, the base before the
// inserted or deleted bases), with REF and ALT in upper case. Indels are shifted as far left as the flank's
// reference sequence allows.
struct FlankingVariant
{
    int64_t position = 0;
    std::string refAllele;
    std::string altAllele;

    bool operator<(const FlankingVariant& other) const
    {
        return std::tie(position, refAllele, altAllele)
            < std::tie(other.position, other.refAllele, other.altAllele);
    }
    bool operator==(const FlankingVariant& other) const
    {
        return position == other.position && refAllele == other.refAllele && altAllele == other.altAllele;
    }
};

struct FlankingVariantReadCounts
{
    FlankingVariant variant;
    int refReads = 0;  // reads of this repeat allele that show the REF allele
    int altReads = 0;  // reads of this repeat allele that show the ALT allele
};

// The FlankingVariants entry of one repeat allele
struct AlleleFlankingVariants
{
    int alleleSize = 0;  // repeat units, as in the genotype
    // From the outermost left-flank position to the outermost right-flank position covered by at least
    // reviewer::kMinReadsToCoverFlankPosition of this allele's reads, including the repeat in between. When only one flank
    // has such a position the region stops at the repeat edge on the other side. Unset when neither flank has one.
    std::optional<GenomicRegion> haplotypeRegion;
    // Every variant listed for the locus, in the same order under each allele
    std::vector<FlankingVariantReadCounts> variants;
};

}  // namespace ehunter
