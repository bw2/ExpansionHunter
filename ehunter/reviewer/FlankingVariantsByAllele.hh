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

// --output-flanking-variants: SNVs and indels in the two flanks of a locus with a single repeat, with the reads of
// each repeat allele that show each variant's REF and ALT allele.
//
// The reads come from the REViewer step, which aligns each read pair to the haplotype paths of the called genotype.
// A read pair counts toward an allele when its alignment score is best on that allele's path alone, and then both of
// its reads count, so the variants are phased with the repeat alleles by the reads themselves. Read pairs that score
// equally well on both paths are not used, even though REViewer then picks a path for them by fragment length: when
// the two alleles differ by less than the spread of fragment lengths, that pick is close to random. These include
// pairs that lie entirely in a flank and pairs with one read in each flank. When the genotype is homozygous or
// hemizygous there is only one allele and every read pair counts.

#pragma once

#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "graphcore/Graph.hh"

#include "genotyping/RepeatGenotype.hh"
#include "locus/FlankingVariants.hh"
#include "locus/LocusSpecification.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/GenotypePaths.hh"

namespace ehunter
{
namespace reviewer
{

// A variant is listed when at least this many reads of one repeat allele show its ALT allele. One read is not
// enough because sequencing errors would then be listed as variants.
constexpr int kMinAltReadsToListFlankingVariant = 2;

// HaplotypeRegion extends to the outermost flank positions covered by at least this many reads of the allele
constexpr int kMinReadsToCoverFlankPosition = 2;

// A read's differences from the reference that are less than this many bases apart (2 SNVs at offsets 10 and 14 are
// 4 bases apart) are merged into one variant: the reference stretch they cover with the read's sequence in its
// place. A change that the aligner splits into several SNVs and indels, differently in different reads, is then
// reported as one variant. Two real SNVs this close together are reported as one variant too.
constexpr int kMergeDifferencesCloserThan = 5;

/// Finds SNVs and indels in the two flanks of a locus with a single repeat and counts, for each repeat allele, the
/// reads that show each variant's REF and ALT allele.
///
/// The left flank is the graph's first node and the right flank its last node. Each read's differences from the
/// flank's reference sequence are merged into variants as described for kMergeDifferencesCloserThan. Mismatched
/// bases that are low quality (lower case) or N, and insertions that contain one, are not differences, so they
/// neither give a variant nor merge with one. A read shows the ALT allele of a variant when it
/// carries that variant. It shows the REF allele of an SNV or a longer substitution when it has the high-quality
/// (upper case) reference bases there, and the REF allele of an indel when it aligns through every position the
/// indel could be placed at without an insertion or deletion there.
///
/// @param graph Locus graph
/// @param referenceRegionOfEachNode Reference region of each graph node; must include both flank nodes
/// @param referenceBaseAfterLeftFlank Reference base just after the left flank ('N' when unknown)
/// @param referenceBaseBeforeRightFlank Reference base just before the right flank, used as the VCF anchor base of
///        an indel at the very start of the right flank ('N' when unknown, in which case such indels are skipped).
///        These two bases also show which indels at the edge of a flank could equally well be placed in the
///        repeat; such indels change the repeat sequence as much as the flank, and are not reported.
/// @param repeatNode The repeat's graph node
/// @param genotype The repeat's genotype
/// @param diplotype The haplotype paths chosen by the REViewer step
/// @param fragById Read pairs by fragment id
/// @param fragPathAlignsById Each read pair's alignments to the haplotype paths
/// @param fragAssignment The alignment the REViewer step picked for each read pair
/// @param fragIdsPlacedOnOnePathByAlignmentScore The read pairs whose alignment score is best on one path alone
/// @return One entry per distinct allele, shorter allele first. Empty when the genotype is heterozygous but the reads
///         can't be split by allele because both allele paths hold the same number of repeat units (both alleles
///         were capped to the same path length).
std::optional<std::vector<AlleleFlankingVariants>> findFlankingVariantsByAllele(
    const graphtools::Graph& graph, const NodeToRegionAssociation& referenceRegionOfEachNode,
    char referenceBaseAfterLeftFlank, char referenceBaseBeforeRightFlank, graphtools::NodeId repeatNode,
    const RepeatGenotype& genotype, const Diplotype& diplotype, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById, const FragAssignment& fragAssignment,
    const std::set<std::string>& fragIdsPlacedOnOnePathByAlignmentScore);

/// Shifts a deletion of reference[start, end) left for as long as the deleted bases rotated by one base still match
/// the reference, which gives the leftmost of the equivalent positions.
/// @return The new start; the deletion keeps its length
int shiftDeletionLeft(const std::string& reference, int start, int end);

/// Shifts an insertion of insertedBases placed just before reference[start] left in the same way.
/// @return The new start and the inserted bases as they read at that position
std::pair<int, std::string> shiftInsertionLeft(const std::string& reference, int start, std::string insertedBases);

}  // namespace reviewer
}  // namespace ehunter
