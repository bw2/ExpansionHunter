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

#include "reviewer/FlankingVariantsByAllele.hh"

#include <algorithm>
#include <cctype>
#include <map>
#include <set>

#include "graphalign/GraphAlignment.hh"
#include "graphalign/Operation.hh"

using std::string;
using std::vector;

namespace ehunter
{
namespace reviewer
{

namespace
{

// Recorded at the repeat position next to a flank when a read aligns across the boundary between that flank and the
// repeat, so that an insertion at the boundary can be checked against reads that cross it. It is never a REF base.
constexpr char kBaseAlignedOutsideFlank = '*';

bool isUpperCaseNucleotide(char base) { return base == 'A' || base == 'C' || base == 'G' || base == 'T'; }

bool areAllUpperCaseNucleotides(const string& bases)
{
    return !bases.empty() && std::all_of(bases.begin(), bases.end(), isUpperCaseNucleotide);
}

struct Flank
{
    graphtools::NodeId node;
    string sequence;  // reference sequence of the flank, upper case
    GenomicRegion region;
    bool isRightFlank;
    // The reference base next to the flank on the repeat's side: the base just after the left flank, or just before
    // the right flank. 'N' when unknown.
    char referenceBaseOnRepeatSide;
};

// The positions a read must align, with nothing inserted between them, to show a variant's REF allele. For an SNV,
// or a substitution of several bases, these are the positions it covers, where the read must also have the
// reference bases (in upper case, so that a low-quality base never counts). For an indel they are every position
// the indel can be shifted across without changing the sequence (the same indel can be placed anywhere in a
// homopolymer or short tandem repeat), plus, for an insertion, the base on each side of that stretch.
struct PositionsShowingRefAllele
{
    int64_t first;
    int64_t last;
    string referenceBases;  // the reference bases at [first, last] that the read must have; empty for an indel
};

// One difference between a read and a flank's reference sequence, at flank offsets [start, end): a mismatched base
// (end = start + 1, altBases = the read's base), a deletion (altBases empty) or an insertion just before offset
// `start` (end = start). Mismatches and insertions with a low-quality or N base are not recorded as differences.
struct FlankDifference
{
    int start;
    int end;
    string altBases;
};

// What one read shows in the two flanks, in genomic coordinates
struct ReadFlankObservations
{
    // The read's base at each flank position it aligns with a match or mismatch (lower case when low quality), plus
    // kBaseAlignedOutsideFlank at the repeat position next to a flank when the read aligns across that boundary
    std::map<int64_t, char> alignedBaseByPosition;
    // Positions that an insertion in the read sits just before
    std::set<int64_t> positionsWithInsertionBefore;
    std::set<FlankingVariant> variants;
    // First and last offsets the read covers on each flank node (matches, mismatches and deletions)
    std::optional<std::pair<int, int>> extentOnLeftFlank;
    std::optional<std::pair<int, int>> extentOnRightFlank;
};

int shiftDeletionRight(const string& reference, int start, int end)
{
    while (end < static_cast<int>(reference.size()) && reference[start] == reference[end])
    {
        ++start;
        ++end;
    }
    return start;
}

// Shifts an insertion of insertedBases placed just before reference[start] right, the mirror image of
// shiftInsertionLeft. Returns the new start and the inserted bases as they read at that position.
std::pair<int, string> shiftInsertionRight(const string& reference, int start, string insertedBases)
{
    while (start < static_cast<int>(reference.size()) && reference[start] == insertedBases.front())
    {
        insertedBases = insertedBases.substr(1) + insertedBases.front();
        ++start;
    }
    return { start, insertedBases };
}

// The VCF anchor base for an indel whose leftmost position is just before flank offset `start`, with the anchor's
// genomic position. Empty when there is no anchor base: an indel at the very start of the left flank, at the very
// start of the right flank when the base before it is unknown, or next to an N in the reference.
std::optional<std::pair<int64_t, char>> getAnchorBase(const Flank& flank, int start)
{
    if (start > 0)
    {
        const char base = flank.sequence[start - 1];
        if (!isUpperCaseNucleotide(base))
        {
            return std::nullopt;
        }
        return std::make_pair(flank.region.start() + start - 1, base);
    }
    if (flank.isRightFlank && isUpperCaseNucleotide(flank.referenceBaseOnRepeatSide))
    {
        return std::make_pair(flank.region.start() - 1, flank.referenceBaseOnRepeatSide);
    }
    return std::nullopt;
}

// Records a deletion of flank.sequence[start, end), together with the positions that show its REF allele. A deletion
// that could equally well be placed in the repeat, such as the deletion of a G at the start of the right flank next to
// a CAG repeat, is not recorded: it changes the repeat sequence as much as the flank, and reads that place it in the
// repeat would show neither allele correctly.
void addDeletion(
    const Flank& flank, int start, int end, ReadFlankObservations& observations,
    std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    const int flankLength = static_cast<int>(flank.sequence.size());
    if (end > flankLength)
    {
        return;
    }
    const int length = end - start;
    const int leftmostStart = shiftDeletionLeft(flank.sequence, start, end);
    const int rightmostStart = shiftDeletionRight(flank.sequence, leftmostStart, leftmostStart + length);
    const bool canShiftIntoRepeat = flank.isRightFlank
        ? leftmostStart == 0 && flank.sequence[length - 1] == flank.referenceBaseOnRepeatSide
        : rightmostStart + length == flankLength && flank.sequence[rightmostStart] == flank.referenceBaseOnRepeatSide;
    const auto anchor = getAnchorBase(flank, leftmostStart);
    const string deletedBases = flank.sequence.substr(leftmostStart, length);
    if (canShiftIntoRepeat || !anchor || !areAllUpperCaseNucleotides(deletedBases))
    {
        return;
    }

    const FlankingVariant variant { anchor->first, anchor->second + deletedBases, string(1, anchor->second) };
    observations.variants.insert(variant);
    positionsShowingRefAlleleByVariant.emplace(
        variant,
        PositionsShowingRefAllele { flank.region.start() + leftmostStart,
                                    flank.region.start() + rightmostStart + length - 1, string() });
}

// Records an insertion of insertedBases just before flank.sequence[start], together with the positions that show its
// REF allele. As with deletions, an insertion that could equally well be placed in the repeat is not recorded.
void addInsertion(
    const Flank& flank, int start, const string& insertedBases, ReadFlankObservations& observations,
    std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    if (!areAllUpperCaseNucleotides(insertedBases))
    {
        return;
    }
    const auto leftmostStartAndBases = shiftInsertionLeft(flank.sequence, start, insertedBases);
    const int leftmostStart = leftmostStartAndBases.first;
    const auto rightmostStartAndBases
        = shiftInsertionRight(flank.sequence, leftmostStart, leftmostStartAndBases.second);
    const int rightmostStart = rightmostStartAndBases.first;
    const bool canShiftIntoRepeat = flank.isRightFlank
        ? leftmostStart == 0 && leftmostStartAndBases.second.back() == flank.referenceBaseOnRepeatSide
        : rightmostStart == static_cast<int>(flank.sequence.size())
            && rightmostStartAndBases.second.front() == flank.referenceBaseOnRepeatSide;
    const auto anchor = getAnchorBase(flank, leftmostStart);
    if (canShiftIntoRepeat || !anchor)
    {
        return;
    }

    const FlankingVariant variant { anchor->first, string(1, anchor->second),
                                    anchor->second + leftmostStartAndBases.second };
    observations.variants.insert(variant);
    positionsShowingRefAlleleByVariant.emplace(
        variant,
        PositionsShowingRefAllele { flank.region.start() + leftmostStart - 1, flank.region.start() + rightmostStart,
                                    string() });
}

// Records a substitution of flank.sequence[start, start + refAllele.size()) by altAllele: an SNV when both are one
// base long
void addSubstitution(
    const Flank& flank, int start, const string& refAllele, const string& altAllele,
    ReadFlankObservations& observations,
    std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    if (!areAllUpperCaseNucleotides(refAllele) || !areAllUpperCaseNucleotides(altAllele))
    {
        return;
    }
    const int64_t position = flank.region.start() + start;
    const FlankingVariant variant { position, refAllele, altAllele };
    observations.variants.insert(variant);
    positionsShowingRefAlleleByVariant.emplace(
        variant,
        PositionsShowingRefAllele { position, position + static_cast<int64_t>(refAllele.size()) - 1, refAllele });
}

// Records the variant that a group of nearby differences between a read and a flank's reference sequence add up to:
// the reference stretch they cover with the read's sequence in its place, trimmed of the bases the two share at
// either end. Reads that split the same change into different SNVs and indels then still give the same variant. A
// single difference gives an SNV or an indel.
void addVariantFromDifferences(
    const Flank& flank, const vector<FlankDifference>& differences, ReadFlankObservations& observations,
    std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    // The differences are in reference order and don't overlap, so the read's sequence over [start, end) is the
    // reference between them with each one's bases in its place
    const int start = differences.front().start;
    const int end = differences.back().end;
    if (end > static_cast<int>(flank.sequence.size()))
    {
        return;
    }
    string altAllele;
    int cursor = start;
    for (const auto& difference : differences)
    {
        altAllele.append(flank.sequence, cursor, difference.start - cursor);
        altAllele += difference.altBases;
        cursor = difference.end;
    }
    string refAllele = flank.sequence.substr(start, end - start);

    // Trim the bases REF and ALT share at the end, then at the start
    while (!refAllele.empty() && !altAllele.empty() && refAllele.back() == altAllele.back())
    {
        refAllele.pop_back();
        altAllele.pop_back();
    }
    size_t sharedStart = 0;
    while (sharedStart < refAllele.size() && sharedStart < altAllele.size()
           && refAllele[sharedStart] == altAllele[sharedStart])
    {
        ++sharedStart;
    }
    refAllele.erase(0, sharedStart);
    altAllele.erase(0, sharedStart);
    const int trimmedStart = start + static_cast<int>(sharedStart);

    if (refAllele.empty() && altAllele.empty())
    {
        return;
    }
    if (refAllele.empty())
    {
        addInsertion(flank, trimmedStart, altAllele, observations, positionsShowingRefAlleleByVariant);
    }
    else if (altAllele.empty())
    {
        addDeletion(
            flank, trimmedStart, trimmedStart + static_cast<int>(refAllele.size()), observations,
            positionsShowingRefAlleleByVariant);
    }
    else
    {
        addSubstitution(flank, trimmedStart, refAllele, altAllele, observations, positionsShowingRefAlleleByVariant);
    }
}

// Splits a read's differences on one flank into groups of differences less than kMergeDifferencesCloserThan bases
// apart, and records the variant each group adds up to
void addVariantsFromDifferences(
    const Flank& flank, const vector<FlankDifference>& differences, ReadFlankObservations& observations,
    std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    size_t groupStart = 0;
    while (groupStart < differences.size())
    {
        size_t groupEnd = groupStart + 1;
        int groupEndOffset = differences[groupStart].end;
        // The distance from the group's last reference base to the next difference, so that 2 SNVs at offsets 10 and
        // 14 are 4 bases apart
        while (groupEnd < differences.size()
               && differences[groupEnd].start - groupEndOffset + 1 < kMergeDifferencesCloserThan)
        {
            groupEndOffset = std::max(groupEndOffset, differences[groupEnd].end);
            ++groupEnd;
        }
        const vector<FlankDifference> group(differences.begin() + groupStart, differences.begin() + groupEnd);
        addVariantFromDifferences(flank, group, observations, positionsShowingRefAlleleByVariant);
        groupStart = groupEnd;
    }
}

void extendExtent(std::optional<std::pair<int, int>>& extent, int offset)
{
    if (!extent)
    {
        extent = std::make_pair(offset, offset);
        return;
    }
    extent->first = std::min(extent->first, offset);
    extent->second = std::max(extent->second, offset);
}

bool isMatchOrMismatch(const graphtools::Operation& operation)
{
    return operation.type() == graphtools::OperationType::kMatch
        || operation.type() == graphtools::OperationType::kMismatch;
}

// Records how a read crosses the boundary between a flank and the repeat. `operationNextToFlank` is the operation of
// the repeat-side node that touches the flank, `positionOutsideFlank` the repeat position next to the flank, and
// `positionAfterBoundary` the position just right of the boundary (where an insertion at the boundary is recorded).
void recordBoundaryWithRepeat(
    const graphtools::Operation& operationNextToFlank, int64_t positionOutsideFlank, int64_t positionAfterBoundary,
    ReadFlankObservations& observations)
{
    if (isMatchOrMismatch(operationNextToFlank))
    {
        observations.alignedBaseByPosition.emplace(positionOutsideFlank, kBaseAlignedOutsideFlank);
    }
    else if (operationNextToFlank.type() == graphtools::OperationType::kInsertionToRef)
    {
        // An insertion placed on the repeat's side of the boundary. Recorded so that this read doesn't count as
        // showing the REF allele of an insertion at the boundary. The insertion itself is not reported, since it
        // is in the repeat.
        observations.positionsWithInsertionBefore.insert(positionAfterBoundary);
    }
}

ReadFlankObservations observeReadInFlanks(
    const string& readBases, const graphtools::GraphAlignment& alignment, const Flank& leftFlank,
    const Flank& rightFlank, std::map<FlankingVariant, PositionsShowingRefAllele>& positionsShowingRefAlleleByVariant)
{
    ReadFlankObservations observations;
    const auto& nodeIds = alignment.path().nodeIds();
    const auto& nodeAlignments = alignment.alignments();

    int queryPosition = 0;
    for (size_t nodeIndex = 0; nodeIndex != nodeIds.size(); ++nodeIndex)
    {
        const graphtools::Alignment& nodeAlignment = nodeAlignments[nodeIndex];
        const bool isLeftFlank = nodeIds[nodeIndex] == leftFlank.node;
        const bool isRightFlank = nodeIds[nodeIndex] == rightFlank.node;
        if (!isLeftFlank && !isRightFlank)
        {
            queryPosition += static_cast<int>(nodeAlignment.queryLength());
            continue;
        }

        const Flank& flank = isLeftFlank ? leftFlank : rightFlank;
        auto& extent = isLeftFlank ? observations.extentOnLeftFlank : observations.extentOnRightFlank;
        const int flankLength = static_cast<int>(flank.sequence.size());

        if (isRightFlank && nodeIndex > 0 && nodeAlignments[nodeIndex - 1].numOperations() > 0)
        {
            recordBoundaryWithRepeat(
                nodeAlignments[nodeIndex - 1].operations().back(), flank.region.start() - 1, flank.region.start(),
                observations);
        }

        vector<FlankDifference> differences;
        int offset = static_cast<int>(nodeAlignment.referenceStart());
        for (const graphtools::Operation& operation : nodeAlignment.operations())
        {
            const int length = static_cast<int>(operation.length());
            switch (operation.type())
            {
            case graphtools::OperationType::kMatch:
            case graphtools::OperationType::kMismatch:
                for (int index = 0; index != length; ++index, ++offset, ++queryPosition)
                {
                    if (offset >= flankLength || queryPosition >= static_cast<int>(readBases.size()))
                    {
                        continue;
                    }
                    const char readBase = readBases[queryPosition];
                    const char referenceBase = flank.sequence[offset];
                    observations.alignedBaseByPosition[flank.region.start() + offset] = readBase;
                    if (isUpperCaseNucleotide(readBase) && isUpperCaseNucleotide(referenceBase)
                        && readBase != referenceBase)
                    {
                        differences.push_back({ offset, offset + 1, string(1, readBase) });
                    }
                    extendExtent(extent, offset);
                }
                break;
            case graphtools::OperationType::kInsertionToRef:
            {
                const string insertedBases = queryPosition < static_cast<int>(readBases.size())
                    ? readBases.substr(queryPosition, length) : string();
                observations.positionsWithInsertionBefore.insert(flank.region.start() + offset);
                if (areAllUpperCaseNucleotides(insertedBases))
                {
                    differences.push_back({ offset, offset, insertedBases });
                }
                queryPosition += length;
                break;
            }
            case graphtools::OperationType::kDeletionFromRef:
                if (offset + length <= flankLength)
                {
                    differences.push_back({ offset, offset + length, string() });
                }
                extendExtent(extent, offset);
                extendExtent(extent, std::min(offset + length, flankLength) - 1);
                offset += length;
                break;
            case graphtools::OperationType::kSoftclip:
                queryPosition += length;
                break;
            case graphtools::OperationType::kMissingBases:
                queryPosition += length;
                offset += length;
                break;
            }
        }
        addVariantsFromDifferences(flank, differences, observations, positionsShowingRefAlleleByVariant);

        if (isLeftFlank && offset == flankLength && nodeIndex + 1 < nodeIds.size()
            && nodeAlignments[nodeIndex + 1].numOperations() > 0)
        {
            recordBoundaryWithRepeat(
                nodeAlignments[nodeIndex + 1].operations().front(), flank.region.end(), flank.region.end(),
                observations);
        }
    }

    return observations;
}

bool showsRefAllele(const ReadFlankObservations& observations, const PositionsShowingRefAllele& positions)
{
    for (int64_t position = positions.first; position <= positions.last; ++position)
    {
        const auto alignedBase = observations.alignedBaseByPosition.find(position);
        if (alignedBase == observations.alignedBaseByPosition.end())
        {
            return false;
        }
        if (position > positions.first && observations.positionsWithInsertionBefore.count(position))
        {
            return false;
        }
        if (!positions.referenceBases.empty()
            && alignedBase->second != positions.referenceBases[position - positions.first])
        {
            return false;
        }
    }
    return true;
}

// The region from the outermost left-flank position to the outermost right-flank position covered by at least
// kMinReadsToCoverFlankPosition reads
std::optional<GenomicRegion> computeHaplotypeRegion(
    const vector<ReadFlankObservations>& reads, const Flank& leftFlank, const Flank& rightFlank)
{
    // Read depth at each offset of a flank, from the first and last offset each read covers
    auto computeDepth = [&reads](const Flank& flank) {
        const int flankLength = static_cast<int>(flank.sequence.size());
        vector<int> depthChange(flankLength + 1, 0);
        for (const auto& read : reads)
        {
            const auto& extent = flank.isRightFlank ? read.extentOnRightFlank : read.extentOnLeftFlank;
            if (extent)
            {
                ++depthChange[extent->first];
                --depthChange[extent->second + 1];
            }
        }
        vector<int> depth(flankLength, 0);
        int runningDepth = 0;
        for (int offset = 0; offset != flankLength; ++offset)
        {
            runningDepth += depthChange[offset];
            depth[offset] = runningDepth;
        }
        return depth;
    };

    const vector<int> leftFlankDepth = computeDepth(leftFlank);
    const vector<int> rightFlankDepth = computeDepth(rightFlank);

    std::optional<int> leftmostCoveredOffset;
    for (int offset = 0; offset != static_cast<int>(leftFlankDepth.size()); ++offset)
    {
        if (leftFlankDepth[offset] >= kMinReadsToCoverFlankPosition)
        {
            leftmostCoveredOffset = offset;
            break;
        }
    }
    std::optional<int> rightmostCoveredOffset;
    for (int offset = static_cast<int>(rightFlankDepth.size()) - 1; offset >= 0; --offset)
    {
        if (rightFlankDepth[offset] >= kMinReadsToCoverFlankPosition)
        {
            rightmostCoveredOffset = offset;
            break;
        }
    }

    if (!leftmostCoveredOffset && !rightmostCoveredOffset)
    {
        return std::nullopt;
    }
    const int64_t start
        = leftmostCoveredOffset ? leftFlank.region.start() + *leftmostCoveredOffset : leftFlank.region.end();
    const int64_t end
        = rightmostCoveredOffset ? rightFlank.region.start() + *rightmostCoveredOffset + 1 : rightFlank.region.start();
    return GenomicRegion(leftFlank.region.contigIndex(), start, end);
}

string toUpperCase(string sequence)
{
    std::transform(sequence.begin(), sequence.end(), sequence.begin(), [](unsigned char base) {
        return static_cast<char>(std::toupper(base));
    });
    return sequence;
}

}  // namespace

int shiftDeletionLeft(const string& reference, int start, int end)
{
    while (start > 0 && end <= static_cast<int>(reference.size()) && reference[start - 1] == reference[end - 1])
    {
        --start;
        --end;
    }
    return start;
}

std::pair<int, string> shiftInsertionLeft(const string& reference, int start, string insertedBases)
{
    while (start > 0 && !insertedBases.empty() && start <= static_cast<int>(reference.size())
           && reference[start - 1] == insertedBases.back())
    {
        insertedBases = insertedBases.back() + insertedBases.substr(0, insertedBases.size() - 1);
        --start;
    }
    return { start, insertedBases };
}

std::optional<vector<AlleleFlankingVariants>> findFlankingVariantsByAllele(
    const graphtools::Graph& graph, const NodeToRegionAssociation& referenceRegionOfEachNode,
    char referenceBaseAfterLeftFlank, char referenceBaseBeforeRightFlank, graphtools::NodeId repeatNode,
    const RepeatGenotype& genotype, const Diplotype& diplotype, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById, const FragAssignment& fragAssignment,
    const std::set<string>& fragIdsPlacedOnOnePathByAlignmentScore)
{
    const graphtools::NodeId leftFlankNode = 0;
    const graphtools::NodeId rightFlankNode = graph.numNodes() - 1;
    const auto leftFlankRegion = referenceRegionOfEachNode.find(leftFlankNode);
    const auto rightFlankRegion = referenceRegionOfEachNode.find(rightFlankNode);
    if (graph.numNodes() < 3 || leftFlankRegion == referenceRegionOfEachNode.end()
        || rightFlankRegion == referenceRegionOfEachNode.end())
    {
        return std::nullopt;
    }
    const Flank leftFlank { leftFlankNode, toUpperCase(graph.nodeSeq(leftFlankNode)), leftFlankRegion->second, false,
                            referenceBaseAfterLeftFlank };
    const Flank rightFlank { rightFlankNode, toUpperCase(graph.nodeSeq(rightFlankNode)), rightFlankRegion->second,
                             true, referenceBaseBeforeRightFlank };

    // Which entry (0 = shorter allele) each haplotype path's reads count toward. A homozygous or hemizygous genotype
    // has a single entry that every read pair counts toward.
    const bool hasOneDistinctAllele = genotype.numAlleles() == 1 || genotype.isHomozygous();
    vector<int> entryIndexOfPath;
    if (!hasOneDistinctAllele)
    {
        if (diplotype.size() != 2)
        {
            return std::nullopt;
        }
        vector<long> repeatUnitsOnPath;
        for (const auto& path : diplotype)
        {
            repeatUnitsOnPath.push_back(std::count(path.nodeIds().begin(), path.nodeIds().end(), repeatNode));
        }
        if (repeatUnitsOnPath[0] == repeatUnitsOnPath[1])
        {
            return std::nullopt;
        }
        entryIndexOfPath = repeatUnitsOnPath[0] < repeatUnitsOnPath[1] ? vector<int> { 0, 1 } : vector<int> { 1, 0 };
    }

    vector<AlleleFlankingVariants> entries(hasOneDistinctAllele ? 1 : 2);
    entries[0].alleleSize = genotype.shortAlleleSizeInUnits();
    if (!hasOneDistinctAllele)
    {
        entries[1].alleleSize = genotype.longAlleleSizeInUnits();
    }

    vector<vector<ReadFlankObservations>> readsOfEachEntry(entries.size());
    std::map<FlankingVariant, PositionsShowingRefAllele> positionsShowingRefAlleleByVariant;
    for (size_t fragIndex = 0; fragIndex != fragAssignment.fragIds.size(); ++fragIndex)
    {
        const string& fragId = fragAssignment.fragIds[fragIndex];
        const int alignIndex = fragAssignment.alignIndexByFrag[fragIndex];
        const auto alignsIt = fragPathAlignsById.find(fragId);
        const auto fragIt = fragById.find(fragId);
        if (alignsIt == fragPathAlignsById.end() || fragIt == fragById.end() || alignIndex < 0
            || alignIndex >= static_cast<int>(alignsIt->second.size()))
        {
            continue;
        }

        int entryIndex = 0;
        if (!hasOneDistinctAllele)
        {
            if (!fragIdsPlacedOnOnePathByAlignmentScore.count(fragId))
            {
                continue;
            }
            // All of this read pair's alignments are on one path
            entryIndex = entryIndexOfPath[alignsIt->second[alignIndex].readAlign.pathIndex];
        }

        const FragPathAlign& fragPathAlign = alignsIt->second[alignIndex];
        const Frag& frag = fragIt->second;
        readsOfEachEntry[entryIndex].push_back(observeReadInFlanks(
            frag.read.bases(), *fragPathAlign.readAlign.align, leftFlank, rightFlank,
            positionsShowingRefAlleleByVariant));
        readsOfEachEntry[entryIndex].push_back(observeReadInFlanks(
            frag.mate.bases(), *fragPathAlign.mateAlign.align, leftFlank, rightFlank,
            positionsShowingRefAlleleByVariant));
    }

    // List the variants that enough reads of at least one allele show
    std::set<FlankingVariant> listedVariants;
    for (const auto& reads : readsOfEachEntry)
    {
        std::map<FlankingVariant, int> altReadsByVariant;
        for (const auto& read : reads)
        {
            for (const auto& variant : read.variants)
            {
                ++altReadsByVariant[variant];
            }
        }
        for (const auto& variantAndAltReads : altReadsByVariant)
        {
            if (variantAndAltReads.second >= kMinAltReadsToListFlankingVariant)
            {
                listedVariants.insert(variantAndAltReads.first);
            }
        }
    }

    for (size_t entryIndex = 0; entryIndex != entries.size(); ++entryIndex)
    {
        const auto& reads = readsOfEachEntry[entryIndex];
        entries[entryIndex].haplotypeRegion = computeHaplotypeRegion(reads, leftFlank, rightFlank);
        for (const auto& variant : listedVariants)
        {
            const PositionsShowingRefAllele& positions = positionsShowingRefAlleleByVariant.at(variant);
            FlankingVariantReadCounts counts;
            counts.variant = variant;
            for (const auto& read : reads)
            {
                if (read.variants.count(variant))
                {
                    ++counts.altReads;
                }
                else if (showsRefAllele(read, positions))
                {
                    ++counts.refReads;
                }
            }
            entries[entryIndex].variants.push_back(std::move(counts));
        }
    }

    return entries;
}

}  // namespace reviewer
}  // namespace ehunter
