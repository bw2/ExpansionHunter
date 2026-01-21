#include "reviewer/ConsensusSequence.hh"

#include <unordered_map>

#include "spdlog/spdlog.h"

#include "graphalign/GraphAlignment.hh"
#include "graphalign/Operation.hh"
#include "graphcore/Graph.hh"
#include "reviewer/Projection.hh"

namespace ehunter
{
namespace reviewer
{

namespace
{

/// Convert a nucleotide base to an array index
/// A=0, C=1, G=2, T=3, other=-1
int baseToIndex(char base)
{
    switch (base)
    {
    case 'A':
    case 'a':
        return 0;
    case 'C':
    case 'c':
        return 1;
    case 'G':
    case 'g':
        return 2;
    case 'T':
    case 't':
        return 3;
    default:
        return -1;
    }
}

/// Convert an array index to a nucleotide base
/// 0=A, 1=C, 2=G, 3=T
char indexToBase(int index)
{
    switch (index)
    {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        return 'N';
    }
}

/// Check if all alignments for a fragment have the same pathIndex
/// Returns true if the fragment aligns to only one allele (anchor candidate)
/// Returns false if the fragment aligns to multiple alleles (ambiguous)
bool alignsToSinglePath(const std::vector<FragPathAlign>& fragPathAligns)
{
    if (fragPathAligns.empty())
    {
        return false;
    }

    int firstPathIndex = fragPathAligns.front().readAlign.pathIndex;
    for (const auto& fragPathAlign : fragPathAligns)
    {
        if (fragPathAlign.readAlign.pathIndex != firstPathIndex)
        {
            return false;
        }
    }

    return true;
}

/// Core consensus building logic shared by buildAlleleConsensus and buildAlleleConsensusHomozygous
/// @param haplotypePath The haplotype path for this allele
/// @param anchors Vector of (fragId, AnchorClassification) pairs
/// @param fragById Map of all fragments by ID
/// @param fragPathAlignsById Map of fragment path alignments
/// @param alignIndexByFragId Map from fragment ID to selected alignment index
/// @param alleleFilter Optional filter: if set, only process anchors with matching assignedAllele
/// @return AlleleConsensus (caller must set alleleIndex afterward)
AlleleConsensus buildConsensusCore(
    const graphtools::Path& haplotypePath,
    const std::vector<std::pair<std::string, AnchorClassification>>& anchors, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById, const std::unordered_map<std::string, int>& alignIndexByFragId,
    std::optional<int> alleleFilter)
{
    AlleleConsensus result;
    result.alleleIndex = 0; // Caller should set this appropriately
    result.anchorReadCount = 0;

    // Compute repeat length from haplotype path (sum of repeat node sequences)
    result.repeatLength = computeRepeatRegionLength(haplotypePath);

    // Initialize positions vector with nullopt for each position in the repeat region
    result.positions.resize(result.repeatLength);

    // Process each anchor
    for (const auto& anchorPair : anchors)
    {
        const std::string& fragId = anchorPair.first;
        const AnchorClassification& classification = anchorPair.second;

        // Apply allele filter if specified
        if (alleleFilter.has_value() && classification.assignedAllele != alleleFilter.value())
        {
            continue;
        }

        // Get the original fragment
        auto fragIt = fragById.find(fragId);
        if (fragIt == fragById.end())
        {
            continue;
        }
        const Frag& frag = fragIt->second;

        // Get the fragment's path alignments
        auto alignsIt = fragPathAlignsById.find(fragId);
        if (alignsIt == fragPathAlignsById.end())
        {
            continue;
        }

        // Use the same alignment index that the visualization uses
        auto alignIndexIt = alignIndexByFragId.find(fragId);
        if (alignIndexIt == alignIndexByFragId.end())
        {
            continue;
        }
        int alignIndex = alignIndexIt->second;
        if (alignIndex < 0 || static_cast<size_t>(alignIndex) >= alignsIt->second.size())
        {
            continue;
        }
        const FragPathAlign* fragAlign = &alignsIt->second[alignIndex];

        // Extract repeat-aligned bases from the read
        RepeatAlignedBases readBases = extractRepeatAlignedBases(frag.read, fragAlign->readAlign, haplotypePath);

        // Add votes from the read using position-accurate repeat coordinates
        if (readBases.valid && !readBases.repeatPositions.empty())
        {
            for (size_t i = 0; i < readBases.sequence.size(); ++i)
            {
                int pos = readBases.repeatPositions[i];
                if (pos >= 0 && pos < result.repeatLength)
                {
                    if (!result.positions[pos].has_value())
                    {
                        result.positions[pos] = PositionVotes();
                    }
                    result.positions[pos]->addVote(
                        readBases.sequence[i], (i < readBases.weights.size()) ? readBases.weights[i] : 0.99);
                }
            }
        }

        // Extract repeat-aligned bases from the mate
        RepeatAlignedBases mateBases = extractRepeatAlignedBases(frag.mate, fragAlign->mateAlign, haplotypePath);

        // Add votes from the mate using position-accurate repeat coordinates
        if (mateBases.valid && !mateBases.repeatPositions.empty())
        {
            for (size_t i = 0; i < mateBases.sequence.size(); ++i)
            {
                int pos = mateBases.repeatPositions[i];
                if (pos >= 0 && pos < result.repeatLength)
                {
                    if (!result.positions[pos].has_value())
                    {
                        result.positions[pos] = PositionVotes();
                    }
                    result.positions[pos]->addVote(
                        mateBases.sequence[i], (i < mateBases.weights.size()) ? mateBases.weights[i] : 0.99);
                }
            }
        }

        // Count this anchor if we got any valid bases from it
        if ((readBases.valid && !readBases.sequence.empty()) || (mateBases.valid && !mateBases.sequence.empty()))
        {
            ++result.anchorReadCount;
        }
    }

    return result;
}

} // anonymous namespace

void PositionVotes::addVote(char base, double qualityWeight)
{
    int index = baseToIndex(base);
    if (index >= 0)
    {
        weights[index] += qualityWeight;
        counts[index] += 1;
        ++coverage;
    }
}

char PositionVotes::consensusBase() const
{
    if (coverage == 0)
    {
        return 'N';
    }

    int bestIndex = 0;
    double bestWeight = weights[0];

    for (int i = 1; i < 4; ++i)
    {
        if (weights[i] > bestWeight)
        {
            bestWeight = weights[i];
            bestIndex = i;
        }
    }

    return indexToBase(bestIndex);
}

int PositionVotes::consensusBaseCount() const
{
    if (coverage == 0)
    {
        return 0;
    }

    // Find the index of the consensus base (same logic as consensusBase)
    int bestIndex = 0;
    double bestWeight = weights[0];

    for (int i = 1; i < 4; ++i)
    {
        if (weights[i] > bestWeight)
        {
            bestWeight = weights[i];
            bestIndex = i;
        }
    }

    return counts[bestIndex];
}

double PositionVotes::totalWeight() const { return weights[0] + weights[1] + weights[2] + weights[3]; }

std::string AlleleConsensus::toString() const
{
    std::string result;
    result.reserve(positions.size());

    for (const auto& position : positions)
    {
        if (position.has_value() && position->hasVotes())
        {
            result += position->consensusBase();
        }
        else
        {
            result += 'N';
        }
    }

    return result;
}

std::string AlleleConsensus::toReadSupportString() const
{
    std::string result;
    result.reserve(positions.size());

    for (const auto& position : positions)
    {
        if (position.has_value() && position->hasVotes())
        {
            // Count only reads that match the consensus base, capped at 9 for single-digit display
            int supportCount = std::min(position->consensusBaseCount(), 9);
            result += static_cast<char>('0' + supportCount);
        }
        else
        {
            result += '0';
        }
    }

    return result;
}

int AlleleConsensus::knownPositions() const
{
    int count = 0;
    for (const auto& position : positions)
    {
        if (position.has_value() && position->hasVotes())
        {
            ++count;
        }
    }
    return count;
}

int computeRepeatRegionLength(const graphtools::Path& haplotypePath)
{
    const graphtools::Graph* graph = haplotypePath.graphRawPtr();
    const auto& nodeIds = haplotypePath.nodeIds();
    int totalLength = 0;

    for (graphtools::NodeId nodeId : nodeIds)
    {
        // Check if this is a repeat node (has self-loop edge)
        if (graph->hasEdge(nodeId, nodeId))
        {
            totalLength += static_cast<int>(graph->nodeSeq(nodeId).length());
        }
    }

    return totalLength;
}

bool hasContiguousRepeatRegion(const graphtools::Path& haplotypePath)
{
    const graphtools::Graph* graph = haplotypePath.graphRawPtr();
    const auto& nodeIds = haplotypePath.nodeIds();

    // Track whether we've seen repeat nodes and whether we've left the repeat region
    bool inRepeatRegion = false;
    bool exitedRepeatRegion = false;

    for (graphtools::NodeId nodeId : nodeIds)
    {
        bool isRepeatNode = graph->hasEdge(nodeId, nodeId);

        if (isRepeatNode)
        {
            // If we see a repeat node after having exited the repeat region,
            // then the repeat region is non-contiguous
            if (exitedRepeatRegion)
            {
                return false;
            }
            inRepeatRegion = true;
        }
        else if (inRepeatRegion)
        {
            // We're in a non-repeat node after having seen repeat nodes
            exitedRepeatRegion = true;
        }
    }

    return true;
}

AnchorClassification classifyAsAnchor(const std::string& fragId, const FragPathAlignsById& fragPathAlignsById)
{
    AnchorClassification result;
    result.isAnchor = false;
    result.assignedAllele = -1;

    // Look up all alignments for this fragment
    auto alignsIt = fragPathAlignsById.find(fragId);
    if (alignsIt == fragPathAlignsById.end() || alignsIt->second.empty())
    {
        return result;
    }

    const std::vector<FragPathAlign>& fragPathAligns = alignsIt->second;

    // A fragment is an anchor if all its alignments have the same pathIndex
    // This matches the visualization opacity logic in LanePlot.cpp
    if (alignsToSinglePath(fragPathAligns))
    {
        result.isAnchor = true;
        result.assignedAllele = fragPathAligns.front().readAlign.pathIndex;
    }

    return result;
}

RepeatAlignedBases extractRepeatAlignedBases(
    const ReadWithAlign& readWithAlign, const ReadPathAlign& pathAlign, const graphtools::Path& haplotypePath)
{
    RepeatAlignedBases result;
    result.valid = false;
    result.repeatStartPos = -1;

    const graphtools::Graph* graph = haplotypePath.graphRawPtr();
    const auto& hapNodeIds = haplotypePath.nodeIds();

    // Calculate the start and end positions of the repeat region on the linearized haplotype path
    int repeatRegionStart = -1;
    int repeatRegionEnd = -1;
    int currentPos = 0;

    for (graphtools::NodeId nodeId : hapNodeIds)
    {
        int nodeLen = static_cast<int>(graph->nodeSeq(nodeId).length());
        bool isRepeatNode = graph->hasEdge(nodeId, nodeId);

        if (isRepeatNode)
        {
            if (repeatRegionStart < 0)
            {
                repeatRegionStart = currentPos;
            }
            repeatRegionEnd = currentPos + nodeLen;
        }

        currentPos += nodeLen;
    }

    // If no repeat region in the path, return invalid result
    if (repeatRegionStart < 0 || repeatRegionEnd < 0)
    {
        return result;
    }

    // Check if the read overlaps the repeat region at all
    // pathAlign.begin and pathAlign.end are positions on the linearized haplotype
    if (pathAlign.end <= repeatRegionStart || pathAlign.begin >= repeatRegionEnd)
    {
        // Read doesn't overlap repeat region
        return result;
    }

    // Now we need to extract the bases from the read that align to the repeat region
    // We'll walk through the projected alignment and track which bases map to repeat nodes

    const graphtools::GraphAlignment& projectedAlign = *pathAlign.align;
    const std::string& readSequence = readWithAlign.bases();
    const auto& alignNodeIds = projectedAlign.path().nodeIds();
    const auto& alignments = projectedAlign.alignments();

    // Default quality weight (Q20 equivalent, as base qualities aren't readily available)
    constexpr double defaultQualityWeight = 0.99;

    // Track our position on the haplotype path and in the read sequence
    int hapPathPos = pathAlign.begin; // Position on linearized haplotype
    int queryPos = 0; // Position in read sequence

    // Calculate the prefix length on the haplotype up to where alignment starts
    // This accounts for the startIndexOnPath
    int prefixLen = 0;
    for (int i = 0; i < pathAlign.startIndexOnPath; ++i)
    {
        prefixLen += static_cast<int>(haplotypePath.graphRawPtr()->nodeSeq(hapNodeIds[i]).length());
    }

    // The alignment's start position is within the first aligned node
    hapPathPos = prefixLen + static_cast<int>(projectedAlign.path().startPosition());

    // Process each node in the projected alignment
    for (size_t nodeIdx = 0; nodeIdx < alignNodeIds.size(); ++nodeIdx)
    {
        graphtools::NodeId nodeId = alignNodeIds[nodeIdx];
        const graphtools::Alignment& nodeAlignment = alignments[nodeIdx];

        // Check if this node is a repeat node
        bool isRepeatNode = graph->hasEdge(nodeId, nodeId);

        // Process each operation in this node's alignment
        for (const auto& op : nodeAlignment.operations())
        {
            uint32_t opLen = op.length();
            graphtools::OperationType opType = op.type();

            switch (opType)
            {
            case graphtools::OperationType::kMatch:
            case graphtools::OperationType::kMismatch:
            {
                // Both match and mismatch consume both query and reference bases
                for (uint32_t i = 0; i < opLen; ++i)
                {
                    // Check if this position is within the repeat region
                    if (isRepeatNode && hapPathPos >= repeatRegionStart && hapPathPos < repeatRegionEnd)
                    {
                        // This base is in the repeat region
                        int repeatCoord = hapPathPos - repeatRegionStart;

                        // Extract the base from the read sequence
                        if (queryPos < static_cast<int>(readSequence.size()))
                        {
                            result.sequence += readSequence[queryPos];
                            result.weights.push_back(defaultQualityWeight);
                            result.repeatPositions.push_back(repeatCoord);
                        }
                    }

                    ++queryPos;
                    ++hapPathPos;
                }
                break;
            }
            case graphtools::OperationType::kInsertionToRef:
            {
                // Insertion: query bases that don't map to reference
                // These bases don't have a position on the reference, so we skip them
                // for consensus building purposes (they would mess up position mapping)
                queryPos += static_cast<int>(opLen);
                break;
            }
            case graphtools::OperationType::kDeletionFromRef:
            {
                // Deletion: reference bases that aren't in the query
                // We advance the reference position but don't consume query bases
                // For consensus, we could add 'N' or '-' for these positions,
                // but we'll just skip them and let other reads fill in
                hapPathPos += static_cast<int>(opLen);
                break;
            }
            case graphtools::OperationType::kSoftclip:
            {
                // Soft clip: query bases that aren't aligned
                // Don't consume reference positions, just query positions
                queryPos += static_cast<int>(opLen);
                break;
            }
            case graphtools::OperationType::kMissingBases:
            {
                // Missing bases: skip both query and reference
                queryPos += static_cast<int>(opLen);
                hapPathPos += static_cast<int>(opLen);
                break;
            }
            }
        }
    }

    // If we extracted any bases, mark as valid and set repeatStartPos from first position
    if (!result.sequence.empty() && !result.repeatPositions.empty())
    {
        result.repeatStartPos = result.repeatPositions.front();
        result.valid = true;
    }

    return result;
}

AlleleConsensus buildAlleleConsensus(
    int alleleIndex, const graphtools::Path& haplotypePath,
    const std::vector<std::pair<std::string, AnchorClassification>>& anchors, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById,
    const std::unordered_map<std::string, int>& alignIndexByFragId)
{
    AlleleConsensus result = buildConsensusCore(
        haplotypePath, anchors, fragById, fragPathAlignsById, alignIndexByFragId, alleleIndex);
    result.alleleIndex = alleleIndex;
    return result;
}

AlleleConsensus buildAlleleConsensusHomozygous(
    const graphtools::Path& haplotypePath,
    const std::vector<std::pair<std::string, AnchorClassification>>& anchors, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById,
    const std::unordered_map<std::string, int>& alignIndexByFragId)
{
    // For homozygous: no allele filter (std::nullopt), include all anchors
    AlleleConsensus result = buildConsensusCore(
        haplotypePath, anchors, fragById, fragPathAlignsById, alignIndexByFragId, std::nullopt);
    result.alleleIndex = 0;
    return result;
}

ConsensusResult buildConsensusFromAnchors(
    const Diplotype& diplotype, const FragById& fragById, const FragPathAlignsById& fragPathAlignsById,
    const FragAssignment& fragAssignment)
{
    ConsensusResult result;
    result.totalFragments = static_cast<int>(fragPathAlignsById.size());
    result.totalAnchors = 0;

    // Check for non-contiguous repeat regions - consensus building doesn't support these yet
    for (size_t i = 0; i < diplotype.size(); ++i)
    {
        if (!hasContiguousRepeatRegion(diplotype[i]))
        {
            spdlog::debug("Consensus building: Skipping locus with non-contiguous repeat region in allele {}", i);
            return result;  // Return empty result
        }
    }

    // Build a map from fragId -> alignIndex for efficient lookup
    // This uses the same alignment selection as the visualization
    std::unordered_map<std::string, int> alignIndexByFragId;
    for (size_t i = 0; i < fragAssignment.fragIds.size(); ++i)
    {
        alignIndexByFragId[fragAssignment.fragIds[i]] = fragAssignment.alignIndexByFrag[i];
    }

    // Check if this is effectively a homozygous genotype (both alleles have identical paths)
    // For homozygous, we should only produce one consensus sequence with all anchors contributing.
    // We require path identity (same nodes, same positions), not just equal repeat length,
    // to avoid merging equal-length but sequence-divergent alleles (e.g., different SNPs).
    bool isHomozygous = false;
    if (diplotype.size() == 2)
    {
        isHomozygous = (diplotype[0] == diplotype[1]);
    }

    // Classify all fragments and collect anchors
    std::vector<std::pair<std::string, AnchorClassification>> allAnchors;

    for (const auto& fragAlignsPair : fragPathAlignsById)
    {
        const std::string& fragId = fragAlignsPair.first;

        AnchorClassification classification;

        if (isHomozygous)
        {
            // For homozygous genotypes (identical paths), all fragments with valid
            // alignments are anchors since both alleles are the same sequence.
            // Normal classification would fail because fragments align to both
            // pathIndex 0 and 1 equally well.
            classification.isAnchor = true;
            classification.assignedAllele = 0;
        }
        else
        {
            classification = classifyAsAnchor(fragId, fragPathAlignsById);
        }

        if (classification.isAnchor && classification.assignedAllele >= 0)
        {
            // Only count as anchor if fragment is in fragAssignment (has a selected alignment)
            if (alignIndexByFragId.find(fragId) != alignIndexByFragId.end())
            {
                allAnchors.emplace_back(fragId, classification);
                ++result.totalAnchors;
            }
        }
    }

    if (isHomozygous)
    {
        // For homozygous genotypes, produce one consensus using all anchors
        // Anchors may have assignedAllele 0 or 1 (both are same repeat length)
        AlleleConsensus alleleConsensus = buildAlleleConsensusHomozygous(
            diplotype[0], allAnchors, fragById, fragPathAlignsById, alignIndexByFragId);

        result.alleleConsensuses.push_back(std::move(alleleConsensus));
    }
    else
    {
        // For heterozygous genotypes, separate anchors by allele
        std::vector<std::vector<std::pair<std::string, AnchorClassification>>> anchorsByAllele(diplotype.size());

        for (const auto& anchorPair : allAnchors)
        {
            int alleleIdx = anchorPair.second.assignedAllele;
            if (alleleIdx >= 0 && static_cast<size_t>(alleleIdx) < diplotype.size())
            {
                anchorsByAllele[alleleIdx].push_back(anchorPair);
            }
        }

        // Build consensus for each allele
        for (size_t alleleIdx = 0; alleleIdx < diplotype.size(); ++alleleIdx)
        {
            AlleleConsensus alleleConsensus = buildAlleleConsensus(
                static_cast<int>(alleleIdx), diplotype[alleleIdx], anchorsByAllele[alleleIdx], fragById,
                fragPathAlignsById, alignIndexByFragId);

            result.alleleConsensuses.push_back(std::move(alleleConsensus));
        }
    }

    return result;
}

} // namespace reviewer
} // namespace ehunter
