//
// Consensus sequence generation for tandem repeat alleles.
//
// This module builds consensus sequences from "anchor reads" - reads that can be
// confidently assigned to a specific allele based on unambiguous alignment to a
// single haplotype path. The consensus represents the quality-weighted base
// composition at each position in the repeat region.
//
// Key concepts:
// - Anchor reads: Fragments that align to only one allele (all alignments have the
//   same pathIndex) - matching the visualization opacity criteria in LanePlot.cpp
// - PositionVotes: Quality-weighted base counts at each repeat position
// - AlleleConsensus: The final consensus sequence with 'N' for uncovered positions
//

#pragma once

#include <array>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "graphalign/GraphAlignment.hh"
#include "graphcore/Graph.hh"
#include "reviewer/Aligns.hh"
#include "reviewer/FragLenFilter.hh"
#include "reviewer/GenotypePaths.hh"

namespace ehunter
{
namespace reviewer
{

/// Quality-weighted base counts at a single position
struct PositionVotes
{
    std::array<double, 4> weights = { 0, 0, 0, 0 }; // A=0, C=1, G=2, T=3
    std::array<int, 4> counts = { 0, 0, 0, 0 };     // Raw counts per base: A=0, C=1, G=2, T=3
    int coverage = 0; // Number of reads covering this position

    /// Add a vote for a base with given quality weight
    void addVote(char base, double qualityWeight);

    /// Returns the most-voted base (A, C, G, or T)
    char consensusBase() const;

    /// Returns the count of reads supporting the consensus base
    int consensusBaseCount() const;

    /// Returns the sum of all quality weights
    double totalWeight() const;

    /// Returns true if any votes have been recorded
    bool hasVotes() const { return coverage > 0; }
};

/// Consensus sequence for one allele
struct AlleleConsensus
{
    int alleleIndex; // 0 or 1
    int repeatLength; // Length in base pairs
    std::vector<std::optional<PositionVotes>> positions;
    int anchorReadCount = 0; // Number of anchor reads used

    /// Returns the consensus sequence string with 'N' for unknown positions
    std::string toString() const;

    /// Returns a string where each character is a digit (0-9) representing the read support
    /// at that position. Coverage values >= 10 are capped to '9', coverage 0 is '0'.
    std::string toReadSupportString() const;

    /// Returns the count of non-N (known) positions
    int knownPositions() const;
};

/// Classification of whether a fragment is an anchor
struct AnchorClassification
{
    bool isAnchor = false;
    int assignedAllele = -1; // 0 or 1, or -1 if ambiguous
};

/// Result of consensus building
struct ConsensusResult
{
    std::vector<AlleleConsensus> alleleConsensuses; // One per allele
    int totalAnchors = 0;
    int totalFragments = 0;
};

/// Data extracted from a read's alignment to the repeat region
struct RepeatAlignedBases
{
    std::string sequence; // Bases from read that align to repeat
    std::vector<double> weights; // Quality weights for each base
    std::vector<int> repeatPositions; // Repeat coordinate for each base (handles deletions correctly)
    int repeatStartPos = 0; // Start position in repeat coordinates (0 = first base of repeat)
    bool valid = false; // Whether extraction was successful
};

/// Compute the total length of the repeat region in a haplotype path (in base pairs)
/// Returns 0 if no repeat nodes are found
int computeRepeatRegionLength(const graphtools::Path& haplotypePath);

/// Check if the repeat region in a haplotype path is contiguous
/// Returns true if all repeat nodes are adjacent (no non-repeat nodes between them)
/// Returns false if there are non-repeat nodes between repeat nodes
bool hasContiguousRepeatRegion(const graphtools::Path& haplotypePath);

/// Classify a fragment as an anchor read for consensus building
/// A fragment is an anchor if all its alignments have the same pathIndex (i.e., it
/// aligns to only one allele). This matches the visualization opacity criteria.
/// @param fragId Fragment ID for lookup
/// @param fragPathAlignsById Map of all fragment path alignments
/// @return AnchorClassification with isAnchor and assignedAllele
AnchorClassification classifyAsAnchor(const std::string& fragId, const FragPathAlignsById& fragPathAlignsById);

/// Extract the portion of a read that aligns to repeat nodes, along with position information
/// @param readWithAlign The read with its original graph alignment
/// @param pathAlign The read's alignment projected onto a haplotype path
/// @param haplotypePath The haplotype path the read is aligned to
/// @return RepeatAlignedBases with sequence, quality weights, and position information
RepeatAlignedBases extractRepeatAlignedBases(
    const ReadWithAlign& readWithAlign, const ReadPathAlign& pathAlign, const graphtools::Path& haplotypePath);

/// Build consensus sequence for one allele from anchor reads
/// @param alleleIndex The allele index (0 or 1)
/// @param haplotypePath The haplotype path for this allele
/// @param anchors Vector of (fragId, AnchorClassification) pairs for anchors assigned to this allele
/// @param fragById Map of all fragments by ID
/// @param fragPathAlignsById Map of fragment path alignments (to get projected alignments)
/// @param alignIndexByFragId Map from fragment ID to selected alignment index (from visualization)
/// @return AlleleConsensus with positions filled from anchor read votes
AlleleConsensus buildAlleleConsensus(
    int alleleIndex, const graphtools::Path& haplotypePath,
    const std::vector<std::pair<std::string, AnchorClassification>>& anchors, const FragById& fragById,
    const FragPathAlignsById& fragPathAlignsById, const std::unordered_map<std::string, int>& alignIndexByFragId);

/// Build consensus sequences for all alleles from anchor reads
/// @param diplotype The diplotype (vector of haplotype paths, one per allele)
/// @param fragById Map of all fragments by ID
/// @param fragPathAlignsById Map of fragment path alignments (filtered by fragment length)
/// @param fragAssignment Fragment assignment with selected alignment indices (same as visualization)
/// @return ConsensusResult with consensus sequences for each allele and statistics
ConsensusResult buildConsensusFromAnchors(
    const Diplotype& diplotype, const FragById& fragById, const FragPathAlignsById& fragPathAlignsById,
    const FragAssignment& fragAssignment);

} // namespace reviewer
} // namespace ehunter
