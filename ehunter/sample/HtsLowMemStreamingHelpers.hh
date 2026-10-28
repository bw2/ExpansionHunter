#pragma once

#include <algorithm>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "core/CountTable.hh"
#include "core/GenomicRegion.hh"
#include "core/Parameters.hh"
#include "core/Read.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/IterativeJsonWriter.hh"
#include "io/IterativeVcfWriter.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"


namespace ehunter
{

struct FastReadAnalysisResult {
    int repeat_sequence_size_in_base_pairs = 0;
    bool overlaps_repeats = false;
    bool is_spanning_read = false;
    bool soft_clipped_bases_contain_repetitive_sequence = false;
    int inserted_bases_within_repeat = 0;  // CIGAR insertions whose position falls within the locus
    int deleted_bases_within_repeat = 0;   // CIGAR deletions overlapping the locus
    // repeat-region read bases matching a perfect consecutive repeat sequence of the motif (ReadRepeatPurity)
    int matched_bases_within_repeat = 0;
    int repeat_read_bases = 0;             // total repeat-region read bases (ReadRepeatPurity denominator)
    // The read's repeat tract in read-sequence coordinates: the read bases lying between the two locus
    // edges. Used to build per-allele consensus sequences. Set when both edges can be placed in read
    // coordinates and the tract is non-empty and within the read sequence -- an edge is placed by an
    // aligned match operation, by a whole-motif insertion abutting the edge, or by a deletion containing
    // the edge. repeat_tract_length < 0 means "not usable". For a read that does not span the locus that
    // is simply an edge the read never reaches; for a spanning read -- the only kind the consensus uses
    // -- both edges are always placed, so the one residual case is a zero-length tract, i.e. an alignment
    // that deletes the whole locus. See processRead() for the reasoning behind each case.
    int repeat_tract_start = -1;
    int repeat_tract_length = -1;
};


// Check if region A is to the left of region B
bool isAToTheLeftOfB(int32_t contigIdA, int64_t posA, int32_t contigIdB, int64_t posB, bool allowEquals = false);

// Check if region A fully contains region B
bool doesAcontainB(int32_t contigIdA, int64_t startA, int64_t endA, int32_t contigIdB, int64_t startB, int64_t endB);


// 64-bit FNV-1a hash, used to derive a per-locus reservoir-sampling seed from the (stable) locusId so the
// sample a locus retains is independent of catalog subsetting / thread count, keeping output byte-identical
// across --threads.
inline uint64_t fnv1a64(const std::string& s)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) { h ^= c; h *= 0x100000001b3ULL; }
    return h;
}

// Says whether the alignment pipeline flagged the read as a duplicate (SAM flag 0x400) or as failing quality checks
// (0x200). Once a locus reaches its --max-depth cap, read pairs with such a read give way to pairs without one (see
// LocusCache::offerReadPair).
inline bool isFlaggedAsDuplicateOrFailingQualityChecks(const LinearAlignmentStats& alignmentStats)
{
    return alignmentStats.isFlaggedAsDuplicate || alignmentStats.isFlaggedAsFailingQualityChecks;
}

// The read pairs collected for one locus in optimized-streaming mode. Exposed here so it can be unit tested.
struct LocusCache {
    // maxReadPairs_ is the per-locus reservoir size, derived from --max-depth so the retained reads bound the
    // average base-level depth over the locus window (reference repeat regions + flanks) to maxDepth. For a
    // window of width W and read length L, ~D*(W+L)/L reads achieve depth D — a read overlaps the window if it
    // starts anywhere in a span of W+L bases — so the pair cap is that read count halved (one reservoir entry
    // is a read pair). The cap therefore scales with locus width, bounding narrow and wide loci to the same
    // depth instead of the same absolute read count. maxDepth <= 0 disables the cap (maxReadPairs_ == 0).
    LocusCache(int locusIndex, const GenomicRegion& locusInterval, const std::string& locusId,
               int maxDepth, int typicalReadLength)
    : locusIndex(locusIndex), locusInterval(locusInterval), rngState_(fnv1a64(locusId))
    {
        if (maxDepth <= 0) {
            maxReadPairs_ = 0;
        } else {
            const int64_t readLength = std::max<int64_t>(1, typicalReadLength);
            const int64_t capReads = static_cast<int64_t>(maxDepth) * (locusInterval.length() + readLength) / readLength;
            maxReadPairs_ = static_cast<std::size_t>(std::max<int64_t>(1, capReads / 2));
        }
    }

    // Stores all reads + mates that will be analyzed for a given locus.
    int locusIndex;

    // the interval spanned by locus reference regions, including flanks
    GenomicRegion locusInterval;

    //a vector of pointers to FullReadPair objects
    std::vector<std::shared_ptr<FullReadPair>> readPairs;

    // fragment ids already offered to this locus, used to avoid double-counting a far-away pair whose
    // two ends are both contained in this locus window (each end is streamed in a separate iteration).
    // Bounded by the reservoir: past the cap, admitFragment() only adds pairs that take a flagged pair's
    // slot, so this set never holds more than ~2 x maxReadPairs entries (see admitFragment).
    std::unordered_set<FragmentId> seenFragmentIds;

    // Decide whether `fragmentId` should be offered to this locus, applying per-locus duplicate suppression.
    // Returns true (admit) for a fragment not seen before. `isFlaggedPair` says whether either read is flagged as a
    // duplicate or as failing quality checks. Once the reservoir is full (readPairs at maxReadPairs), a fragment
    // seen earlier is still rejected (earlier builds admitted it again, so over-cap loci can keep a slightly
    // different sample than before, even without flagged reads), but new ids are tracked only for the unflagged
    // pairs offered while flagged slots remain, at most one per flagged slot. Tracking every fragment id would
    // make seenFragmentIds grow to O(distinct reads) at the very pathological loci the cap exists to bound. So a
    // far-away pair whose two ends both land in this window, and whose first end arrives past the cap untracked,
    // may be offered twice, which only perturbs already-sampled over-cap loci and stays deterministic across
    // --threads (the cap is reached at the same offer under any thread count). When the cap is disabled
    // (maxReadPairs == 0) this is exactly the previous unbounded dedup.
    bool admitFragment(const FragmentId& fragmentId, bool isFlaggedPair, std::size_t maxReadPairs)
    {
        if (maxReadPairs == 0 || readPairs.size() < maxReadPairs)
        {
            return seenFragmentIds.insert(fragmentId).second;
        }
        if (seenFragmentIds.count(fragmentId) != 0)
        {
            return false;
        }
        if (!isFlaggedPair && !slotsHoldingFlaggedPairs_.empty())
        {
            seenFragmentIds.insert(fragmentId);
        }
        return true;
    }

    // Add a (deduped) read pair to this locus, keeping at most `maxReadPairs` entries. maxReadPairs == 0 means no
    // cap. While the locus is below the cap this is a plain push_back, so loci that never exceed the cap end up
    // with exactly the same readPairs (in the same order) as without the cap — only over-cap loci (e.g.
    // centromeric/satellite pileups) get sampled. Past the cap, a pair with neither read flagged as a duplicate or
    // as failing quality checks takes the slot of a random flagged pair while any remain, so flagged pairs are kept
    // only where there aren't enough unflagged ones. Within each group the retained pairs are a uniform sample
    // (Algorithm R): an unflagged pair offered once no flagged slots remain replaces a random slot with probability
    // maxReadPairs / (unflagged pairs offered), and a flagged pair replaces a random flagged slot with probability
    // (flagged slots) / (flagged pairs offered). With no flagged pairs this is plain Algorithm R. The replacement
    // choices use a per-locus RNG seeded from the locusId, and reads are offered in coordinate-stream order (the
    // same order under any --threads), so the retained sample is deterministic and identical across thread counts.
    void offerReadPair(const std::shared_ptr<FullReadPair>& readPair, bool isFlaggedPair, std::size_t maxReadPairs)
    {
        // number of pairs of this pair's group offered so far (before the current one)
        const uint64_t offeredCountInGroup
            = isFlaggedPair ? offeredFlaggedCount_ : offeredCount_ - offeredFlaggedCount_;
        ++offeredCount_;
        if (isFlaggedPair)
        {
            ++offeredFlaggedCount_;
        }

        if (maxReadPairs == 0 || readPairs.size() < maxReadPairs)
        {
            if (isFlaggedPair)
            {
                slotsHoldingFlaggedPairs_.push_back(readPairs.size());
            }
            readPairs.push_back(readPair);
        }
        else if (!isFlaggedPair && !slotsHoldingFlaggedPairs_.empty())
        {
            const uint64_t k = nextRandom(slotsHoldingFlaggedPairs_.size());
            readPairs[slotsHoldingFlaggedPairs_[k]] = readPair;
            slotsHoldingFlaggedPairs_[k] = slotsHoldingFlaggedPairs_.back();
            slotsHoldingFlaggedPairs_.pop_back();
        }
        else
        {
            // The group's slots hold a uniform sample of the group's pairs offered so far: every slot once no flagged
            // slots remain (until then no unflagged pair was dropped), or the remaining flagged slots.
            const uint64_t groupSlotCount = isFlaggedPair ? slotsHoldingFlaggedPairs_.size() : maxReadPairs;
            const uint64_t j = nextRandom(offeredCountInGroup + 1);  // uniform in [0, offeredCountInGroup]
            if (j < groupSlotCount)
            {
                readPairs[isFlaggedPair ? slotsHoldingFlaggedPairs_[j] : j] = readPair;
            }
        }
    }

    // Per-locus reservoir size in read pairs; 0 = unlimited. Derived from --max-depth (see constructor).
    std::size_t maxReadPairs() const { return maxReadPairs_; }

    // True if this locus exceeded its reservoir cap — more read pairs were offered than maxReadPairs_, so the
    // retained set is a reservoir sample rather than every read. False when the cap is disabled
    // (maxReadPairs_ == 0) or the locus stayed at/under the cap (its read set is complete, byte-identical to no cap).
    bool reservoirSampled() const { return maxReadPairs_ != 0 && offeredCount_ > maxReadPairs_; }

    // Total read pairs offered to this locus (admitted + dropped). Exceeds maxReadPairs() exactly when
    // reservoir sampling was triggered.
    uint64_t offeredCount() const { return offeredCount_; }

private:
    // number of read pairs offered so far (before the current one); drives the reservoir replacement odds
    uint64_t offeredCount_ = 0;
    // how many of them had a read flagged as a duplicate or as failing quality checks
    uint64_t offeredFlaggedCount_ = 0;
    // indexes into readPairs of the slots that hold such flagged pairs
    std::vector<std::size_t> slotsHoldingFlaggedPairs_;
    uint64_t rngState_;
    std::size_t maxReadPairs_ = 0;

    // splitmix64, returning a value uniform in [0, bound). The modulo bias is negligible here (bound is at
    // most a few million while the generator is full 64-bit), and avoiding <random> keeps the per-locus
    // state at 8 bytes — important since one LocusCache exists per active locus.
    uint64_t nextRandom(uint64_t bound)
    {
        rngState_ += 0x9E3779B97F4A7C15ULL;
        uint64_t z = rngState_;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        z = z ^ (z >> 31);
        return z % bound;
    }
};


// Check if a RepeatGenotype is homozygous reference
bool isRepeatGenotypeHomRef(const RepeatGenotype& genotype, int referenceSizeInUnits);

// Decide whether a locus should be excluded from the VCF/JSON output given the --skip-hom-ref and
// --skip-missing-genotypes flags. The locus is filtered out only if it has at least one variant and EVERY
// variant is "skippable": homozygous reference (when skipHomRef) or missing/no-genotype (when skipMissing).
// A locus with any non-reference call, or any unclassifiable variant, is always kept. With skipMissing=false
// this reproduces the original --skip-hom-ref behavior (skip iff all variants are called homozygous ref);
// with both flags set, loci whose variants are any mixture of hom-ref and missing are also filtered out.
bool shouldFilterLocus(
    const LocusSpecification& locusSpec, const LocusFindings& locusFindings, bool skipHomRef, bool skipMissing);

// Analyze a single read against a repeat locus given in 0-based half-open coordinates
// (start inclusive, end exclusive). Exposed here so it can be unit tested directly.
FastReadAnalysisResult processRead(
    const FullRead& read, int64_t locus_start_0based, int64_t locus_end_0based, const std::string& locus_motif);

// Says whether at least two of the reads (both mates of every pair, wherever they're aligned) are in-repeat reads for
// the motif, i.e. evidence of an allele longer than a read, which processLocusFast can't size. Exposed for unit tests.
bool hasEnoughInRepeatReadsForFullGenotyping(
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, const std::string& locusMotif);

// Says whether some read-length window of the reference sequence is nearly pure repeat of the motif (weighted purity
// at least 0.75), in which case in-repeat reads at the locus aren't evidence of an expansion. Exposed for unit tests.
bool hasReadLengthStretchOfInRepeatSequence(
    const std::string& referenceSequence, const std::string& locusMotif, int readLength);

// Converts read counts keyed by repeat size in bp into a table keyed by repeat size in whole units (bp / motif size,
// rounded down, the way processLocusFast sizes alleles). Counts of bp sizes that round to the same number of units are
// added up. Exposed for unit tests.
CountTable countReadsByRepeatSizeInUnits(const std::map<int, int>& readCountsByRepeatSizeInBasePairs, int motifSize);

bool processLocusFast(
	const ProgramParameters& params, Reference& reference, LocusDescription& locusDescription,
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, bool reservoirSampled, int typicalReadLength,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter);

// When --output-motif-composition is specified, computes the motif composition of the locus's repeat from its cached reads and
// attaches it to the repeat's findings. Applies only to loci with a single repeat variant (small variants may sit
// alongside it) whose motif is between 2 bp and a third of typicalReadLength; other loci are left unchanged. In
// loci-with-known-motifs mode, loci whose catalog record has no KnownMotifs are also left unchanged. With no reads,
// all-loci and loci-with-known-motifs modes attach an empty composition (no reference access), so every eligible locus
// written to the JSON carries the field. offtargetReadPairs are the locus's off-target read pairs (see
// computeMotifComposition); only loci genotyped from seeking-mode reads have them.
void addMotifComposition(
    const ProgramParameters& params, Reference& reference, const LocusSpecification& locusSpec,
    LocusFindings& locusFindings, const std::vector<std::shared_ptr<FullReadPair>>& readPairs, int typicalReadLength,
    const std::vector<std::shared_ptr<FullReadPair>>& offtargetReadPairs = {});

// Emit a no-call record for a locus with zero coverage, identical to what seeking/streaming mode produce,
// without running genotyping. For single-region loci it builds a flankless stub LocusSpecification
// (extendFlanks=false) so the reference flanks are not read and the only reference access is the single
// left-flank base the VCF REF column requires. Multi-region loci instead use the flank-extending decode
// (extendFlanks=true), which reads ~regionExtensionLength reference flank bases per side; see the
// implementation comment in HtsLowMemStreamingHelpers.cpp for why. It then writes empty findings whose
// per-variant values mirror what RepeatAnalyzer/SmallVariantAnalyzer return on low depth (empty read
// counts, no genotype, LowDepth filter). The stub spec is local and freed as soon as the record is
// written, preserving the streaming mode's low-memory profile.
// Returns true if a no-call record was written, false if the locus was suppressed (--skip-missing-genotypes).
bool writeZeroCoverageRecord(
    const ProgramParameters& params, Reference& reference, const LocusDescription& locusDescription,
    int typicalReadLength, IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter);
}
