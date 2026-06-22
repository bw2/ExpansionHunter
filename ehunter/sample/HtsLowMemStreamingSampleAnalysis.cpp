/**

TODO: profiling using heaptrack and time profiling using gprof and perf or callgrind

TODO: stop after all loci have been processed and all mates have been collected
TODO: downsample reads?

TODO: implement off-target regions (using stripy to infer off-targets?)
TODO: implement bamlet writing?
TODO: speed up alignment:  lower coverage, narrower reference region, only process reads with signal
TODO: add multithreading
TODO: add reviewer

Core data structures:

mateExtractor cache
UnpairedReadsCache: a hashset of FullRead objects that have been read but whose mate has not yet been found
LocusCache: a hashmap of LocusCache objects, each of which contains a vector of FullReadPair objects
    that will be analyzed for a given locus

*/

#include "sample/HtsStreamingSampleAnalysis.hh"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <exception>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "spdlog/spdlog.h"
#include <boost/optional.hpp>
#include <boost/timer/progress_display.hpp>

#include "core/HtsHelpers.hh"
#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "io/BamletWriter.hh"
#include "io/CatalogLoading.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/IterativeJsonWriter.hh"
#include "io/IterativeVcfWriter.hh"
#include "io/ParameterLoading.hh"
#include "io/StreamingOutputMerge.hh"
#include "io/SampleStats.hh"
#include "io/StringUtils.hh"
#include "io/VcfWriter.hh"
#include "locus/VariantFindings.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusSpecification.hh"
#include "locus/LocusAnalyzerUtil.hh"
#include "sample/AnalyzerFinder.hh"
#include "sample/GenomeQueryCollection.hh"
#include "sample/HtsFileStreamer.hh"
#include "sample/HtsLowMemStreamingHelpers.hh"
#include "sample/HtsStreamingReadPairQueue.hh"
#include "sample/MateExtractor.hh"

#include "locus/RepeatAnalyzer.hh" // TODO drop this
#include "genotyping/AlignMatrix.hh" // TODO drop this

using ehunter::locus::initializeLocusAnalyzers;
using ehunter::locus::LocusAnalyzer;
using std::string;
using std::vector;
using std::unordered_map;
using std::shared_ptr;

namespace ehunter
{


FullRead decodeRead(const htshelpers::HtsFileStreamer& readStreamer) {
    // fully parse the current read
    LinearAlignmentStats readAlignmentStats;
    Read read = readStreamer.decodeRead(readAlignmentStats);
    return {std::move(read), std::move(readAlignmentStats)};
}


// 64-bit FNV-1a hash, used to derive a per-locus reservoir-sampling seed from the (stable) locusId so the
// sample a locus retains is independent of catalog subsetting / thread count, keeping output byte-identical
// across --threads.
inline uint64_t fnv1a64(const std::string& s)
{
    uint64_t h = 0xcbf29ce484222325ULL;
    for (unsigned char c : s) { h ^= c; h *= 0x100000001b3ULL; }
    return h;
}

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
    std::vector<shared_ptr<FullReadPair>> readPairs;

    // fragment ids already added to readPairs, used to avoid double-counting a far-away pair whose
    // two ends are both contained in this locus window (each end is streamed in a separate iteration).
    // Bounded by the reservoir: once readPairs reaches the cap, admitFragment() stops inserting here, so
    // this set never holds more than ~maxReadPairs entries (see admitFragment).
    std::unordered_set<FragmentId> seenFragmentIds;

    // Decide whether `fragmentId` should be offered to this locus, applying per-locus duplicate suppression.
    // Returns true (admit) for a fragment not seen before. Once the reservoir is full (readPairs at maxReadPairs)
    // we stop tracking fragment ids — otherwise seenFragmentIds would grow to O(distinct reads) at the very
    // pathological loci the cap exists to bound. Past that point a far-away pair whose two ends both land in
    // this window may be offered twice, which only perturbs already-sampled over-cap loci and stays
    // deterministic across --threads (the cap is reached at the same offer under any thread count). When the
    // cap is disabled (maxReadPairs == 0) this is exactly the previous unbounded dedup.
    bool admitFragment(const FragmentId& fragmentId, std::size_t maxReadPairs)
    {
        if (maxReadPairs != 0 && readPairs.size() >= maxReadPairs)
        {
            return true;
        }
        return seenFragmentIds.insert(fragmentId).second;
    }

    // Add a (deduped) read pair to this locus, keeping at most `maxReadPairs` entries via reservoir sampling
    // (Algorithm R). maxReadPairs == 0 means no cap. While the locus is below the cap this is a plain
    // push_back, so loci that never exceed the cap end up with exactly the same readPairs (in the same
    // order) as without the cap — only over-cap loci (e.g. centromeric/satellite pileups) get sampled.
    // The replacement choice uses a per-locus RNG seeded from the locusId, and reads are offered in
    // coordinate-stream order (the same order under any --threads), so the retained sample is deterministic
    // and identical across thread counts.
    void offerReadPair(const shared_ptr<FullReadPair>& readPair, std::size_t maxReadPairs)
    {
        if (maxReadPairs == 0 || readPairs.size() < maxReadPairs)
        {
            readPairs.push_back(readPair);
        }
        else
        {
            const uint64_t j = nextRandom(offeredCount_ + 1);  // uniform in [0, offeredCount_]
            if (j < maxReadPairs)
            {
                readPairs[j] = readPair;
            }
        }
        ++offeredCount_;
    }

    // Per-locus reservoir size in read pairs; 0 = unlimited. Derived from --max-depth (see constructor).
    std::size_t maxReadPairs() const { return maxReadPairs_; }

    // True if this locus exceeded its reservoir cap — more read pairs were offered than maxReadPairs_, so the
    // retained set is a reservoir sample (Algorithm R) rather than every read. False when the cap is disabled
    // (maxReadPairs_ == 0) or the locus stayed at/under the cap (its read set is complete, byte-identical to no cap).
    bool reservoirSampled() const { return maxReadPairs_ != 0 && offeredCount_ > maxReadPairs_; }

private:
    // number of read pairs offered so far (before the current one); drives the reservoir replacement odds
    uint64_t offeredCount_ = 0;
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



// the UnpairedReadsCache is used to store reads for which the mate is yet to be found
using UnpairedReadsCache = std::unordered_map<FragmentId, FullRead>;


// Result of genotyping one locus, produced by genotypeLocusFull() (possibly on a worker thread) and
// consumed by the single-threaded writer. Move-only because LocusFindings owns unique_ptrs. The
// LocusAnalyzer is kept alive via shared_ptr so its locusSpec() backs the json/vcf writes without
// copying the heavy graph-bearing LocusSpecification across the thread boundary.
struct LocusOutput {
    enum class Kind { kGenotyped, kFilteredOut, kNoCoverage, kHeuristicOnlySkip, kError };

    unsigned locusIndex = 0;
    Kind kind = Kind::kNoCoverage;
    std::shared_ptr<LocusAnalyzer> analyzer;
    LocusFindings findings;
    std::string locusId;
    std::string message;  // error text when kind == kError
};


// True if a read alignment that spans [alignmentStart, alignmentEnd) on contig `contigId` is fully
// contained within at least one of the given read-extraction regions. Defers to the same containment test
// (isAlignmentContainedInInterval) AnalyzerFinder::query (sample/AnalyzerFinder.cpp) applies when deciding
// which locus a read or its mate belongs to in seeking and streaming modes.
bool readIsContainedInAnyExtractionRegion(
    const std::vector<GenomicRegion>& extractionRegions, int32_t contigId, int64_t alignmentStart,
    int64_t alignmentEnd)
{
    for (const auto& extractionRegion : extractionRegions)
    {
        if (extractionRegion.contigIndex() == contigId
            && isAlignmentContainedInInterval(
                   extractionRegion.start(), extractionRegion.end(), alignmentStart, alignmentEnd))
        {
            return true;
        }
    }
    return false;
}


// Genotype a single locus and return the result instead of writing it inline, so the call can run on a
// worker thread while the json/vcf writes happen later on the single-threaded writer (preserving output
// order). Exceptions are captured into the result (kError) rather than thrown, matching the per-locus
// error isolation of the serial path.
LocusOutput genotypeLocusFull(const ProgramParameters& params, Reference& reference, unsigned locusIndex,
                  const LocusDescription& locusDescription,
                  const std::vector<shared_ptr<FullReadPair>>& readPairs, bool reservoirSampled,
                  graphtools::AlignerSelector& alignerSelector,
                  BamletWriterPtr bamletWriter) {

    LocusOutput out;
    out.locusIndex = locusIndex;
    out.locusId = locusDescription.locusId();

    try {
        // Per-locus full-genotyping timing (thread-CPU clock) starts here so it covers graph decode + read
        // alignment + analyze(); only sampled under --output-genotype-timing.
        const bool recordGenotypingTime = params.outputGenotypeTiming();
        timespec genotypingStart{};
        if (recordGenotypingTime) { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &genotypingStart); }

        const Sex& sampleSex = params.sample().sex();
        // convert the lite-weight LocusDescription object into the heavier LocusSpecification
        LocusSpecification locusSpec = decodeLocusSpecification(locusDescription, reference, params.heuristics(), true);

        // Apply plot policy from CLI flags
        if (params.disableAllPlots())
        {
            locusSpec.setPlotPolicy(PlotPolicy::kNone);
        }
        else if (params.plotAll())
        {
            locusSpec.setPlotPolicy(PlotPolicy::kAll);
        }

        // analyze the read data
        //
        // TODO(off-target IRRs): every read pair is dispatched with RegionType::kTarget here, and
        // the per-locus assignment loop in `doTheAnalysis` only checks reads against
        // locusAndFlanksStart()/End() — never against locusDescription.offtargetRegions().
        // As a result, IRRs that mismap to off-target regions (e.g. C9ORF72, FMR1) are not
        // routed to the RepeatAnalyzer in low-mem-streaming / optimized-streaming. Users
        // needing accurate calls at off-target-dependent loci should currently use seeking
        // or streaming mode. See docs/03_Usage.md "Known limitations" section.
        // Heap-allocate the analyzer and keep it alive in the result so its locusSpec() backs the
        // json/vcf writes performed later (possibly on the writer thread) without copying the heavy,
        // graph-bearing LocusSpecification across the thread boundary. The constructor already
        // std::move's its by-value parameter.
        out.analyzer = std::make_shared<LocusAnalyzer>(std::move(locusSpec), params.heuristics(), bamletWriter,
                                    params.enableAlleleQualityMetrics(), params.enableConsensusSequences());
        // Route each cached read pair to the analyzer the same way seeking and streaming modes do via
        // AnalyzerFinder::query (sample/AnalyzerFinder.cpp). When the two ends are nearby but only ONE of
        // them is contained in the locus's target extraction region, that pair is processed single-ended,
        // so the other end — which lies outside the locus — is never aligned to the repeat graph. Passing
        // such a pair as both reads (the previous behavior here) let the out-of-locus end align spuriously
        // into the repeat node and inflate the genotype, which over-called homopolymers in low-mem and
        // optimized streaming full genotyping even though seeking and streaming called them correctly.
        // Pairs with both ends contained, or with far-apart ends, still pass both reads so that genuine
        // in-repeat reads from a real expansion continue to be counted.
        const std::vector<GenomicRegion>& targetExtractionRegions
            = out.analyzer->locusSpec().targetReadExtractionRegions();
        const std::vector<GenomicRegion>& offtargetExtractionRegions
            = out.analyzer->locusSpec().offtargetReadExtractionRegions();
        for (auto const& readPair : readPairs) {
            // Copy the read + mate locally before processMates: align() mutates the Read in place via
            // reverseComplement(), and the same FullReadPair is shared across multiple loci's caches, so
            // genotyping the shared object directly would be a data race between parallel workers.
            Read read = readPair->firstMate->r;
            const LinearAlignmentStats& readStats = readPair->firstMate->s;
            const int64_t readEnd = readStats.pos + static_cast<int64_t>(read.sequence().length());
            const bool readIsInTargetRegion = readIsContainedInAnyExtractionRegion(
                targetExtractionRegions, readStats.chromId, readStats.pos, readEnd);

            // Single-ended entry (the mate was unmapped/unavailable when the read was cached): analyze the
            // read on its own, the same way seeking mode handles such singletons via processMates(read,
            // nullptr). Without a mate there is no out-of-locus end that could align spuriously into the
            // repeat graph, so the containment-based single/paired routing below does not apply.
            if (!readPair->secondMate) {
                if (readIsInTargetRegion) {
                    out.analyzer->processMates(read, nullptr, locus::RegionType::kTarget, alignerSelector);
                }
                continue;
            }

            Read mate = readPair->secondMate->r;
            const LinearAlignmentStats& mateStats = readPair->secondMate->s;
            const int64_t mateEnd = mateStats.pos + static_cast<int64_t>(mate.sequence().length());

            const bool mateIsInTargetRegion = readIsContainedInAnyExtractionRegion(
                targetExtractionRegions, mateStats.chromId, mateStats.pos, mateEnd);
            const bool readIsInAnyRegion = readIsInTargetRegion
                || readIsContainedInAnyExtractionRegion(
                       offtargetExtractionRegions, readStats.chromId, readStats.pos, readEnd);
            const bool mateIsInAnyRegion = mateIsInTargetRegion
                || readIsContainedInAnyExtractionRegion(
                       offtargetExtractionRegions, mateStats.chromId, mateStats.pos, mateEnd);

            const bool matesAreNearby
                = areMatesNearby(readStats.chromId, readStats.pos, mateStats.chromId, mateStats.pos);

            if (matesAreNearby && readIsInTargetRegion && !mateIsInAnyRegion) {
                out.analyzer->processMates(read, nullptr, locus::RegionType::kTarget, alignerSelector);
            } else if (matesAreNearby && mateIsInTargetRegion && !readIsInAnyRegion) {
                out.analyzer->processMates(mate, nullptr, locus::RegionType::kTarget, alignerSelector);
            } else {
                out.analyzer->processMates(read, &mate, locus::RegionType::kTarget, alignerSelector);
            }
        }

        // genotype the locus
        out.findings = out.analyzer->analyze(sampleSex, boost::none, params.outputPaths().outputPrefix());

        // Thread-CPU genotyping time for this locus (immune to contention from the other --threads workers),
        // measured only under --output-genotype-timing.
        boost::optional<double> genotypingTimeMillis;
        if (recordGenotypingTime) {
            timespec genotypingEnd{};
            clock_gettime(CLOCK_THREAD_CPUTIME_ID, &genotypingEnd);
            genotypingTimeMillis = (genotypingEnd.tv_sec - genotypingStart.tv_sec) * 1e3
                + (genotypingEnd.tv_nsec - genotypingStart.tv_nsec) / 1e6;
        }

        // Annotate each repeat variant's findings: ReservoirSampling when this locus's read set was capped by
        // --max-depth (mirrors the fast path's setReservoirSampling), and GenotypingTimeMillis when timing is on.
        if (reservoirSampled || genotypingTimeMillis) {
            for (auto& variantIdAndFindings : out.findings.findingsForEachVariant) {
                if (auto* repeatFindings = dynamic_cast<RepeatFindings*>(variantIdAndFindings.second.get())) {
                    if (reservoirSampled) { repeatFindings->setReservoirSampling(true); }
                    if (genotypingTimeMillis) { repeatFindings->setGenotypingTimeMillis(*genotypingTimeMillis); }
                }
            }
        }

        // --skip-hom-ref / --skip-missing-genotypes: genotyped but emit no record (still counts as full-genotyped)
        out.kind = shouldFilterLocus(out.analyzer->locusSpec(), out.findings, params.skipHomRef(), params.skipMissingGenotypes())
            ? LocusOutput::Kind::kFilteredOut
            : LocusOutput::Kind::kGenotyped;
    } catch (const std::exception& e) {
        out.kind = LocusOutput::Kind::kError;
        out.message = e.what();
        out.analyzer.reset();
    }

    return out;
}


// ---------------------------------------------------------------------------------------------------------
// Chromosome-stride parallel helpers (shared by the --threads 1 direct-write path and the --threads >1
// per-contig path; see htsLowMemStreamingSampleAnalysis below).
// ---------------------------------------------------------------------------------------------------------

// Read-streamer region for one contig's loci: [minStart - halo, maxEnd + halo] clamped to [0, contigSize].
// halo MUST equal farAwayMateDistanceThreshold (== the locus-activation reach); a narrower halo would
// silently drop reads that activate loci at the contig's locus edges.
GenomicRegion clampedContigRegion(int32_t contigIndex, int64_t minStart, int64_t maxEnd, int64_t contigSize, int halo)
{
    const int64_t begin = std::max<int64_t>(0, minStart - halo);
    const int64_t end = std::min<int64_t>(contigSize, maxEnd + halo);
    return GenomicRegion(contigIndex, begin, end);
}

// Whole-contig regions owned by stride worker `firstContig` under stride `strideT`: the contigs
// c = firstContig, firstContig+strideT, ... < numContigs, each as [0, contigSize), in ascending index order
// (satisfies the region streamer's sorted/non-overlapping precondition). Used by the parallel prepass, whose
// stride scans EVERY contig (cross-contig mates make every contig a potential mate home).
std::vector<GenomicRegion> wholeContigRegionsForStride(const ReferenceContigInfo& contigInfo, int firstContig, int strideT)
{
    std::vector<GenomicRegion> regions;
    const int numContigs = contigInfo.numContigs();
    for (int c = firstContig; c < numContigs; c += strideT)
    {
        regions.emplace_back(static_cast<int32_t>(c), 0, contigInfo.getContigSize(c));
    }
    return regions;
}

// One per loci-bearing contig: the [catalogBegin, catalogEnd) range of the position-sorted catalog holding
// that contig's loci, plus the read-streamer region covering them (with the activation halo). Requires the
// catalog to be position-sorted (--sort-catalog-by position), so each contig's loci form one contiguous run
// and the slices come out in ascending contig-index (== merge) order.
struct ContigGenotypingSlice
{
    int32_t contigIndex;
    std::size_t catalogBegin;
    std::size_t catalogEnd;
    GenomicRegion region;
};

std::vector<ContigGenotypingSlice> perContigGenotypingSlices(
    const LocusDescriptionCatalog& catalog, const ReferenceContigInfo& contigInfo, int halo)
{
    std::vector<ContigGenotypingSlice> slices;
    std::size_t k = 0;
    while (k < catalog.size())
    {
        const int32_t contigIndex = catalog[k].locusContigIndex();
        int64_t minStart = catalog[k].locusAndFlanksStart();
        int64_t maxEnd = catalog[k].locusAndFlanksEnd();
        std::size_t m = k;
        while (m < catalog.size() && catalog[m].locusContigIndex() == contigIndex)
        {
            minStart = std::min(minStart, catalog[m].locusAndFlanksStart());
            maxEnd = std::max(maxEnd, catalog[m].locusAndFlanksEnd());
            ++m;
        }
        slices.push_back({ contigIndex, k, m,
            clampedContigRegion(contigIndex, minStart, maxEnd, contigInfo.getContigSize(contigIndex), halo) });
        k = m;
    }
    return slices;
}

// Genome-wide typical read length: the max sequence length over the first kProbeReadCount primary aligned
// reads, computed ONCE at startup and held constant for the whole run (and across thread counts). This makes
// the prepass cache contents and the genotyping admission predicate thread-count-independent. It replaces the
// previous per-sweep first-read latch that grew mid-traversal and was therefore read-order-dependent.
int probeTypicalReadLength(const InputPaths& inputPaths)
{
    const int kProbeReadCount = 1000;
    const int kFallbackReadLength = 150;
    htshelpers::HtsFileStreamer probeStreamer(inputPaths.htsFile(), inputPaths.reference(), 1u);
    int maxReadLength = 0;
    int readCount = 0;
    while (readCount < kProbeReadCount && probeStreamer.tryReadingNextPrimaryAlignment()
           && probeStreamer.isStreamingAlignedReads())
    {
        maxReadLength = std::max(maxReadLength, static_cast<int>(probeStreamer.decodeRead().sequence().length()));
        ++readCount;
    }
    return maxReadLength > 0 ? maxReadLength : kFallbackReadLength;
}


// Scans the read file for read pairs whose mates are far away (different contig, or >= the threshold apart)
// and that are near a catalog locus, caching them in `mateExtractor`'s localCache_. With an empty
// `streamerRegions` it streams the WHOLE file (used at --threads 1). With a non-empty `streamerRegions` it
// index-seeks only those regions (used by each chromosome-stride prepass worker over its owned whole
// contigs), so several workers can build disjoint private shards concurrently. The caller logs the cache
// size (per-call counts would interleave across workers).
void prepareCache(
    const ProgramParameters& params, Reference& reference, const GenomeQueryCollection& genomeQuery,
    htshelpers::MateExtractor& mateExtractor, const int farAwayMateDistanceThreshold, const int typicalReadLength,
    const std::vector<GenomicRegion>& streamerRegions)
{
    const InputPaths& inputPaths = params.inputPaths();

    // keep track of the previous chromosome id to make sure that the read input file has the same chromosome ordering
    // as the reference fasta
    int previousReadContigId = -1;

    // Whole-file mode (streamerRegions empty) uses all decompression threads and validates the read file's
    // contig order against the reference FASTA globally. In chromosome-stride mode (non-empty regions) each
    // worker index-seeks only its owned whole contigs with a single decompression thread (to avoid T-way bgzf
    // oversubscription across concurrent workers); there the contig-order check only asserts the worker's own
    // stride is ascending — the global FASTA-vs-BAM order is not revalidated in the parallel prepass.
    const unsigned streamerDecompressionThreads = streamerRegions.empty() ? std::min(params.threadCount, 12) : 1u;
    std::unique_ptr<htshelpers::HtsFileStreamer> readStreamer = streamerRegions.empty()
        ? std::make_unique<htshelpers::HtsFileStreamer>(
              inputPaths.htsFile(), inputPaths.reference(), streamerDecompressionThreads)
        : std::make_unique<htshelpers::HtsFileStreamer>(
              inputPaths.htsFile(), inputPaths.reference(), inputPaths.htsIndexFile(), streamerRegions,
              streamerDecompressionThreads);

    while (readStreamer->tryReadingNextPrimaryAlignment() && readStreamer->isStreamingAlignedReads())
    {
        //skip unpaired reads
        if (!readStreamer->currentIsPaired()) {
            continue;
        }

        const int32_t readContigId = readStreamer->currentReadContigId();
        const int64_t readStart = readStreamer->currentReadPosition();
        int64_t readEnd = readStart + typicalReadLength;

        // handle readStreamer transitioning from one chromosome to the next
        if (previousReadContigId != readContigId) {
            const auto& chromosomeName = reference.contigInfo().getContigName(readContigId);

            // make sure chromosome ordering is consistent between the read input file and the reference fasta
            if (previousReadContigId != -1 && previousReadContigId > readContigId) {
                const auto& previousChromosomeName = reference.contigInfo().getContigName(previousReadContigId);
                throw std::runtime_error(
                    "The input read file contains reads from " + chromosomeName + " before " + previousChromosomeName
                    + ", which is a different chromosome order than the reference input FASTA file.");
            }
            previousReadContigId = readContigId;
        }

        const int32_t mateContigId = readStreamer->currentMateContigId();
        const int64_t mateStart = readStreamer->currentMatePosition();
        int64_t mateEnd = mateStart + typicalReadLength;

        // since prepareCache only caches reads with far-away mates, skip reads with nearby mates
        const bool isMateFarAway = readContigId != mateContigId || std::abs(readStart - mateStart) >= farAwayMateDistanceThreshold;
        if (!isMateFarAway) {
            continue;
        }

        //if read & mate don't overlap any locus, skip this read pair
        const bool isReadNearTargetRegion = genomeQuery.targetRegionMask.query(readContigId, readStart);
        const bool isMateNearTargetRegion = genomeQuery.targetRegionMask.query(mateContigId, mateStart);
        if (!isReadNearTargetRegion && !isMateNearTargetRegion) {
            continue;
        }

        if (genomeQuery.analyzerFinder.query(readContigId, readStart, readEnd, mateContigId, mateStart, mateEnd).empty()) {
            continue;
        }

        // the read and/or the the mate are near a locus and are not close to each other, so cache the read in the mateExtractor

        FullRead fullRead = decodeRead(*readStreamer);
        // Copy the ReadId out before the std::move below: passing both `fullRead.r.readId()` and
        // `std::move(fullRead)` as function args is use-after-move in C++17 (parameter init order is
        // indeterminate, and unlike std::pair's mem-init-list there is no ordering guarantee here).
        ReadId mateReadId = fullRead.r.readId();
        mateExtractor.addMateToCache(mateReadId, std::move(fullRead));
    }
}


// Running tally of how each locus was genotyped, summed across all doTheAnalysis calls in a run so the
// "Processed N loci ..." summary is emitted ONCE for the whole run (at --threads >1 there is one call per
// loci-bearing contig, so a per-call summary would print many partial, interleaved lines).
struct GenotypingCounts {
    int fast = 0;
    int full = 0;
    int skipped = 0;
    int zeroCoverage = 0;
};

void logGenotypingSummary(const ProgramParameters& params, const GenotypingCounts& counts)
{
    const int totalProcessed = counts.fast + counts.full + counts.skipped + counts.zeroCoverage;
    if (totalProcessed > 0) {
        const float fastPct = 100.0f * counts.fast / totalProcessed;
        const float fullPct = 100.0f * counts.full / totalProcessed;
        if (counts.skipped > 0) {
            const float skippedPct = 100.0f * counts.skipped / totalProcessed;
            spdlog::info("Processed {} ({:.1f}%) loci via fast genotyping, {} ({:.1f}%) loci via full genotyping, and skipped {} ({:.1f}%) loci",
                add_commas_at_thousands(counts.fast), fastPct,
                add_commas_at_thousands(counts.full), fullPct,
                add_commas_at_thousands(counts.skipped), skippedPct);
        } else {
            spdlog::info("Processed {} ({:.1f}%) loci via fast genotyping, and {} ({:.1f}%) loci via full genotyping",
                add_commas_at_thousands(counts.fast), fastPct,
                add_commas_at_thousands(counts.full), fullPct);
        }
        if (counts.zeroCoverage > 0) {
            const float zeroCoveragePct = 100.0f * counts.zeroCoverage / totalProcessed;
            if (params.skipMissingGenotypes()) {
                spdlog::info("Excluded {} ({:.1f}%) loci with zero coverage from output (--skip-missing-genotypes)",
                    add_commas_at_thousands(counts.zeroCoverage), zeroCoveragePct);
            } else {
                spdlog::info("Emitted {} ({:.1f}%) no-call records for loci with zero coverage",
                    add_commas_at_thousands(counts.zeroCoverage), zeroCoveragePct);
            }
        }
    }
}

// Accumulates this call's per-locus genotyping tallies into `counts` (summed across calls by the caller).
// The "Processed N loci ..." summary is logged here only in whole-file mode (showProgress, i.e. --threads 1);
// at --threads >1 the per-contig calls stay silent and the orchestrator logs one aggregate summary.
void doTheAnalysis(
    const ProgramParameters& params, Reference& reference, LocusDescriptionCatalog& locusDescriptionCatalog,
    const GenomeQueryCollection& genomeQuery, htshelpers::MateExtractor& mateExtractor,
    BamletWriterPtr bamletWriter, const int farAwayMateDistanceThreshold, const int typicalReadLength,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter,
    const std::vector<GenomicRegion>& streamerRegions, GenotypingCounts& counts)
{
    if (locusDescriptionCatalog.empty()) {
        return;
    }

    // The per-locus reservoir-sampling cap is derived from --max-depth and the locus window width, and is
    // computed and stored inside each LocusCache at construction (see LocusCache); the offer/admit calls below
    // read it back via locusCache->maxReadPairs(). 0 (from --max-depth 0) disables the cap.

    int fastGenotypedCount = 0;
    int fullGenotypedCount = 0;
    int skippedCount = 0;
    int zeroCoverageCount = 0;

    // Constructed below only in whole-file mode (showProgress). With a non-empty streamerRegions multiple
    // per-contig chromosome-stride workers run this function concurrently and a progress_display writes its
    // layout to std::cout on construction, so one per worker would garble the terminal. Left null when
    // showProgress is false.
    std::unique_ptr<boost::timer::progress_display> progressBar;

    const int threadCount = params.threadCount;
    const unsigned htsDecompressionThreads(std::min(threadCount, 12));

    const InputPaths& inputPaths = params.inputPaths();

    // jsonWriter and vcfWriter are injected by the caller. Each sweep finalizes its own writer via the
    // jsonWriter.close()/vcfWriter.close() calls at the end of this function.

    // When streamerRegions is non-empty the read stream is restricted to those regions (one per-contig
    // chromosome-stride genotyping worker); gating the progress bar on whole-file mode avoids garbled
    // interleaved output across concurrent sweeps.
    const bool showProgress = streamerRegions.empty();
    if (showProgress) {
        progressBar = std::make_unique<boost::timer::progress_display>(locusDescriptionCatalog.size());
    }

    // nextLocusDescriptionIndex points to the first locusDescription in the locusDescriptionCatalog whose
    // reference regions + flanks are entirely to the right of the current read
    unsigned nextLocusDescriptionIndex = 0;
    LocusDescription* nextLocusDescription = &locusDescriptionCatalog[nextLocusDescriptionIndex];

    // list of active locus indices. An active locus L is any locus from the catalog for which the current read overlaps:
    // (L.locusAndFlanksStart - farAwayMateDistanceThreshold, L.locusAndFlanksEnd + farAwayMateDistanceThreshold)
    std::list<unsigned> locusIndicesCollectingReads;

    // list of locus indices that are entirely to the left of the current read position, meaning that all relevant reads
    // have been added to their locusCaches, and so they are ready to be analyzed. As reads are parsed, locusIndices
    // move from the locusIndicesCollectingReads list to this list.
    std::list<unsigned> locusIndicesReadyForAnalysis;

    // map of locusIndex to LocusCache object which stores all the read pairs needed to analyze a locus
    unordered_map<int, shared_ptr<LocusCache> > locusCachesMap;

    // unpaired read cache
    UnpairedReadsCache unpairedReadsCache;

    // keep track of the previous chromosome id to make sure that the read input file has the same chromosome ordering
    // as the reference fasta
    int previousReadContigId = -1;

    // A single AlignerSelector (and thus a single PinnedDagAligner with its large DP score matrices)
    // reused across all full-genotyped loci on the serial path. The aligner is graph-agnostic — the graph
    // is supplied per align() call via the seed path — so reuse is output-neutral and keeps the matrix
    // buffers warm instead of reallocating ~hundreds of KB per locus.
    graphtools::AlignerSelector alignerSelector(params.heuristics().alignerType());

    // Genotyping runs serially within this call (one contig / whole-file sweep per call); cross-contig
    // parallelism is handled by the caller. genotypeLocusFull copies each Read before processMates instead
    // of passing the shared FullReadPair members by reference. Since align() mutates a Read in place
    // (reverseComplement) and one read pair can belong to multiple loci, this fixes a prior cross-locus
    // mutation of shared reads and so may change output for catalogs with overlapping loci relative to
    // earlier releases.

    // Write one finished locus result, inline as each locus is genotyped, in catalog order. Counter
    // semantics: a locus reaching full genotyping is counted as full-genotyped even if it produces no record
    // because of --skip-hom-ref / --skip-missing-genotypes or an error.
    auto writeOutput = [&](LocusOutput& out) {
        switch (out.kind) {
            case LocusOutput::Kind::kNoCoverage:
                zeroCoverageCount++;
                writeZeroCoverageRecord(params, reference, locusDescriptionCatalog[out.locusIndex], jsonWriter, vcfWriter);
                break;
            case LocusOutput::Kind::kHeuristicOnlySkip:
                skippedCount++;
                jsonWriter.addSkippedRecord(out.locusId, "heuristic_only_mode");
                break;
            case LocusOutput::Kind::kError:
                fullGenotypedCount++;
                spdlog::error("Error while processing {}: {}", out.locusId, out.message);
                break;
            case LocusOutput::Kind::kFilteredOut:
                fullGenotypedCount++;
                break;
            case LocusOutput::Kind::kGenotyped:
                fullGenotypedCount++;
                jsonWriter.addRecord(out.analyzer->locusSpec(), out.findings);
                vcfWriter.addRecords(out.analyzer->locusSpec(), out.findings);
                break;
        }
    };

    // Analyze every locus index in locusIndicesReadyForAnalysis, then clear the queue.
    // Defined as a lambda so the same code path runs both during streaming and during the post-EOF drain
    // (otherwise loci still in flight when EOF is reached would silently disappear from the output).
    auto processReadyLoci = [&]() {
        if (locusIndicesReadyForAnalysis.empty()) {
            return;
        }
        for (auto const& locusIndex : locusIndicesReadyForAnalysis) {
            if (locusCachesMap.find(locusIndex) == locusCachesMap.end()) {
                // Locus has zero coverage; emit a no-call record identical to seeking/streaming mode (see
                // writeZeroCoverageRecord) so the output is consistent across analysis modes.
                zeroCoverageCount++;
                writeZeroCoverageRecord(params, reference, locusDescriptionCatalog[locusIndex], jsonWriter, vcfWriter);
                continue;
            }
            const shared_ptr<LocusCache> locusCache = locusCachesMap[locusIndex];

            // TODO check if # of reads >= minLocusCoverage and skip locus if not
            bool needToProcessSlowly = true;
            if (params.analysisMode() == AnalysisMode::kOptimizedStreaming) {
                const bool doneGenotyping = processLocusFast(params, reference,
                    locusDescriptionCatalog[locusIndex], locusCache->readPairs, locusCache->reservoirSampled(),
                    jsonWriter, vcfWriter);

                needToProcessSlowly = !doneGenotyping;
                if (doneGenotyping) {
                    fastGenotypedCount++;
                }
            }

            if (needToProcessSlowly) {
                if (params.heuristicGenotypingOnly()) {
                    skippedCount++;
                    jsonWriter.addSkippedRecord(locusDescriptionCatalog[locusIndex].locusId(), "heuristic_only_mode");
                } else {
                    LocusOutput out = genotypeLocusFull(params, reference, locusIndex,
                        locusDescriptionCatalog[locusIndex], locusCache->readPairs, locusCache->reservoirSampled(),
                        alignerSelector, bamletWriter);
                    writeOutput(out);
                }
            }

            try {
                locusCachesMap.erase(locusIndex);
            } catch (const std::exception& e) {
                spdlog::error("Unexpected error while erasing locus {}: {}. Skipping this locus.", locusDescriptionCatalog[locusIndex].locusId(), e.what());
            }
        }
        locusIndicesReadyForAnalysis.clear();
    };

    if (showProgress) {
        spdlog::info("Starting to process {} loci", add_commas_at_thousands(locusDescriptionCatalog.size()));
    }
    // In whole-file mode use all decompression threads; in region-restricted mode use a single one to avoid
    // worker-count x htsDecompressionThreads bgzf thread oversubscription across concurrent stride workers.
    const unsigned streamerDecompressionThreads = streamerRegions.empty() ? htsDecompressionThreads : 1u;
    std::unique_ptr<htshelpers::HtsFileStreamer> readStreamer = streamerRegions.empty()
        ? std::make_unique<htshelpers::HtsFileStreamer>(
              inputPaths.htsFile(), inputPaths.reference(), streamerDecompressionThreads)
        : std::make_unique<htshelpers::HtsFileStreamer>(
              inputPaths.htsFile(), inputPaths.reference(), inputPaths.htsIndexFile(), streamerRegions,
              streamerDecompressionThreads);
    while (readStreamer->tryReadingNextPrimaryAlignment() && readStreamer->isStreamingAlignedReads())
    {

        //skip unpaired reads
        if (!readStreamer->currentIsPaired()) {
            continue;
        }

        const int32_t readContigId = readStreamer->currentReadContigId();
        const int64_t readStart = readStreamer->currentReadPosition();
        int64_t readEnd = readStart + typicalReadLength;

        // When the mate is unmapped, htslib reports mtid=-1/mpos=-1; fall back to the read's own
        // coordinates (as seeking mode does) so the mate is treated as nearby and we never attempt a
        // region jump to an invalid contig, which would throw in extractMate.
        const bool isMateMapped = readStreamer->currentMateContigId() >= 0;
        const int32_t mateContigId = isMateMapped ? readStreamer->currentMateContigId() : readContigId;
        const int64_t mateStart = isMateMapped ? readStreamer->currentMatePosition() : readStart;
        int64_t mateEnd = mateStart + typicalReadLength;

        const bool isMateFarAway = readContigId != mateContigId || std::abs(readStart - mateStart) >= farAwayMateDistanceThreshold;

        // handle readStreamer transitioning from one chromosome to the next
        if (previousReadContigId != readContigId) {
            const auto& chromosomeName = reference.contigInfo().getContigName(readContigId);

            // make sure chromosome ordering is consistent between the read input file and the reference fasta
            if (previousReadContigId != -1 && previousReadContigId > readContigId) {
                const auto& previousChromosomeName = reference.contigInfo().getContigName(previousReadContigId);
                throw std::runtime_error(
                    "The input read file contains reads from " + chromosomeName + " before " + previousChromosomeName
                    + ", which is a different chromosome order than the reference input FASTA file.");
            }
            previousReadContigId = readContigId;

            // load the next chromosome into the cache
            if (showProgress) {
                spdlog::info("Processing loci on {}", chromosomeName);
            }
            reference.clearContigCache();
            try {
                reference.loadContigIntoCache(chromosomeName);
            } catch (const MissingContigError&) {
                // The read file's header contains a contig absent from the reference FASTA (e.g. _fix
                // patch contigs in assembly38 BAMs aligned against an hg38.fa without patches). No catalog
                // locus can target such a contig, so warn and skip its reads rather than aborting the run.
                // The per-read target/analyzer checks below already discard reads on no-locus contigs, so
                // leaving the contig cache empty is safe (no getSequence call is reached).
                spdlog::warn("Skipping reads on {}: contig not present in the reference FASTA", chromosomeName);
                continue;
            }
        }

        //if read & mate don't overlap any locus, skip this read pair
        const bool isReadNearTargetRegion = genomeQuery.targetRegionMask.query(readContigId, readStart);
        const bool isMateNearTargetRegion = genomeQuery.targetRegionMask.query(mateContigId, mateStart);
        if (!isReadNearTargetRegion && !isMateNearTargetRegion) {
            continue;
        }

        if (genomeQuery.analyzerFinder.query(readContigId, readStart, readEnd, mateContigId, mateStart, mateEnd).empty()) {
            continue;
        }

        // check if some windows in the active locus list are now to the left of the current read. If yes, move
        // them out of the active list since read collection is done, and add them to the list of loci to analyze

        //bool activeLociRemoved = false;
        auto i = locusIndicesCollectingReads.begin();
        while (i != locusIndicesCollectingReads.end())
        {
            const LocusDescription* locusDescription = &locusDescriptionCatalog[*i];
            if (isAToTheLeftOfB(
                locusDescription->locusContigIndex(), locusDescription->locusAndFlanksEnd() + farAwayMateDistanceThreshold,
                readContigId, readStart))
            {
                // transfer this locusIndex from locusIndicesCollectingReads to locusIndicesReadyForAnalysis
                locusIndicesReadyForAnalysis.push_back(*i);
                locusIndicesCollectingReads.erase(i++);
            } else {
                i++;
            }
        }


        //bool activeLociAdded = false;
        // Increment nextLocusDescriptionIndex until nextLocusDescription window is to the right of the current read.
        // Append locusIndices to the locusIndicesCollectingReads list if their window overlaps the current read start/end.
        while (nextLocusDescriptionIndex < locusDescriptionCatalog.size() && isAToTheLeftOfB(
               nextLocusDescription->locusContigIndex(), nextLocusDescription->locusAndFlanksStart(), readContigId, readEnd, true))
        {
            // if the current read is contained within the nextLocusDescription window or overlaps its right edge, add
            // the locusIndex to the locusIndicesCollectingReads list
            if (not isAToTheLeftOfB(nextLocusDescription->locusContigIndex(), nextLocusDescription->locusAndFlanksEnd(), readContigId, readStart)) {
                //std::cout << "Adding locus " << nextLocusDescriptionIndex << " (" << nextLocusDescription.locusId << ") to locusIndicesCollectingReads" << std::endl;
                locusIndicesCollectingReads.push_back(nextLocusDescriptionIndex);

                //activeLociAdded = true;
            } else {
                // this locus has no coverage so skip the locusIndicesCollectingReads list and add it directly to the list of loci to analyze.
                //spdlog::warn("Locus {} has zero coverage", nextLocusDescription.locusId);

                /*
                spdlog::warn("Next locus {}:{}-{} description {} overlaps end of read {}:{}-{}, but it's start is to the left of the read start",
                                             nextLocusDescription.locusContigIndex+1, nextLocusDescription.locusAndFlanksStart, nextLocusDescription.locusAndFlanksEnd,
                                             nextLocusDescription.locusId, readContigId + 1, readStart, readEnd,
                                );
                */
                locusIndicesReadyForAnalysis.push_back(nextLocusDescriptionIndex);
            }
            //std::cout << "Incrementing nextLocusDescriptionIndex from " << nextLocusDescriptionIndex << " to " << nextLocusDescriptionIndex + 1 << std::endl;
            ++nextLocusDescriptionIndex;
            if (showProgress) { ++(*progressBar); }
            if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size()) {
                break;
            }
            nextLocusDescription = &locusDescriptionCatalog[nextLocusDescriptionIndex];
        }

        // For a far-away pair whose ONLY in-locus end is the upstream one, that upstream end (streamed
        // first in coordinate order) already prefetched its mate via mateExtractor and pushed the pair
        // into every active locus cache. The downstream end then arrives later and would prefetch and
        // re-push the same pair, so we skip it. We must NOT skip when the downstream end is itself near
        // a target region: in that case its locus only becomes active once the stream reaches the
        // downstream position, so the upstream end could not have pushed the pair there, and skipping
        // would silently drop the pair for that locus. (Pairs whose two ends are both in one locus are
        // de-duplicated below via LocusCache::seenFragmentIds.)
        if (isMateFarAway && !isAToTheLeftOfB(readContigId, readStart, mateContigId, mateStart, false)
            && !isReadNearTargetRegion) {
            continue;
        }

        // fully parse the current read
        FullRead fullRead = decodeRead(*readStreamer);

        // A read whose mate is unmapped (htslib reports mtid<0) has no mate record to pair with: the
        // unmapped mate either never appears in the coordinate-sorted stream or arrives with contig -1,
        // which ends streaming (HtsFileStreamer::isStreamingAlignedReads). Parking the read in
        // unpairedReadsCache to await a pair that never completes would silently drop it (only completed
        // pairs are added to the locus caches below). Instead add it now as single-ended evidence to every
        // active locus whose window contains it, mirroring seeking mode's processMates(read, nullptr)
        // handling of such singletons.
        if (!isMateMapped) {
            shared_ptr<FullReadPair> singleEndedRead = nullptr;
            for (auto const& locusIndex : locusIndicesCollectingReads) {
                const LocusDescription* locusDescription = &locusDescriptionCatalog[locusIndex];
                if (!doesAcontainB(
                        locusDescription->locusContigIndex(), locusDescription->locusAndFlanksStart(),
                        locusDescription->locusAndFlanksEnd(), readContigId, readStart, readEnd)) {
                    continue;
                }
                if (locusCachesMap.find(locusIndex) == locusCachesMap.end()) {
                    locusCachesMap[locusIndex] = std::make_shared<LocusCache>(
                        locusIndex,
                        GenomicRegion {
                            locusDescription->locusContigIndex(),
                            locusDescription->locusAndFlanksStart(),
                            locusDescription->locusAndFlanksEnd()
                        },
                        locusDescription->locusId(), params.maxDepth(), typicalReadLength);
                }
                shared_ptr<LocusCache> locusCache = locusCachesMap[locusIndex];
                // de-duplicate in case the same fragment is encountered more than once for this locus
                if (!locusCache->admitFragment(fullRead.r.fragmentId(), locusCache->maxReadPairs())) {
                    continue;
                }
                if (singleEndedRead == nullptr) {
                    singleEndedRead = std::make_shared<FullReadPair>();
                    singleEndedRead->firstMate = fullRead;
                }
                locusCache->offerReadPair(singleEndedRead, locusCache->maxReadPairs());
            }
            continue;
        }

        // if mate if far away, retrieve it using the MateExtractor and add it to the unpairedReadsCache so that
        // the read pair is available immediately. This way we don't have to wait until the read position reaches
        // the mate locus, allowing the locus to be processed earlier.
        if (isMateFarAway && unpairedReadsCache.find(fullRead.r.fragmentId()) == unpairedReadsCache.end()) {

            //retrieve far-away mate and add it to the unpairedReadsCache
            //auto mateGenomicRegion = GenomicRegion(mateContigId, mateStart, mateStart + 1);
            auto mateNumber = fullRead.r.readId().mateNumber() == MateNumber::kFirstMate ? MateNumber::kSecondMate : MateNumber::kFirstMate;
            auto mateReadId = ReadId(fullRead.r.fragmentId(), mateNumber);
            // NOTE: For smaller catalogs, the "prepare cache" step is skipped, and far-away mates are retrieved on the fly
            // using the mateExtractor. In that case, when the mateExtractor retrieves a far-away mate, it's useful to fetch
            // reads not just at the exact genomic position of the mate, but also in a window around it, whose size is set by
            // mateExtractorWindowSize. The reason is that far-away mates from different repetitive loci often mismap to the
            // same regions of the genome. This window often captures multiple far-away mates that will be needed later. The
            // same reasoning underlies the --cache-mates option which improves performance (ie. speed) in "seeking" mode.
            const int mateExtractorWindowSize = 500;
            auto mateGenomicRegion = GenomicRegion(mateContigId, mateStart - mateExtractorWindowSize, mateStart + mateExtractorWindowSize + 1);
            boost::optional<FullRead> extractedMate = mateExtractor.extractMate(mateReadId, mateGenomicRegion);
            if (extractedMate == boost::none) {
                spdlog::warn("No mate found for read {}", fullRead.r.fragmentId());
                continue;
            }
            unpairedReadsCache.emplace(fullRead.r.fragmentId(), std::move(*extractedMate));
        }

        // get the mate of the current read from the unpairedReadsCache
        // TODO make the unpaired read cache keep both the read and the mate for far away mates? or skip reads where mate overlaps locus but read doesn't?
        auto unpairedReadIter = unpairedReadsCache.find(fullRead.r.fragmentId());
        if (unpairedReadIter == unpairedReadsCache.end())
        {
            // confirm that the mate is to the right of the read and is < farAwayMateDistanceThreshold away
            if (isMateFarAway) {
                spdlog::error("Mate not found in unpairedReadsCache even though it is far away from the read");
            } else if (not isAToTheLeftOfB(readContigId, readStart, mateContigId, mateStart, true)) {
                spdlog::error("Mate not found in unpairedReadsCache even though it is to the left of the read");
            }
            // here we know the mate is to the right of the read and is < farAwayMateDistanceThreshold away, so we can
            // wait for the read position to reach it and all active locusDescriptions that contain this read will
            // still be in the locusIndicesCollectingReads list so this read and its mate will be added to their LocusCaches
            // when the mate is reached
            unpairedReadsCache.emplace(fullRead.r.fragmentId(), std::move(fullRead));
            continue;
        }

        FullRead fullMate = std::move(unpairedReadIter->second);
        unpairedReadsCache.erase(unpairedReadIter);

        // Both the read & mate are now available, so recompute readEnd/mateEnd from their ACTUAL sequence
        // lengths (the readEnd/mateEnd above used the genome-wide typicalReadLength estimate). typicalReadLength
        // is a constant seeded once at startup (max over the first primary reads) and is deliberately NOT grown
        // here, so the cache-admission and locus-activation decisions are thread-count-independent.
        const int readSequenceLength = (int) fullRead.r.sequence().length();
        if (readSequenceLength != typicalReadLength) {
            readEnd = readStart + readSequenceLength;
        }
        const int mateSequenceLength = (int) fullMate.r.sequence().length();
        if (mateSequenceLength != typicalReadLength) {
            mateEnd = mateStart + mateSequenceLength;
        }

        // iterate over locusIndicesCollectingReads and, for each locus, add the current read pair to its LocusCache
        shared_ptr<FullReadPair> readPair = nullptr;
        for (auto const& locusIndex : locusIndicesCollectingReads) {
            const LocusDescription* locusDescription = &locusDescriptionCatalog[locusIndex];

            // make sure this locusDescription window fully contains the read or its mate
            if ( !doesAcontainB(
                        locusDescription->locusContigIndex(), locusDescription->locusAndFlanksStart(), locusDescription->locusAndFlanksEnd(),
                        readContigId, readStart, readEnd)
                && !doesAcontainB(
                        locusDescription->locusContigIndex(), locusDescription->locusAndFlanksStart(), locusDescription->locusAndFlanksEnd(),
                        mateContigId, mateStart, mateEnd) ) {
                continue;
            }

            //if the LocusCache object doesn't exist yet for the current locusIndex, create it
            shared_ptr<LocusCache> locusCache;
            if (locusCachesMap.find(locusIndex) == locusCachesMap.end()) {
                //create a new LocusCache object add it to the locusCachesMap and locusCaches list
                locusCache = std::make_shared<LocusCache>(locusIndex,
                    GenomicRegion {
                        locusDescription->locusContigIndex(),
                        locusDescription->locusAndFlanksStart(),
                        locusDescription->locusAndFlanksEnd()
                    },
                    locusDescription->locusId(), params.maxDepth(), typicalReadLength);
                locusCachesMap[locusIndex] = locusCache;
            }

            locusCache = locusCachesMap[locusIndex];

            // skip if this fragment's pair was already added to this locus (e.g. a far-away pair whose
            // upstream and downstream ends are both contained in this locus window are streamed, and
            // processed, in two separate iterations)
            if (!locusCache->admitFragment(fullRead.r.fragmentId(), locusCache->maxReadPairs())) {
                continue;
            }

            if (readPair == nullptr) {
                readPair = std::make_shared<FullReadPair>(fullRead, fullMate);
            }

            locusCache->offerReadPair(readPair, locusCache->maxReadPairs());
        }

		// if all loci have been processed, break out of the loop
		if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size() && locusIndicesCollectingReads.empty() && locusIndicesReadyForAnalysis.empty()) {
            // done processing the locusDescriptionCatalog
            break;
        }

        // TODO clean up unpairedReadsCache, locusCaches, and mateExtractorCache?
        // every 500 bases, process and then destroy caches that are kExtensionLength to the left of the current read

        processReadyLoci();

        if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size() && locusIndicesCollectingReads.empty() && locusIndicesReadyForAnalysis.empty()) {
            // done processing the locusDescriptionCatalog
            break;
        }
    }

    // EOF reached. Drain any loci that were still in flight so they appear in the output rather than
    // being silently dropped — both loci still collecting reads, and catalog loci the read stream never reached.
    for (auto locusIndex : locusIndicesCollectingReads) {
        locusIndicesReadyForAnalysis.push_back(locusIndex);
    }
    locusIndicesCollectingReads.clear();
    while (nextLocusDescriptionIndex < locusDescriptionCatalog.size()) {
        locusIndicesReadyForAnalysis.push_back(nextLocusDescriptionIndex);
        ++nextLocusDescriptionIndex;
        if (showProgress) { ++(*progressBar); }
    }
    processReadyLoci();

    reference.clearContigCache();

    if (showProgress) {
        spdlog::info("Writing JSON and VCF output files");
    }

    jsonWriter.close();
    vcfWriter.close();

    GenotypingCounts callCounts;
    callCounts.fast = fastGenotypedCount;
    callCounts.full = fullGenotypedCount;
    callCounts.skipped = skippedCount;
    callCounts.zeroCoverage = zeroCoverageCount;
    // The genotyping summary is logged once by the orchestrator (htsLowMemStreamingSampleAnalysis) after this
    // returns -- in whole-file mode from the counts filled here, in --threads >1 mode from the summed per-worker
    // tallies -- so it always prints exactly once at the end, independent of showProgress (the per-contig progress
    // gate). This call only accumulates its tallies into the caller's counts.
    counts.fast += callCounts.fast;
    counts.full += callCounts.full;
    counts.skipped += callCounts.skipped;
    counts.zeroCoverage += callCounts.zeroCoverage;
}



void htsLowMemStreamingSampleAnalysis(
    LocusDescriptionCatalog& locusDescriptionCatalog,
    const ProgramParameters& programParams,
    Reference& reference,
    BamletWriterPtr bamletWriter)
{
    //initialize the genomeQuery object
    GenomeQueryCollection genomeQuery(locusDescriptionCatalog);

    // initialize the mateExtractor for caching far-away mates
    //
    // NOTE: For large catalogs, two-pass analysis mode does a two-pass traversal of the bam/cram file.
    // On the first pass, it reads the file from start to finnish, looking for read pairs where mates are
    // aligned far away from each other. It caches these mates in memory within the mateExtractor.
    // On the second pass, it performs the core ExpansionHunter analysis without having to do any seek i/o operations
    // since all distant mates are already in the cache.
    // For small catalogs, it skips the first pass and just does a single-pass analysis just like the original
    // ExpansionHunter streaming mode, retrieving distant mates on the fly via seek i/o operations using mateExtractor.
    // There is a trade-off where increasing the farAwayMateDistanceThreshold will reduce memory requirements of the
    // "prepare cache" step for large catalogs, and will reduce i/o seeks for small catalogs during the main traversal.
    // The downside is that it may increase memory requirements for both large and small catalogs during the main
    // traversal because it will lead to a larger active region which extends from the start of the flank of the
    // nearest-to-the-left locusDescription (eg. L) minus farAwayMateDistanceThreshold to the the start of the current
    // read, (eg. (current read position - (L.locusAndFlanksStart - farAwayMateDistanceThreshold))

    const int farAwayMateDistanceThreshold = 1000;  // base pairs
    const InputPaths& inputPaths = programParams.inputPaths();
    const int threadCount = std::max(1, programParams.threadCount);

    // Phase 0: seed a single genome-wide typicalReadLength and build the far-away-mate admission predicate
    // from the FULL catalog (so the prepass predicate is a superset of every per-contig genotyping predicate).
    const int typicalReadLength = probeTypicalReadLength(inputPaths);

    // --threads 1: serial whole-file far-away-mate prepass, then the original whole-file single-pass analysis
    // writing straight to the final JSON/VCF — no temp files, no merge, no region restriction, no cache freeze.
    // This IS the golden serial baseline that every higher thread count must reproduce byte-for-byte.
    if (threadCount == 1)
    {
        htshelpers::MateExtractor mateExtractor(
            inputPaths.htsFile(), inputPaths.htsIndexFile(), inputPaths.reference(), true, farAwayMateDistanceThreshold);
        spdlog::info("Caching read pairs with distant mates");
        prepareCache(programParams, reference, genomeQuery, mateExtractor, farAwayMateDistanceThreshold,
                     typicalReadLength, {});
        spdlog::info("Added {} reads to the mate cache", add_commas_at_thousands(mateExtractor.mateCacheSize()));

        IterativeJsonWriter jsonWriter(programParams.sample(), reference.contigInfo(), programParams.outputPaths().json(), programParams.copyCatalogFields());
        IterativeVcfWriter vcfWriter(programParams.sample().id(), reference, programParams.outputPaths().vcf());
        GenotypingCounts counts;  // filled by doTheAnalysis; the summary is logged here, always, after it returns
        doTheAnalysis(programParams, reference, locusDescriptionCatalog, genomeQuery, mateExtractor, bamletWriter,
                      farAwayMateDistanceThreshold, typicalReadLength, jsonWriter, vcfWriter, {}, counts);
        logGenotypingSummary(programParams, counts);
        return;
    }

    const ReferenceContigInfo& contigInfo = reference.contigInfo();
    const int numContigs = std::max(1, contigInfo.numContigs());

    // PARALLEL PREPASS (--threads > 1): T workers each scan their owned whole contigs (chromosome stride
    // c % T == w; the prepass stride scans EVERY contig because a cross-contig mate makes any contig a
    // potential mate home) via the index, inserting far-away mates into a PRIVATE per-worker shard cache —
    // no shared writes, no locks. The shards are key-disjoint (each read's single primary alignment is
    // streamed by exactly one worker), so after every worker joins they are spliced into one immutable cache
    // on the main thread by mergeAndFreeze (which throws on the never-expected duplicate key).
    const int prepassWorkerCount = std::min(threadCount, numContigs);
    std::vector<std::unique_ptr<htshelpers::MateExtractor>> prepassShards;
    prepassShards.reserve(prepassWorkerCount);
    for (int w = 0; w < prepassWorkerCount; ++w)
    {
        prepassShards.push_back(std::make_unique<htshelpers::MateExtractor>(
            inputPaths.htsFile(), inputPaths.htsIndexFile(), inputPaths.reference(), true, farAwayMateDistanceThreshold));
    }
    spdlog::info("Caching read pairs with distant mates");
    {
        std::vector<std::exception_ptr> prepassExceptions(prepassWorkerCount);
        std::vector<std::thread> prepassWorkers;
        prepassWorkers.reserve(prepassWorkerCount);
        for (int w = 0; w < prepassWorkerCount; ++w)
        {
            prepassWorkers.emplace_back([&, w]() {
                try
                {
                    prepareCache(programParams, reference, genomeQuery, *prepassShards[w],
                        farAwayMateDistanceThreshold, typicalReadLength,
                        wholeContigRegionsForStride(contigInfo, w, threadCount));
                }
                catch (...)
                {
                    prepassExceptions[w] = std::current_exception();
                }
            });
        }
        for (std::thread& worker : prepassWorkers)
        {
            worker.join();
        }
        for (int w = 0; w < prepassWorkerCount; ++w)
        {
            if (prepassExceptions[w])
            {
                spdlog::error("Chromosome-stride prepass worker {} failed", w);
                std::rethrow_exception(prepassExceptions[w]);
            }
        }
    }

    // Freeze: splice the disjoint shards into one immutable cache on the main thread (happens-after every
    // prepass join via std::thread::join), then publish it read-only to the genotyping workers by value.
    std::vector<htshelpers::MateExtractor*> shardPtrs;
    shardPtrs.reserve(prepassShards.size());
    for (const std::unique_ptr<htshelpers::MateExtractor>& shard : prepassShards)
    {
        shardPtrs.push_back(shard.get());
    }
    std::shared_ptr<const htshelpers::MateCache> frozenCache = htshelpers::MateExtractor::mergeAndFreeze(shardPtrs);
    spdlog::info("Added {} reads to the mate cache", add_commas_at_thousands(frozenCache->size()));

    // PARALLEL GENOTYPING (--threads > 1): split the position-sorted catalog into one slice per loci-bearing
    // contig, genotype each on a chromosome-stride worker (its own FastaReference / MateExtractor over the
    // frozen cache / region-restricted streamer) into a per-contig temp JSON+VCF, then merge the temps in
    // ascending-contig-index order. Output is byte-identical to the threadCount==1 direct-write path: each
    // contig's loci are emitted in coordinate order into one temp file, concatenated in contig (== coordinate)
    // order.
    const std::vector<ContigGenotypingSlice> slices =
        perContigGenotypingSlices(locusDescriptionCatalog, contigInfo, farAwayMateDistanceThreshold);

    std::vector<std::string> jsonPaths;
    std::vector<std::string> vcfPaths;
    jsonPaths.reserve(slices.size());
    vcfPaths.reserve(slices.size());
    for (const ContigGenotypingSlice& slice : slices)
    {
        const std::string contigSuffix = ".contig" + std::to_string(slice.contigIndex);
        jsonPaths.push_back(programParams.outputPaths().outputPrefix() + contigSuffix + ".json");
        vcfPaths.push_back(programParams.outputPaths().outputPrefix() + contigSuffix + ".vcf");
    }

    // Remove the per-contig temp files on any exit path (success or exception). Ignores ENOENT.
    struct TempFileRemover
    {
        const std::vector<std::string>& jsonPaths;
        const std::vector<std::string>& vcfPaths;
        ~TempFileRemover()
        {
            for (const std::vector<std::string>* paths : { &jsonPaths, &vcfPaths })
            {
                for (const std::string& path : *paths)
                {
                    if (std::remove(path.c_str()) != 0 && errno != ENOENT)
                    {
                        spdlog::warn("Could not remove temporary file {}", path);
                    }
                }
            }
        }
    } tempFileRemover{ jsonPaths, vcfPaths };

    // Genotype the per-contig slices on chromosome-stride worker threads: worker w owns the contigs whose
    // index ≡ w (mod threadCount), so each loci-bearing contig is processed by exactly one worker, and each
    // worker writes only its own contigs' temp files (distinct paths → no file contention). Spawn at most
    // min(threadCount, numContigs) workers; a worker that owns no loci-bearing contig simply does nothing.
    // Each worker owns a private FastaReference + MateExtractor (the frozen cache is shared read-only by
    // value; localCache_ and the htsFile/bam1_t seek buffers are per-worker mutable and MUST NOT be shared).
    const int workerCount = std::min(threadCount, std::max(1, contigInfo.numContigs()));
    std::vector<std::exception_ptr> workerExceptions(workerCount);
    std::vector<GenotypingCounts> workerCounts(workerCount);  // per-worker tallies, summed after join
    std::vector<std::thread> workers;
    workers.reserve(workerCount);
    for (int w = 0; w < workerCount; ++w)
    {
        workers.emplace_back([&, w]() {
            try
            {
                for (std::size_t i = 0; i < slices.size(); ++i)
                {
                    const ContigGenotypingSlice& slice = slices[i];
                    if (slice.contigIndex % threadCount != w)
                    {
                        continue;  // owned by another stride worker
                    }
                    LocusDescriptionCatalog subsetCatalog(locusDescriptionCatalog.begin() + slice.catalogBegin,
                        locusDescriptionCatalog.begin() + slice.catalogEnd);
                    GenomeQueryCollection subsetGenomeQuery(subsetCatalog);
                    FastaReference workerReference(inputPaths.reference(), contigInfo);
                    htshelpers::MateExtractor workerMateExtractor(inputPaths.htsFile(), inputPaths.htsIndexFile(),
                        inputPaths.reference(), true, frozenCache, farAwayMateDistanceThreshold);
                    IterativeJsonWriter jsonWriter(programParams.sample(), workerReference.contigInfo(),
                        jsonPaths[i], programParams.copyCatalogFields());
                    IterativeVcfWriter vcfWriter(programParams.sample().id(), workerReference, vcfPaths[i]);
                    const std::vector<GenomicRegion> streamerRegions{ slice.region };
                    doTheAnalysis(programParams, workerReference, subsetCatalog, subsetGenomeQuery,
                        workerMateExtractor, bamletWriter, farAwayMateDistanceThreshold, typicalReadLength,
                        jsonWriter, vcfWriter, streamerRegions, workerCounts[w]);
                }
            }
            catch (...)
            {
                workerExceptions[w] = std::current_exception();
            }
        });
    }
    for (std::thread& worker : workers)
    {
        worker.join();
    }
    // Rethrow the first worker failure (after all joins) so no merge runs on partial output.
    for (int w = 0; w < workerCount; ++w)
    {
        if (workerExceptions[w])
        {
            spdlog::error("Chromosome-stride genotyping worker {} failed", w);
            std::rethrow_exception(workerExceptions[w]);
        }
    }

    spdlog::info("Merging {} per-contig output files into final JSON and VCF", slices.size());
    mergeRegionJsonFiles(programParams.outputPaths().json(), jsonPaths);
    mergeRegionVcfFiles(programParams.outputPaths().vcf(), vcfPaths);

    // One run-wide summary (each per-contig worker call stayed silent; sum their tallies here).
    GenotypingCounts aggregateCounts;
    for (const GenotypingCounts& workerTally : workerCounts)
    {
        aggregateCounts.fast += workerTally.fast;
        aggregateCounts.full += workerTally.full;
        aggregateCounts.skipped += workerTally.skipped;
        aggregateCounts.zeroCoverage += workerTally.zeroCoverage;
    }
    logGenotypingSummary(programParams, aggregateCounts);
}

}
