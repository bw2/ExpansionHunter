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
#include <deque>
#include <exception>
#include <fstream>
#include <future>
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
#include "core/ThreadPool.hh"
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


std::chrono::system_clock::rep microsecondsSinceEpoch() {
    auto now = std::chrono::system_clock::now().time_since_epoch();
    return std::chrono::duration_cast<std::chrono::microseconds>(now).count();
}


bool readOverlapsRepeat(const int32_t& contigId, const int64_t& readStart, const int64_t& readEnd, const LocusDescription& locusDescription) {
    if (contigId != locusDescription.locusContigIndex()) {
        return false;
    }

    //  read  |-----|
    //          |-----|

    return readStart <= locusDescription.locusWithoutFlanksEnd() && readEnd >= locusDescription.locusWithoutFlanksStart();
}


FullRead decodeRead(const htshelpers::HtsFileStreamer& readStreamer) {
    // fully parse the current read
    LinearAlignmentStats readAlignmentStats;
    Read read = readStreamer.decodeRead(readAlignmentStats);
    return {std::move(read), std::move(readAlignmentStats)};
}


struct LocusCache {
    LocusCache(int locusIndex, const GenomicRegion& locusInterval)
    : locusIndex(locusIndex), locusInterval(locusInterval) {}

    // Stores all reads + mates that will be analyzed for a given locus.
    int locusIndex;

    // the interval spanned by locus reference regions, including flanks
    GenomicRegion locusInterval;

    //a vector of pointers to FullReadPair objects
    std::vector<shared_ptr<FullReadPair>> readPairs;

    // fragment ids already added to readPairs, used to avoid double-counting a far-away pair whose
    // two ends are both contained in this locus window (each end is streamed in a separate iteration)
    std::unordered_set<FragmentId> seenFragmentIds;
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
                  const std::vector<shared_ptr<FullReadPair>>& readPairs,
                  graphtools::AlignerSelector& alignerSelector,
                  BamletWriterPtr bamletWriter) {

    LocusOutput out;
    out.locusIndex = locusIndex;
    out.locusId = locusDescription.locusId();

    try {
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
                                    params.enableAlleleQualityMetrics());
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
            Read mate = readPair->secondMate->r;

            const LinearAlignmentStats& readStats = readPair->firstMate->s;
            const LinearAlignmentStats& mateStats = readPair->secondMate->s;
            const int64_t readEnd = readStats.pos + static_cast<int64_t>(read.sequence().length());
            const int64_t mateEnd = mateStats.pos + static_cast<int64_t>(mate.sequence().length());

            const bool readIsInTargetRegion = readIsContainedInAnyExtractionRegion(
                targetExtractionRegions, readStats.chromId, readStats.pos, readEnd);
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


void prepareCache(
    const ProgramParameters& params, Reference& reference, const GenomeQueryCollection& genomeQuery,
    htshelpers::MateExtractor& mateExtractor, const int farAwayMateDistanceThreshold)
{

    spdlog::info("Caching read pairs with distant mates");
    const InputPaths& inputPaths = params.inputPaths();

    // keep track of the previous chromosome id to make sure that the read input file has the same chromosome ordering
    // as the reference fasta
    int previousReadContigId = -1;
    int typicalReadLength = -1;

    // print the number of threads
    const unsigned htsDecompressionThreads(std::min(params.threadCount, 12));
    htshelpers::HtsFileStreamer readStreamer(inputPaths.htsFile(), inputPaths.reference(), htsDecompressionThreads);
    while (readStreamer.tryReadingNextPrimaryAlignment() && readStreamer.isStreamingAlignedReads())
    {
        //skip unpaired reads
        if (!readStreamer.currentIsPaired()) {
            continue;
        }

        // get read length from the 1st read
        if (typicalReadLength == -1) {
            typicalReadLength = readStreamer.decodeRead().sequence().length();
        }

        const int32_t readContigId = readStreamer.currentReadContigId();
        const int64_t readStart = readStreamer.currentReadPosition();
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

        const int32_t mateContigId = readStreamer.currentMateContigId();
        const int64_t mateStart = readStreamer.currentMatePosition();
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

        FullRead fullRead = decodeRead(readStreamer);
        // Copy the ReadId out before the std::move below: passing both `fullRead.r.readId()` and
        // `std::move(fullRead)` as function args is use-after-move in C++17 (parameter init order is
        // indeterminate, and unlike std::pair's mem-init-list there is no ordering guarantee here).
        ReadId mateReadId = fullRead.r.readId();
        mateExtractor.addMateToCache(mateReadId, std::move(fullRead));
    }

    spdlog::info("Added {} reads to the mate cache",  add_commas_at_thousands(mateExtractor.mateCacheSize()));
}


void doTheAnalysis(
    const ProgramParameters& params, Reference& reference, LocusDescriptionCatalog& locusDescriptionCatalog,
    const GenomeQueryCollection& genomeQuery, htshelpers::MateExtractor& mateExtractor,
    BamletWriterPtr bamletWriter, const int farAwayMateDistanceThreshold,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter,
    const std::vector<GenomicRegion>& streamerRegions)
{
    if (locusDescriptionCatalog.empty()) {
        return;
    }

    int fastGenotypedCount = 0;
    int fullGenotypedCount = 0;
    int skippedCount = 0;
    int zeroCoverageCount = 0;

    // Constructed below only in whole-file mode (showProgress). In region-parallel mode P workers run this
    // function concurrently and a progress_display writes its layout to std::cout on construction, so one per
    // worker would garble the terminal. Left null when showProgress is false.
    std::unique_ptr<boost::timer::progress_display> progressBar;

    const int threadCount = params.threadCount;
    const unsigned htsDecompressionThreads(std::min(threadCount, 12));

    const InputPaths& inputPaths = params.inputPaths();

    // jsonWriter and vcfWriter are injected by the caller. Each sweep finalizes its own writer via the
    // jsonWriter.close()/vcfWriter.close() calls at the end of this function.

    // When streamerRegions is non-empty the read stream is restricted to those regions (region-parallel mode);
    // gating the progress bar on whole-file mode avoids garbled interleaved output across concurrent sweeps.
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
    int typicalReadLength = -1;

    // A single AlignerSelector (and thus a single PinnedDagAligner with its large DP score matrices)
    // reused across all full-genotyped loci on the serial path. The aligner is graph-agnostic — the graph
    // is supplied per align() call via the seed path — so reuse is output-neutral and keeps the matrix
    // buffers warm instead of reallocating ~hundreds of KB per locus.
    graphtools::AlignerSelector alignerSelector(params.heuristics().alignerType());

    // Parallel genotyping is enabled only for low-mem-streaming with >1 thread. optimized-streaming (which
    // also uses processLocusFast / --heuristic-genotyping-only) and single-thread runs keep the serial
    // path below. That path is structurally unchanged (same counters and catalog-ordered output), but it
    // now genotypes via genotypeLocusFull, which copies each Read before processMates instead of passing
    // the shared FullReadPair members by reference. Since align() mutates a Read in place (reverseComplement)
    // and one read pair can belong to multiple loci, this fixes a prior cross-locus mutation of shared reads
    // and so may change serial output for catalogs with overlapping loci relative to earlier releases.
    const bool useParallelGenotyping =
        (params.analysisMode() == AnalysisMode::kLowMemStreaming) && (threadCount > 1);

    // One AlignerSelector per worker thread, since its DP scratch matrices are mutated per align() call.
    std::vector<std::unique_ptr<graphtools::AlignerSelector>> threadAligners;
    if (useParallelGenotyping) {
        threadAligners.reserve(threadCount);
        for (int t = 0; t < threadCount; ++t) {
            threadAligners.push_back(
                std::make_unique<graphtools::AlignerSelector>(params.heuristics().alignerType()));
        }
    }

    // Ordered queue of pending per-locus outputs. FIFO order == ascending locusIndex == catalog/coordinate
    // order, so output is written deterministically regardless of which worker finishes first. An item is
    // either an in-flight future (full genotyping on the pool) or an already-built result (zero-coverage).
    struct PendingItem {
        bool isFuture = false;
        std::future<LocusOutput> future;
        LocusOutput ready;
    };
    std::deque<PendingItem> pending;
    // Cap in-flight + buffered loci to bound memory (low-mem mode); the pool keeps draining as we go.
    const size_t maxPending = static_cast<size_t>(threadCount) * 4;

    // Declared after threadAligners and pending so that on scope exit (including exception unwinding) the
    // pool is destroyed FIRST: ~thread_pool() runs and joins every queued/running task while the aligners
    // those tasks dereference are still alive. std::future destructors from the pool's packaged_tasks do
    // not block, so the pool itself — not the futures buffered in `pending` — is what joins the workers.
    ctpl::thread_pool genotypingPool(useParallelGenotyping ? threadCount : 0);

    // Write one finished locus result. Single-threaded: called inline on the serial path, and only from
    // the main thread while draining the pending queue on the parallel path. Counter semantics match the
    // previous serial code (a locus reaching full genotyping is counted as full-genotyped even if it
    // produces no record because of --skip-hom-ref / --skip-missing-genotypes or an error).
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

    // Drain the pending queue in FIFO order until at most `keep` items remain. keep==0 is a full barrier
    // (used before evicting the reference contig cache and at EOF).
    auto drainPending = [&](size_t keep) {
        while (pending.size() > keep) {
            PendingItem item = std::move(pending.front());
            pending.pop_front();
            LocusOutput out = item.isFuture ? item.future.get() : std::move(item.ready);
            writeOutput(out);
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
                if (useParallelGenotyping) {
                    PendingItem item;
                    item.ready.locusIndex = locusIndex;
                    item.ready.kind = LocusOutput::Kind::kNoCoverage;
                    item.ready.locusId = locusDescriptionCatalog[locusIndex].locusId();
                    pending.push_back(std::move(item));
                } else {
                    zeroCoverageCount++;
                    writeZeroCoverageRecord(params, reference, locusDescriptionCatalog[locusIndex], jsonWriter, vcfWriter);
                }
                continue;
            }
            const shared_ptr<LocusCache> locusCache = locusCachesMap[locusIndex];

            // TODO check if # of reads >= minLocusCoverage and skip locus if not
            bool needToProcessSlowly = true;
            if (params.analysisMode() == AnalysisMode::kOptimizedStreaming) {
                const bool doneGenotyping = processLocusFast(params, reference,
                    locusDescriptionCatalog[locusIndex], locusCache->readPairs, jsonWriter, vcfWriter);

                needToProcessSlowly = !doneGenotyping;
                if (doneGenotyping) {
                    fastGenotypedCount++;
                }
            }

            if (needToProcessSlowly) {
                if (params.heuristicGenotypingOnly()) {
                    if (useParallelGenotyping) {
                        // Route through the ordered pending queue (like kNoCoverage) so this skipped
                        // record is emitted in catalog order relative to other buffered records, rather
                        // than written directly ahead of lower-index records still buffered in `pending`.
                        PendingItem item;
                        item.ready.locusIndex = locusIndex;
                        item.ready.kind = LocusOutput::Kind::kHeuristicOnlySkip;
                        item.ready.locusId = locusDescriptionCatalog[locusIndex].locusId();
                        pending.push_back(std::move(item));
                    } else {
                        skippedCount++;
                        jsonWriter.addSkippedRecord(locusDescriptionCatalog[locusIndex].locusId(), "heuristic_only_mode");
                    }
                } else if (useParallelGenotyping) {
                    // Dispatch full genotyping to the pool. locusCache is captured by value (shared_ptr),
                    // so the worker keeps it alive even after the map entry is erased below.
                    //
                    // Note: the JSON/VCF records are buffered in `pending` and written in catalog order, but
                    // the bamlet (--enable-bamlet-output) records are written directly from the worker via
                    // LocusAligner::align -> BamletWriter::write, so their order in <prefix>_realigned.bam is
                    // worker-completion order, i.e. nondeterministic across runs in parallel mode. This is
                    // intentional: the bamlet is SO:unknown and consumers (e.g. REViewer) match records by
                    // read name + the XG graph-alignment tag, not file order. The serial path writes them in
                    // locus order. To make the bamlet deterministic in parallel mode, buffer each worker's
                    // realigned reads in the LocusOutput and replay them from writeOutput in catalog order.
                    PendingItem item;
                    item.isFuture = true;
                    item.future = genotypingPool.push(
                        [&params, &reference, &locusDescriptionCatalog, locusIndex, locusCache, bamletWriter,
                         &threadAligners](int threadId) {
                            return genotypeLocusFull(params, reference, locusIndex,
                                locusDescriptionCatalog[locusIndex], locusCache->readPairs,
                                *threadAligners[threadId], bamletWriter);
                        });
                    pending.push_back(std::move(item));
                } else {
                    LocusOutput out = genotypeLocusFull(params, reference, locusIndex,
                        locusDescriptionCatalog[locusIndex], locusCache->readPairs, alignerSelector, bamletWriter);
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

        if (useParallelGenotyping) {
            drainPending(maxPending);
        }
    };

    spdlog::info("Starting to process {} loci", add_commas_at_thousands(locusDescriptionCatalog.size()));
    // In whole-file mode use all decompression threads; in region mode use a single one to avoid
    // P-way x htsDecompressionThreads bgzf thread oversubscription across concurrent region sweeps.
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

        // get read length from the 1st read
        if (typicalReadLength == -1) {
            typicalReadLength = readStreamer->decodeRead().sequence().length();
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

            // Before evicting the previous contig's sequence from the reference cache, finish genotyping
            // every locus on it. Workers read the cached contig inside decodeLocusSpecification, so the
            // cache must stay valid until they complete. Once a read on the next contig arrives, all
            // loci still collecting reads belong to the previous contig (coordinate order), so flush them
            // here exactly as the per-read loop below would, then fully drain the pool.
            if (useParallelGenotyping) {
                // Flush unconditionally. Besides loci still collecting reads, locusIndicesReadyForAnalysis
                // can hold loci left unprocessed by a prior per-read `continue`; the old `if (!collecting
                // .empty())` guard skipped those, so they were genotyped only after the contig was evicted
                // — racing in the unsynchronized faidx cache-miss path. (processReadyLoci no-ops on empty.)
                for (auto idx : locusIndicesCollectingReads) {
                    locusIndicesReadyForAnalysis.push_back(idx);
                }
                locusIndicesCollectingReads.clear();
                processReadyLoci();
                drainPending(0);
            }

            // load the next chromosome into the cache
            spdlog::info("Processing loci on {}", chromosomeName);
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

        // Since both the read & mate are now available, double-check that the readEnd and mateEnd coordinates were
        // computed correctly and that the read length parsed from the 1st read in the BAM file was not atypical
        const int readSequenceLength = (int) fullRead.r.sequence().length();
        if (readSequenceLength != typicalReadLength) {
            readEnd = readStart + readSequenceLength;
        }
        const int mateSequenceLength = (int) fullMate.r.sequence().length();
        if (mateSequenceLength != typicalReadLength) {
            mateEnd = mateStart + mateSequenceLength;
        }
        if (typicalReadLength < readSequenceLength || typicalReadLength < mateSequenceLength) {
            typicalReadLength = std::max(readSequenceLength, mateSequenceLength);
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
                    });
                locusCachesMap[locusIndex] = locusCache;
            }

            locusCache = locusCachesMap[locusIndex];

            // skip if this fragment's pair was already added to this locus (e.g. a far-away pair whose
            // upstream and downstream ends are both contained in this locus window are streamed, and
            // processed, in two separate iterations)
            if (!locusCache->seenFragmentIds.insert(fullRead.r.fragmentId()).second) {
                continue;
            }

            if (readPair == nullptr) {
                readPair = std::make_shared<FullReadPair>(fullRead, fullMate);
            }

            locusCache->readPairs.push_back(readPair);
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

    // Wait for all in-flight workers and write their results before tearing down the pool / reference.
    if (useParallelGenotyping) {
        drainPending(0);
    }

    reference.clearContigCache();

    spdlog::info("Writing JSON and VCF output files");

    jsonWriter.close();
    vcfWriter.close();

    const int totalProcessed = fastGenotypedCount + fullGenotypedCount + skippedCount + zeroCoverageCount;
    if (totalProcessed > 0) {
        const float fastPct = 100.0f * fastGenotypedCount / totalProcessed;
        const float fullPct = 100.0f * fullGenotypedCount / totalProcessed;
        if (skippedCount > 0) {
            const float skippedPct = 100.0f * skippedCount / totalProcessed;
            spdlog::info("Processed {} ({:.1f}%) loci via fast genotyping, {} ({:.1f}%) loci via full genotyping, and skipped {} ({:.1f}%) loci",
                add_commas_at_thousands(fastGenotypedCount), fastPct,
                add_commas_at_thousands(fullGenotypedCount), fullPct,
                add_commas_at_thousands(skippedCount), skippedPct);
        } else {
            spdlog::info("Processed {} ({:.1f}%) loci via fast genotyping, and {} ({:.1f}%) loci via full genotyping",
                add_commas_at_thousands(fastGenotypedCount), fastPct,
                add_commas_at_thousands(fullGenotypedCount), fullPct);
        }
        if (zeroCoverageCount > 0) {
            const float zeroCoveragePct = 100.0f * zeroCoverageCount / totalProcessed;
            if (params.skipMissingGenotypes()) {
                spdlog::info("Excluded {} ({:.1f}%) loci with zero coverage from output (--skip-missing-genotypes)",
                    add_commas_at_thousands(zeroCoverageCount), zeroCoveragePct);
            } else {
                spdlog::info("Emitted {} ({:.1f}%) no-call records for loci with zero coverage",
                    add_commas_at_thousands(zeroCoverageCount), zeroCoveragePct);
            }
        }
    }
}



void htsLowMemStreamingSampleAnalysis(
    LocusDescriptionCatalog& locusDescriptionCatalog,
    const ProgramParameters& programParams,
    Reference& reference,
    BamletWriterPtr bamletWriter)
{
    if (locusDescriptionCatalog.empty()) {
        return;
    }

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

    htshelpers::MateExtractor mateExtractor(inputPaths.htsFile(), inputPaths.htsIndexFile(), inputPaths.reference(), true, farAwayMateDistanceThreshold);

    prepareCache(programParams, reference, genomeQuery, mateExtractor, farAwayMateDistanceThreshold);

    //perform the main bam/cram traversal and analysis
    IterativeJsonWriter jsonWriter(programParams.sample(), reference.contigInfo(), programParams.outputPaths().json(), programParams.copyCatalogFields());
    IterativeVcfWriter vcfWriter(programParams.sample().id(), reference, programParams.outputPaths().vcf());
    doTheAnalysis(programParams, reference, locusDescriptionCatalog, genomeQuery, mateExtractor, bamletWriter,
                  farAwayMateDistanceThreshold, jsonWriter, vcfWriter, {});

}


// Region-parallel variant of the low-mem streaming analysis. The coordinate-sorted catalog is split into
// P contiguous, equal-locus-count slices (one per worker). A single serial prepass builds the far-away
// mate cache once; it is then frozen and shared read-only across all workers. Each worker runs a private
// full sweep restricted to its slice's genomic span (+/- a farAwayMateDistanceThreshold halo) with its own
// file handle / FastaReference / MateExtractor, writing a per-region temp JSON+VCF. After all workers join,
// the temp files are merged in genomic (ascending region) order into the final outputs. Intra-worker
// genotyping is serial because the analysis mode is kRegionParallelStreaming (not kLowMemStreaming).
void htsRegionParallelStreamingSampleAnalysis(
    LocusDescriptionCatalog& locusDescriptionCatalog,
    const ProgramParameters& programParams,
    Reference& reference,
    BamletWriterPtr bamletWriter)
{
    if (locusDescriptionCatalog.empty()) {
        return;
    }

    const int farAwayMateDistanceThreshold = 1000;  // base pairs (also used as the per-side streamer halo)
    const int N = (int) locusDescriptionCatalog.size();
    const int P = std::min(std::max(1, programParams.threadCount), N);

    spdlog::info("Running region-parallel-streaming analysis with {} workers over {} loci", P, N);

    // Capture the contig info by const-ref before spawning threads so workers never touch the shared
    // `reference` object; each worker builds its own FastaReference from this contig info instead.
    const ReferenceContigInfo& contigInfo = reference.contigInfo();
    const InputPaths& inputPaths = programParams.inputPaths();

    // Serial prepass: build the shared far-away mate cache once, then freeze it for read-only sharing.
    GenomeQueryCollection fullGenomeQuery(locusDescriptionCatalog);
    htshelpers::MateExtractor cacheBuilder(
        inputPaths.htsFile(), inputPaths.htsIndexFile(), inputPaths.reference(), true, farAwayMateDistanceThreshold);
    prepareCache(programParams, reference, fullGenomeQuery, cacheBuilder, farAwayMateDistanceThreshold);
    std::shared_ptr<const htshelpers::MateCache> sharedCache = cacheBuilder.freezeAndShareCache();

    // Partition the catalog into P contiguous equal-locus slices using floor boundaries. All slices are
    // non-empty because P <= N. For each slice also compute its per-contig streamer regions and temp paths.
    std::vector<LocusDescriptionCatalog> slices;
    std::vector<std::vector<GenomicRegion>> sliceRegions;
    std::vector<std::string> jsonPaths;
    std::vector<std::string> vcfPaths;
    slices.reserve(P);
    sliceRegions.reserve(P);
    jsonPaths.reserve(P);
    vcfPaths.reserve(P);

    for (int i = 0; i < P; ++i) {
        const int lo = (int) ((long long) i * N / P);
        const int hi = (int) ((long long) (i + 1) * N / P);
        slices.emplace_back(locusDescriptionCatalog.begin() + lo, locusDescriptionCatalog.begin() + hi);
        const LocusDescriptionCatalog& slice = slices.back();

        // One GenomicRegion per contig the slice touches. The slice is coordinate-sorted, so same-contig
        // loci are contiguous; collapse each contig's loci into [minStart - halo, maxEnd + halo] clamped to
        // [0, contigLen].
        std::vector<GenomicRegion> regionsForSlice;
        size_t k = 0;
        while (k < slice.size()) {
            const int32_t c = slice[k].locusContigIndex();
            int64_t minStart = slice[k].locusAndFlanksStart();
            int64_t maxEnd = slice[k].locusAndFlanksEnd();
            size_t m = k;
            while (m < slice.size() && slice[m].locusContigIndex() == c) {
                minStart = std::min(minStart, slice[m].locusAndFlanksStart());
                maxEnd = std::max(maxEnd, slice[m].locusAndFlanksEnd());
                ++m;
            }
            const int64_t beg = std::max<int64_t>(0, minStart - farAwayMateDistanceThreshold);
            const int64_t end = std::min<int64_t>(contigInfo.getContigSize(c), maxEnd + farAwayMateDistanceThreshold);
            regionsForSlice.emplace_back(c, beg, end);
            k = m;
        }
        sliceRegions.push_back(std::move(regionsForSlice));

        jsonPaths.push_back(programParams.outputPaths().outputPrefix() + ".region" + std::to_string(i) + ".json");
        vcfPaths.push_back(programParams.outputPaths().outputPrefix() + ".region" + std::to_string(i) + ".vcf");
    }

    // Best-effort removal of the per-region temp files, as an RAII guard so it runs on every exit path: the
    // normal return after the merge below, and stack unwinding when a worker failure is rethrown. ENOENT is
    // ignored because a worker that failed early may never have created its temp file.
    struct TempFileRemover
    {
        const std::vector<std::string>& jsonPaths;
        const std::vector<std::string>& vcfPaths;
        ~TempFileRemover()
        {
            for (const std::vector<std::string>* paths : { &jsonPaths, &vcfPaths }) {
                for (const std::string& path : *paths) {
                    if (std::remove(path.c_str()) != 0 && errno != ENOENT) {
                        spdlog::warn("Could not remove temporary file {}", path);
                    }
                }
            }
        }
    } tempFileRemover{ jsonPaths, vcfPaths };

    // Spawn P worker threads. Each worker only reads its OWN index i of the shared vectors (no aliasing).
    // bamletWriter is internally thread-safe; sharedCache is an immutable shared_ptr<const>; contigInfo is
    // a const reference captured before any threads start. Any exception thrown by a worker is captured into
    // workerExceptions[i] and rethrown on the main thread after join.
    std::vector<std::exception_ptr> workerExceptions(P);
    std::vector<std::thread> workers;
    workers.reserve(P);
    for (int i = 0; i < P; ++i) {
        workers.emplace_back([&, i]() {
            try {
                FastaReference workerReference(inputPaths.reference(), contigInfo);
                htshelpers::MateExtractor workerMateExtractor(
                    inputPaths.htsFile(), inputPaths.htsIndexFile(), inputPaths.reference(), true, sharedCache,
                    farAwayMateDistanceThreshold);
                GenomeQueryCollection sliceGenomeQuery(slices[i]);
                IterativeJsonWriter jsonWriter(
                    programParams.sample(), workerReference.contigInfo(), jsonPaths[i],
                    programParams.copyCatalogFields());
                IterativeVcfWriter vcfWriter(programParams.sample().id(), workerReference, vcfPaths[i]);
                doTheAnalysis(programParams, workerReference, slices[i], sliceGenomeQuery, workerMateExtractor,
                              bamletWriter, farAwayMateDistanceThreshold, jsonWriter, vcfWriter, sliceRegions[i]);
            } catch (...) {
                workerExceptions[i] = std::current_exception();
            }
        });
    }

    for (auto& worker : workers) {
        worker.join();
    }

    // Propagate the first worker failure (if any) to main's catch before merging or cleaning up.
    for (int i = 0; i < P; ++i) {
        if (workerExceptions[i]) {
            spdlog::error("Region-parallel-streaming worker {} failed", i);
            std::rethrow_exception(workerExceptions[i]);
        }
    }

    // Merge per-region temp files in ascending (== genomic) order into the final outputs.
    spdlog::info("Merging {} per-region output files into final JSON and VCF", P);
    mergeRegionJsonFiles(programParams.outputPaths().json(), jsonPaths);
    mergeRegionVcfFiles(programParams.outputPaths().vcf(), vcfPaths);

    // The per-region temp files are removed by tempFileRemover (RAII) on scope exit.
}

}
