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

#include <climits>
#include <cstdlib>
#include <chrono>
#include <ctime>
#include <fstream>
#include <memory>
#include <unordered_map>

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
#include "io/MateCacheIO.hh"
#include "io/ParameterLoading.hh"
#include "io/SampleStats.hh"
#include "io/StringUtils.hh"
#include "io/VcfWriter.hh"
#include "locus/VariantFindings.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusSpecification.hh"
#include "locus/LocusAnalyzerUtil.hh"
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


const FullRead decodeRead(const htshelpers::HtsFileStreamer& readStreamer) {
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
};



// the UnpairedReadsCache is used to store reads for which the mate is yet to be found
using UnpairedReadsCache = std::unordered_map<FragmentId, FullRead>;


// Check if all variants in a locus are homozygous reference
bool isLocusHomRef(const LocusSpecification& locusSpec, const LocusFindings& locusFindings)
{
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        const VariantFindings* findings = variantIdAndFindings.second.get();

        // Check if it's a repeat variant
        const RepeatFindings* repeatFindings = dynamic_cast<const RepeatFindings*>(findings);
        if (repeatFindings != nullptr)
        {
            if (!repeatFindings->optionalGenotype())
            {
                // No genotype means we can't determine hom-ref status
                return false;
            }
            const auto& variantSpec = locusSpec.getVariantSpecById(variantId);
            const auto repeatNodeId = variantSpec.nodes().front();
            const auto& repeatUnit = locusSpec.regionGraph().nodeSeq(repeatNodeId);
            const int referenceSizeInUnits = variantSpec.referenceLocus().length() / repeatUnit.length();

            if (!isRepeatGenotypeHomRef(*repeatFindings->optionalGenotype(), referenceSizeInUnits))
            {
                return false;
            }
        }
        else
        {
            // Check if it's a small variant
            const SmallVariantFindings* smallVariantFindings = dynamic_cast<const SmallVariantFindings*>(findings);
            if (smallVariantFindings != nullptr)
            {
                if (!smallVariantFindings->optionalGenotype())
                {
                    // No genotype means we can't determine hom-ref status
                    return false;
                }
                if (!smallVariantFindings->optionalGenotype()->isHomRef())
                {
                    return false;
                }
            }
        }
    }
    return true;
}

void processLocus(const ProgramParameters& params, Reference& reference, const LocusDescription& locusDescription,
                  const std::vector<shared_ptr<FullReadPair>>& readPairs,
                  graphtools::AlignerSelector& alignerSelector,
                  BamletWriterPtr bamletWriter,
                  IterativeJsonWriter& jsonWriter,
                  IterativeVcfWriter& vcfWriter) {

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
    LocusAnalyzer locusAnalyzer(locusSpec, params.heuristics(), bamletWriter,
                                params.enableAlleleQualityMetrics());
    for (auto const& readPair : readPairs) {
        locusAnalyzer.processMates(readPair->firstMate->r, &readPair->secondMate->r, locus::RegionType::kTarget, alignerSelector);
    }

    // genotype the locus
    LocusFindings locusFindings = locusAnalyzer.analyze(
        sampleSex, boost::none, params.outputPaths().outputPrefix());

    // Check if --skip-hom-ref is enabled and all variants are hom-ref
    if (params.skipHomRef() && isLocusHomRef(locusSpec, locusFindings))
    {
        return;  // Skip output for hom-ref loci
    }

    jsonWriter.addRecord(locusSpec, locusFindings);

    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        vcfWriter.addRecord(variantId, locusSpec, locusFindings);
    }
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
                spdlog::error("The input read file contains reads from {} before {}, which is a different chromosome "
                    "order than the reference input FASTA file. Exiting..", chromosomeName, previousChromosomeName);
                return;
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

        const FullRead& fullRead = decodeRead(readStreamer);
        mateExtractor.addMateToCache(fullRead.r.readId(), fullRead);
    }

    spdlog::info("Added {} reads to the mate cache",  add_commas_at_thousands(mateExtractor.mateCacheSize()));
}


void doTheAnalysis(
    const ProgramParameters& params, Reference& reference, LocusDescriptionCatalog& locusDescriptionCatalog,
    const GenomeQueryCollection& genomeQuery, htshelpers::MateExtractor& mateExtractor,
    BamletWriterPtr bamletWriter, const int farAwayMateDistanceThreshold)
{
    if (locusDescriptionCatalog.empty()) {
        return;
    }

    int fastGenotypedCount = 0;
    int fullGenotypedCount = 0;
    int skippedCount = 0;

    boost::timer::progress_display progressBar(locusDescriptionCatalog.size());

    const int threadCount = params.threadCount;
    const unsigned htsDecompressionThreads(std::min(threadCount, 12));

    const InputPaths& inputPaths = params.inputPaths();
    const OutputPaths& outputPaths = params.outputPaths();
    IterativeJsonWriter jsonWriter(params.sample(), reference.contigInfo(), outputPaths.json(), params.copyCatalogFields());
    IterativeVcfWriter vcfWriter(params.sample().id(), reference, outputPaths.vcf());

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

    spdlog::info("Starting to process {} loci", add_commas_at_thousands(locusDescriptionCatalog.size()));
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

        const int32_t mateContigId = readStreamer.currentMateContigId();
        const int64_t mateStart = readStreamer.currentMatePosition();
        int64_t mateEnd = mateStart + typicalReadLength;

        const bool isMateFarAway = readContigId != mateContigId || std::abs(readStart - mateStart) >= farAwayMateDistanceThreshold;

        // handle readStreamer transitioning from one chromosome to the next
        if (previousReadContigId != readContigId) {
            const auto& chromosomeName = reference.contigInfo().getContigName(readContigId);

            // make sure chromosome ordering is consistent between the read input file and the reference fasta
            if (previousReadContigId != -1 && previousReadContigId > readContigId) {
                const auto& previousChromosomeName = reference.contigInfo().getContigName(previousReadContigId);
                spdlog::error("The input read file contains reads from {} before {}, which is a different chromosome "
                    "order than the reference input FASTA file. Exiting..", chromosomeName, previousChromosomeName);
                return;
            }
            previousReadContigId = readContigId;

            // load the next chromosome into the cache
            spdlog::info("Processing loci on {}", chromosomeName);
            reference.clearContigCache();
            reference.loadContigIntoCache(chromosomeName);
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
            ++progressBar;
            if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size()) {
                break;
            }
            nextLocusDescription = &locusDescriptionCatalog[nextLocusDescriptionIndex];
        }

        // fully parse the current read
        const FullRead& fullRead = decodeRead(readStreamer);

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

            if (readPair == nullptr) {
                readPair = std::make_shared<FullReadPair>(fullRead, fullMate);
            }

            locusCache = locusCachesMap[locusIndex];
            locusCache->readPairs.push_back(readPair);
        }

		// if all loci have been processed, break out of the loop
		if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size() && locusIndicesCollectingReads.empty() && locusIndicesReadyForAnalysis.empty()) {
            // done processing the locusDescriptionCatalog
            break;
        }

        // TODO clean up unpairedReadsCache, locusCaches, and mateExtractorCache?
        // every 500 bases, process and then destroy caches that are kExtensionLength to the left of the current read

        // process locusIndicesReadyForAnalysis
        if (locusIndicesReadyForAnalysis.size() > 0) {

            // process all active loci
            for (auto const& locusIndex : locusIndicesReadyForAnalysis) {
                //std::cout << "Processing locusCache " << locusIndex << std::endl;
                //delete the locusCache from the locusCachesMap

				if (locusCachesMap.find(locusIndex) == locusCachesMap.end()) {
					//a locus with 0 coverage won't have a locusCache. Skip analyzing it.
					//TODO handle missing / empty LocusCache?
					continue;
				}
				const shared_ptr<LocusCache> locusCache = locusCachesMap[locusIndex];

				// TODO check if # of reads >= minLocusCoverage and skip locus if not
				//if fast mode, check for active regions and downsample the reads
				bool needToProcessSlowly = true;
				if (params.analysisMode() == AnalysisMode::kOptimizedStreaming) {
                    const bool doneGenotyping = processLocusFast(params, reference,
                        locusDescriptionCatalog[locusIndex], locusCache->readPairs, jsonWriter, vcfWriter);

                    needToProcessSlowly = !doneGenotyping;
                    if (doneGenotyping) {
                        fastGenotypedCount++;
                    }

					//downsampleReads(locusCache->readPairs, locusDescriptionCatalog[locusIndex])
					//int readPairSignalCounter = 0;
					//for (auto const& readPair : locusCache->readPairs) {
					//     //check for active region
					//	bool readPairHasSignal = readPair->read.hasRepeatVariantSignal() || readPair->mate.hasRepeatVariantSignal();
					//	if (readPairHasSignal) {
					//		++readPairSignalCounter;
					//	}
					//}
					//const auto totalReadPairs = locusCache->readPairs.size();
					//if (readPairSignalCounter <= 2 || (readPairSignalCounter / (float) totalReadPairs) < 0.01) {
						//skip locus
					//	std::cout << "Skipping locus " << locusDescriptionCatalog[locusIndex].locusId << " because " << readPairSignalCounter << " out of " << totalReadPairs << std::setprecision(1) << " (" << 100.0 * readPairSignalCounter / (float) totalReadPairs << "%) reads had signal." << std::endl;
					//		continue;
					//} //else {
					//	std::cout << "Processing locus " << locusDescriptionCatalog[locusIndex].locusId << ": " << readPairSignalCounter << " out of " << totalReadPairs  << std::setprecision(1) << " (" << 100.0 * readPairSignalCounter / (float) totalReadPairs << "%) reads had signal." << std::endl;
					//}

					/*
					if (totalReadPairs > 10) {
						// compute coverage
						const int locusSize = locusDescriptionCatalog[locusIndex].locusAndFlanksEnd - locusDescriptionCatalog[locusIndex].locusAndFlanksStart;
						const int readLength = locusCache->readPairs[0]->read.sequence().length();
						const float coverage = 2 * totalReadPairs * readLength / (float) locusSize;
						if (coverage > 35) {
							std::cout << "Downsampling locus " << locusDescriptionCatalog[locusIndex].locusId << " with " << (int) coverage << "x coverage" << std::endl;
							//downsample to 35x coverage
							for (int k = totalReadPairs - 1; k >= 0; --k) {
								const float r = rand() / (float) RAND_MAX;
								if (r > 35 / coverage) {
									// delete current read pair
									locusCache->readPairs.erase(locusCache->readPairs.begin() + k);
								}
							}
							std::cout << "Kept " << locusCache->readPairs.size() << " out of " << totalReadPairs << " reads." << std::endl;
						}
					}
					*/
				}

				if (needToProcessSlowly) {
					if (params.heuristicGenotypingOnly()) {
						skippedCount++;
						jsonWriter.addSkippedRecord(locusDescriptionCatalog[locusIndex].locusId());
					} else {
						fullGenotypedCount++;
						try {
							graphtools::AlignerSelector alignerSelector(params.heuristics().alignerType());
							processLocus(params, reference, locusDescriptionCatalog[locusIndex], locusCache->readPairs,
								         alignerSelector, bamletWriter, jsonWriter, vcfWriter);
						} catch (const std::exception& e) {
							spdlog::error("Error while processing {}: {}", locusDescriptionCatalog[locusIndex].locusId(), e.what());
						}
					}
				}

				try {
					//TODO clear locusCache?
					locusCachesMap.erase(locusIndex);
				} catch (const std::exception& e) {
					spdlog::error("Unexpected error while erasing locus {}: {}. Skipping this locus.", locusDescriptionCatalog[locusIndex].locusId(), e.what());
				}

            }
            locusIndicesReadyForAnalysis.clear();
        }

        if (nextLocusDescriptionIndex >= locusDescriptionCatalog.size() && locusIndicesCollectingReads.empty() && locusIndicesReadyForAnalysis.empty()) {
            // done processing the locusDescriptionCatalog
            break;
        }
    }

    reference.clearContigCache();

    spdlog::info("Writing JSON and VCF output files");

    jsonWriter.close();
    vcfWriter.close();

    const int totalProcessed = fastGenotypedCount + fullGenotypedCount + skippedCount;
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

    htshelpers::MateExtractor mateExtractor(inputPaths.htsFile(), inputPaths.reference(), true, farAwayMateDistanceThreshold);

    const std::string mateCachePath = programParams.outputPaths().outputPrefix() + ".mate-cache.bam";

    bool cacheLoaded = false;

    // Try to load existing cache if --cache-mates is enabled
    if (programParams.cacheMates()) {
        cacheLoaded = htshelpers::MateCacheIO::readCache(
            mateCachePath,
            reference.contigInfo(),
            mateExtractor.mutableMateCache(),
            locusDescriptionCatalog.size());

        if (cacheLoaded) {
            spdlog::info("Loaded {} cached mates from {}",
                mateExtractor.mateCacheSize(), mateCachePath);
        }
    }

    // Only run prepareCache if we didn't load from file
    if (!cacheLoaded) {
        prepareCache(programParams, reference, genomeQuery, mateExtractor, farAwayMateDistanceThreshold);

        // Save cache if --cache-mates is enabled
        if (programParams.cacheMates()) {
            spdlog::info("Saving mate cache to {}", mateCachePath);
            htshelpers::MateCacheIO::writeCache(
                mateCachePath,
                reference.contigInfo(),
                mateExtractor.mateCache(),
                locusDescriptionCatalog.size());
        }
    }

    //perform the main bam/cram traversal and analysis
    doTheAnalysis(programParams, reference, locusDescriptionCatalog, genomeQuery, mateExtractor, bamletWriter,
                  farAwayMateDistanceThreshold);

}

}
