//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Chris Saunders <csaunders@illumina.com>
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

#include "sample/HtsSeekingSampleAnalysis.hh"

#include <atomic>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

// clang-format off
// Note that spdlog.h must be included before ostr.h
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
// clang-format on

#include "core/Read.hh"
#include "core/ReadPairs.hh"
#include "locus/LocusAnalyzer.hh"
#include "sample/AnalyzerFinder.hh"
#include "sample/HtsFileSeeker.hh"
#include "sample/IndexBasedDepthEstimate.hh"
#include "sample/MateExtractor.hh"

using boost::optional;
using ehunter::htshelpers::HtsFileSeeker;
using ehunter::locus::LocusAnalyzer;
using std::ostream;
using std::string;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

namespace ehunter
{

namespace
{

vector<GenomicRegion>
combineRegions(const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions)
{
    vector<GenomicRegion> combinedRegions(targetRegions);
    combinedRegions.insert(combinedRegions.end(), offtargetRegions.begin(), offtargetRegions.end());
    return combinedRegions;
}

bool checkIfMatesWereMappedNearby(const LinearAlignmentStats& alignmentStats)
{
    const int kMaxMateDistance = 1000;
    if ((alignmentStats.chromId == alignmentStats.mateChromId)
        && (std::abs(alignmentStats.pos - alignmentStats.matePos) < kMaxMateDistance))
    {
        return true;
    }
    return false;
}

bool addMateToRegionIfNearby(
    GenomicRegion& mateGenomicRegion, const ReadId& readId, htshelpers::MateRegionToRecover& mateRegionToRecover)
{
    const int maxDistance = 2000;
    bool isNearby = mateRegionToRecover.genomicRegion.distance(mateGenomicRegion) < maxDistance;

    ReadId mateReadId(readId.fragmentId(),
        (readId.mateNumber() == MateNumber::kFirstMate) ? MateNumber::kSecondMate : MateNumber::kFirstMate);

    if (isNearby) {
        mateRegionToRecover.genomicRegion = GenomicRegion(
            mateRegionToRecover.genomicRegion.contigIndex(),
            std::min(mateRegionToRecover.genomicRegion.start(), mateGenomicRegion.start()),
            std::max(mateRegionToRecover.genomicRegion.end(), mateGenomicRegion.end()));

        mateRegionToRecover.mateReadIds.emplace(mateReadId);
    }

    return isNearby;
}

void recoverMates(
    htshelpers::MateExtractor& mateExtractor, FullReadPairs& readPairs)
{

    // compute a vector of MateRegionToRecover structs, each of which represents a GenomicRegion + a vector of
    // ReadIds that need to be recovered from that region.
    vector<htshelpers::MateRegionToRecover> mateRegionsToRecover;

    for (auto& fragmentIdAndReadPair : readPairs)
    {
        FullReadPair& readPair = fragmentIdAndReadPair.second;

        if (readPair.numMatesSet() == 2)
        {
            continue;
        }

        const FullRead& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;

        if (checkIfMatesWereMappedNearby(read.s)) {
            continue;
        }

        const int32_t mateContigIndex = read.s.isMateMapped ? read.s.mateChromId : read.s.chromId;
        const int32_t matePos = read.s.isMateMapped ? read.s.matePos : read.s.pos;

        // check if mate is close to other mates so they can be fetched together
        bool foundNearbyRegion = false;
        for (auto& mateRegionToRecover : mateRegionsToRecover) {
            GenomicRegion mateGenomicRegion = GenomicRegion(mateContigIndex, matePos, matePos + read.r.sequence().length());

            if (addMateToRegionIfNearby(mateGenomicRegion, read.r.readId(), mateRegionToRecover)) {
                foundNearbyRegion = true;
                break;
            }
        }
        if (!foundNearbyRegion) {
            // create a new entry in mateRegionsToRecover
            GenomicRegion mateGenomicRegion = GenomicRegion(mateContigIndex, matePos, matePos + read.r.sequence().length());
            htshelpers::MateRegionToRecover newMateRegionToRecover = { mateGenomicRegion, std::unordered_set<ReadId, boost::hash<ReadId>>() };
            addMateToRegionIfNearby(mateGenomicRegion, read.r.readId(), newMateRegionToRecover);

            mateRegionsToRecover.push_back(newMateRegionToRecover);
        }
    }

    // fetch mates for the regions
    for (auto& mateRegionToRecover : mateRegionsToRecover) {
        for (auto& mate : mateExtractor.extractMates(mateRegionToRecover)) {
            readPairs.AddMateToExistingRead(mate);
        }
    }
}


FullReadPairs collectCandidateReads(
    const vector<GenomicRegion>& targetRegions, const vector<GenomicRegion>& offtargetRegions,
    HtsFileSeeker& htsFileSeeker, htshelpers::MateExtractor& mateExtractor)
{
    vector<GenomicRegion> regionsWithReads = combineRegions(targetRegions, offtargetRegions);
    FullReadPairs readPairs;

    for (const auto& regionWithReads : regionsWithReads)
    {
        const int numReadsBeforeCollection = readPairs.NumReads();
        htsFileSeeker.setRegion(regionWithReads);
        while (htsFileSeeker.trySeekingToNextPrimaryAlignment())
        {
            LinearAlignmentStats alignmentStats;
            Read read = htsFileSeeker.decodeRead(alignmentStats);
            if (alignmentStats.isPaired)
            {
                FullRead fullRead(std::move(read), std::move(alignmentStats));
                readPairs.Add(std::move(fullRead));
            }
            else
            {
                spdlog::warn("Skipping {} because it is unpaired", read.readId().toString());
            }
        }
        const int numReadsCollected = readPairs.NumReads() - numReadsBeforeCollection;
        spdlog::debug("Collected {} reads from {}", numReadsCollected, regionWithReads.toString());
    }

    const int numReadsBeforeRecovery = readPairs.NumReads();
    recoverMates(mateExtractor, readPairs);
    const int numReadsAfterRecovery = readPairs.NumReads() - numReadsBeforeRecovery;
    spdlog::debug("Recovered {} reads", numReadsAfterRecovery);

    return readPairs;
}

void analyzeReadPair(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, AnalyzerFinder& analyzerFinder, FullRead& read, FullRead& mate,
    graphtools::AlignerSelector& alignerSelector)
{


    const int64_t readEnd = read.s.pos + read.r.sequence().length();
    const int64_t mateEnd = mate.s.pos + mate.r.sequence().length();
    vector<AnalyzerBundle> analyzers
        = analyzerFinder.query(read.s.chromId, read.s.pos, readEnd, mate.s.chromId, mate.s.pos, mateEnd);

    if (analyzers.empty())
    {
        return;
    }

    assert(analyzers.size() == 1);
    auto& analyzer(analyzers.front());
    processAnalyzerBundleReadPair(
        *locusAnalyzers[analyzer.locusIndex], analyzer.regionType, analyzer.inputType, read.r, mate.r, alignerSelector);
}

void analyzeRead(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, AnalyzerFinder& analyzerFinder, FullRead& read,
    graphtools::AlignerSelector& alignerSelector)
{
    const int64_t readEnd = read.s.pos + read.r.sequence().length();

    vector<AnalyzerBundle> analyzers = analyzerFinder.query(read.s.chromId, read.s.pos, readEnd);

    if (analyzers.empty())
    {
        return;
    }

    assert(analyzers.size() == 1);
    auto& analyzer(analyzers.front());
    locusAnalyzers[analyzer.locusIndex]->processMates(read.r, nullptr, analyzer.regionType, alignerSelector);
}

void processReads(
    vector<unique_ptr<LocusAnalyzer>>& locusAnalyzers, FullReadPairs& candidateReadPairs,
    AnalyzerFinder& analyzerFinder, graphtools::AlignerSelector& alignerSelector)
{
    for (auto& fragmentIdAndReads : candidateReadPairs)
    {
        auto& readPair = fragmentIdAndReads.second;
        if (readPair.numMatesSet() == 2)
        {
            analyzeReadPair(
                locusAnalyzers, analyzerFinder, *readPair.firstMate, *readPair.secondMate,
                alignerSelector);
        }
        else
        {
            FullRead& read = readPair.firstMate ? *readPair.firstMate : *readPair.secondMate;
            analyzeRead(locusAnalyzers, analyzerFinder, read, alignerSelector);
        }
    }
}

/// \brief Mutable data shared by all worker threads
///
class LocusThreadSharedData
{
public:
    LocusThreadSharedData()
        : isWorkerThreadException(false)
        , locusIndex(0)
    {
    }

    std::atomic<bool> isWorkerThreadException;
    std::atomic<unsigned> locusIndex;
};

/// \brief Data isolated to each locus-processing thread
///
struct LocusThreadLocalData
{
    std::exception_ptr threadExceptionPtr = nullptr;
};

/// \brief Process a series of loci on one thread
///
void processLocus(
    const int threadIndex, const ProgramParameters& programParams,
    const HeuristicParameters& heuristicParams, const RegionCatalog& regionCatalog,
    BamletWriterPtr alignmentWriter, SampleFindings& sampleFindings, LocusThreadSharedData& locusThreadSharedData,
    std::vector<LocusThreadLocalData>& locusThreadLocalDataPool)
{
    const InputPaths& inputPaths = programParams.inputPaths();
    const Sex& sampleSex = programParams.sample().sex();

    LocusThreadLocalData& locusThreadData(locusThreadLocalDataPool[threadIndex]);
    std::string locusId = "Unknown";

    try
    {
        HtsFileSeeker htsFileSeeker(inputPaths.htsFile(), inputPaths.reference());
        htshelpers::MateExtractor mateExtractor(inputPaths.htsFile(), inputPaths.reference(), programParams.cacheMates());
        graphtools::AlignerSelector alignerSelector(heuristicParams.alignerType());

        const unsigned size(regionCatalog.size());
        while (true)
        {
            if (locusThreadSharedData.isWorkerThreadException.load())
            {
                return;
            }
            const auto locusIndex(locusThreadSharedData.locusIndex.fetch_add(1));
            if (locusIndex >= size)
            {
                return;
            }

            const auto& locusSpec(regionCatalog[locusIndex]);
            locusId = locusSpec.locusId();

            std::string locusMotif = "";
            const auto& graph = locusSpec.regionGraph();
            for (const auto& variantSpec : locusSpec.variantSpecs()) {
                const int repeatNodeId = static_cast<int>(variantSpec.nodes().front());
                const auto& motif = graph.nodeSeq(repeatNodeId);
                if (locusMotif.length() > 0) {
                    locusMotif += "/";
                }
                locusMotif += motif;
            }

            spdlog::info("Analyzing locus {}: {}", locusIndex + 1, locusId);
            vector<unique_ptr<LocusAnalyzer>> locusAnalyzers;
            auto analyzer(std::make_unique<LocusAnalyzer>(locusSpec, heuristicParams, alignmentWriter,
                                                          programParams.enableAlleleQualityMetrics()));
            locusAnalyzers.emplace_back(std::move(analyzer));
            AnalyzerFinder analyzerFinder(locusAnalyzers);

            FullReadPairs readPairs = collectCandidateReads(
                locusSpec.targetReadExtractionRegions(), locusSpec.offtargetReadExtractionRegions(),
                htsFileSeeker, mateExtractor);

            processReads(locusAnalyzers, readPairs, analyzerFinder, alignerSelector);
            sampleFindings[locusIndex] = locusAnalyzers.front()->analyze(
                sampleSex, boost::none,
                programParams.outputPaths().outputPrefix());
        }
    }
    catch (const std::exception& e)
    {
        locusThreadSharedData.isWorkerThreadException = true;
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Exception caught in thread {} while processing locus: {} : {}", threadIndex, locusId, e.what());
        throw;
    }
    catch (...)
    {
        locusThreadSharedData.isWorkerThreadException = true;
        locusThreadData.threadExceptionPtr = std::current_exception();

        spdlog::error("Unknown exception caught in thread {} while processing locus: {}", threadIndex, locusId);
        throw;
    }
}
}

SampleFindings htsSeekingSampleAnalysis(
    const ProgramParameters& programParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, BamletWriterPtr alignmentWriter)
{
    const int threadCount = programParams.threadCount;
    if (threadCount > 1)
    {
        const InputPaths& inputPaths = programParams.inputPaths();
        if (ehunter::isURL(inputPaths.htsFile()))
        {
            // For URL input paths, the index needs to be downloaded in advance if seeking mode is using multiple threads.
            // This is needed because htslib has no protection against the race condition created by multiple threads
            // independently downloading this index to the same file path.
            //
            (void)HtsFileSeeker(inputPaths.htsFile(), inputPaths.reference());
        }
    }

    LocusThreadSharedData locusThreadSharedData;
    std::vector<LocusThreadLocalData> locusThreadLocalDataPool(threadCount);

    const unsigned locusCount(regionCatalog.size());
    SampleFindings sampleFindings(locusCount);

    // Start all locus worker threads
    std::vector<std::thread> locusThreads;
    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        locusThreads.emplace_back(
            processLocus, threadIndex, std::cref(programParams), std::cref(heuristicParams),
            std::cref(regionCatalog), alignmentWriter, std::ref(sampleFindings), std::ref(locusThreadSharedData),
            std::ref(locusThreadLocalDataPool));
    }

    // Rethrow exceptions from worker pool in thread order:
    if (locusThreadSharedData.isWorkerThreadException.load())
    {
        for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
        {
            const auto& locusThreadData(locusThreadLocalDataPool[threadIndex]);
            if (locusThreadData.threadExceptionPtr)
            {
                std::rethrow_exception(locusThreadData.threadExceptionPtr);
            }
        }
    }

    for (int threadIndex(0); threadIndex < threadCount; ++threadIndex)
    {
        locusThreads[threadIndex].join();
    }

    return sampleFindings;
}

}
