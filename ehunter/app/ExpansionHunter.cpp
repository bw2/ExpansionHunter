//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// clang-format off
// Note that spdlog.h must be included before ostr.h
#include "spdlog/spdlog.h"
//#include "thirdparty/spdlog/include/spdlog/fmt/ostr.h"
// clang-format on

#include "app/Version.hh"
#include "core/Common.hh"
#include "core/Parameters.hh"
#include "io/BamletWriter.hh"
#include "io/CatalogLoading.hh"
#include "io/GraphBlueprint.hh"
#include "io/StringUtils.hh"
#include "locus/LocusSpecification.hh"
#include "io/JsonWriter.hh"
#include "io/ParameterLoading.hh"
#include "io/SampleStats.hh"
#include "io/VcfWriter.hh"
#include "locus/VariantFindings.hh"
#include "sample/HtsLowMemStreamingHelpers.hh"
#include "sample/HtsLowMemStreamingSampleAnalysis.hh"
#include "sample/HtsSeekingSampleAnalysis.hh"
#include "sample/HtsStreamingSampleAnalysis.hh"

namespace spd = spdlog;

using namespace ehunter;

// Filter out hom-ref and/or missing-genotype loci (per --skip-hom-ref / --skip-missing-genotypes) from
// regionCatalog and sampleFindings, using the same shouldFilterLocus predicate as the streaming modes.
static void filterSkippedLoci(
    RegionCatalog& regionCatalog, SampleFindings& sampleFindings, bool skipHomRef, bool skipMissing)
{
    RegionCatalog filteredRegionCatalog;
    SampleFindings filteredSampleFindings;

    const size_t locusCount = sampleFindings.size();
    for (size_t locusIndex = 0; locusIndex < locusCount; ++locusIndex)
    {
        const LocusSpecification& locusSpec = regionCatalog[locusIndex];
        const LocusFindings& locusFindings = sampleFindings[locusIndex];

        if (!shouldFilterLocus(locusSpec, locusFindings, skipHomRef, skipMissing))
        {
            filteredRegionCatalog.push_back(locusSpec);
            filteredSampleFindings.push_back(std::move(sampleFindings[locusIndex]));
        }
    }

    regionCatalog = std::move(filteredRegionCatalog);
    sampleFindings = std::move(filteredSampleFindings);
}

// Drop catalog loci that cannot be reliably genotyped from reads of the given length: those whose repeat
// reference region is wider than 2x the read length, or whose repeat motif is longer than half the read
// length. This mirrors the existing ">5 Ns in the flanks" catalog filter (see extendLocusStructure in
// io/LocusSpecDecoding.cpp): each dropped locus is logged and excluded from all output.
static void filterLociByReadLength(LocusDescriptionCatalog& locusDescriptionCatalog, int typicalReadLength)
{
    if (typicalReadLength <= 0)
    {
        return;
    }

    LocusDescriptionCatalog keptLoci;
    keptLoci.reserve(locusDescriptionCatalog.size());
    for (auto& locusDescription : locusDescriptionCatalog)
    {
        const int64_t referenceRegionWidth
            = locusDescription.locusWithoutFlanksEnd() - locusDescription.locusWithoutFlanksStart();
        if (referenceRegionWidth > 2 * static_cast<int64_t>(typicalReadLength))
        {
            spdlog::warn(
                "Skipping locus {}: reference region width ({} bp) is more than 2x the read length ({} bp)",
                locusDescription.locusId(), add_commas_at_thousands(static_cast<long>(referenceRegionWidth)),
                typicalReadLength);
            continue;
        }

        size_t largestMotifLength = 0;
        try
        {
            for (const auto& feature : decodeFeaturesFromRegex(locusDescription.locusStructure()))
            {
                if (feature.type == GraphBlueprintFeatureType::kSkippableRepeat
                    || feature.type == GraphBlueprintFeatureType::kUnskippableRepeat)
                {
                    largestMotifLength = std::max(largestMotifLength, feature.sequences.front().length());
                }
            }
        }
        catch (const std::exception&)
        {
            // Leave structurally-malformed loci for the normal decode path to report/skip rather than
            // dropping them here for an unrelated parse error.
            keptLoci.push_back(std::move(locusDescription));
            continue;
        }

        if (2 * static_cast<int64_t>(largestMotifLength) > typicalReadLength)
        {
            spdlog::warn(
                "Skipping locus {}: repeat motif length ({} bp) is more than half the read length ({} bp)",
                locusDescription.locusId(), add_commas_at_thousands(static_cast<unsigned long>(largestMotifLength)),
                typicalReadLength);
            continue;
        }

        keptLoci.push_back(std::move(locusDescription));
    }

    const size_t numSkipped = locusDescriptionCatalog.size() - keptLoci.size();
    if (numSkipped > 0)
    {
        spdlog::info("Skipped {} of {} loci whose motif or reference region is too large for the {} bp read length",
            add_commas_at_thousands(static_cast<unsigned long>(numSkipped)),
            add_commas_at_thousands(static_cast<unsigned long>(locusDescriptionCatalog.size())), typicalReadLength);
    }
    locusDescriptionCatalog = std::move(keptLoci);
}

template <typename T> static void writeToFile(std::string fileName, T streamable)
{
    std::ofstream outFile;
    boost::iostreams::filtering_ostream outStream;

    outFile.open(fileName, std::ios::out | std::ios::binary);
    if (!outFile.is_open())
    {
        throw std::runtime_error("Failed to open " + fileName + " for writing (" + strerror(errno) + ")");
    }

    if (fileName.size() > 2 && fileName.substr(fileName.size() - 2) == "gz")
    {
        outStream.push(boost::iostreams::gzip_compressor());
    }

    outStream.push(outFile);
    outStream << streamable;
    outStream.flush();
    outStream.reset();  // Ensures proper flushing of data

    outFile.close();
}

void setLogLevel(LogLevel logLevel)
{
    switch (logLevel)
    {
    case LogLevel::kTrace:
        spdlog::set_level(spdlog::level::trace);
        break;
    case LogLevel::kDebug:
        spdlog::set_level(spdlog::level::debug);
        break;
    case LogLevel::kInfo:
        spdlog::set_level(spdlog::level::info);
        break;
    case LogLevel::kWarn:
        spdlog::set_level(spdlog::level::warn);
        break;
    case LogLevel::kError:
        spdlog::set_level(spdlog::level::err);
        break;
    }
}

// Reconstructs the invoking command line (argv joined with spaces, unquoted) for the "CommandLine"
// RunInfo field, so the output JSON records how a run was invoked.
static std::string joinCommandLine(int argc, char** argv)
{
    std::string commandLine;
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
        {
            commandLine += " ";
        }
        commandLine += argv[i];
    }
    return commandLine;
}

int main(int argc, char** argv)
{
    spdlog::set_pattern("%Y-%m-%dT%H:%M:%S,[%v]");
    const std::time_t startedEpoch = currentEpochSeconds();
    const std::string commandLine = joinCommandLine(argc, argv);

    try
    {
        spdlog::info("Starting {}", kProgramVersion);

        auto optionalProgramParameters = tryLoadingProgramParameters(argc, argv);
        if (!optionalProgramParameters)
        {
            return 0;
        }
        const ProgramParameters& params = *optionalProgramParameters;

        srand(14345);  // For reproducibility

        setLogLevel(params.logLevel());

        const SampleParameters& sampleParams = params.sample();
        const InputPaths& inputPaths = params.inputPaths();

        spdlog::info("Analyzing sample {} from {}", sampleParams.id(), inputPaths.htsFile());
        spdlog::info("Initializing reference from {}", inputPaths.reference());
        FastaReference reference(inputPaths.reference(), extractReferenceContigInfo(inputPaths.htsFile()));

        spdlog::info("Reading catalog from {}", inputPaths.catalog());
        LocusDescriptionCatalog locusDescriptionCatalog = loadLocusDescriptions(params, reference);

        if (locusDescriptionCatalog.empty())
        {
            throw std::runtime_error(
                "No loci to analyze: the variant catalog is empty (check the catalog file and any "
                "--locus / --region / --start-with filters)");
        }

        // Skip loci that cannot be reliably genotyped at this sample's read length: those whose repeat
        // motif is larger than half the read length, or whose reference region is wider than 2x the read
        // length. Mirrors the >5-Ns-in-flanks catalog filter and applies to all analysis modes.
        filterLociByReadLength(locusDescriptionCatalog, probeTypicalReadLength(inputPaths));
        if (locusDescriptionCatalog.empty())
        {
            throw std::runtime_error(
                "No loci to analyze: every catalog locus was filtered out because its motif is larger than "
                "half the read length or its reference region is wider than 2x the read length");
        }

        const HeuristicParameters& heuristicParams = params.heuristics();
        const OutputPaths& outputPaths = params.outputPaths();
        if (params.analysisMode() == AnalysisMode::kLowMemStreaming
            || params.analysisMode() == AnalysisMode::kOptimizedStreaming)
        {
            spdlog::info("Running sample analysis in {} mode", analysisModeToString(params.analysisMode()));
            BamletWriterPtr bamletWriter = params.enableBamletOutput
                ? std::make_shared<BamletWriterImpl>(outputPaths.bamlet(), reference.contigInfo(), RegionCatalog{})
                : std::make_shared<BamletWriter>();
            htsLowMemStreamingSampleAnalysis(locusDescriptionCatalog, params, reference, bamletWriter, startedEpoch, commandLine);
            bamletWriter->finish();  // join the writer thread and surface any deferred bamlet write error
            return 0;
        }

        RegionCatalog regionCatalog = convertLocusDescriptionsToLocusSpecs(
            locusDescriptionCatalog, heuristicParams, reference);

        // Apply plot policy from CLI flags to all loci
        if (params.disableAllPlots())
        {
            for (auto& locusSpec : regionCatalog)
            {
                locusSpec.setPlotPolicy(PlotPolicy::kNone);
            }
        }
        else if (params.plotAll())
        {
            for (auto& locusSpec : regionCatalog)
            {
                locusSpec.setPlotPolicy(PlotPolicy::kAll);
            }
        }

        BamletWriterPtr bamletWriter;
        if (params.enableBamletOutput)
        {
            bamletWriter = std::make_shared<BamletWriterImpl>(outputPaths.bamlet(), reference.contigInfo(), regionCatalog);
        }
        else
        {
            bamletWriter = std::make_shared<BamletWriter>();
        }

        SampleFindings sampleFindings;
        if (params.analysisMode() == AnalysisMode::kSeeking)
        {
            spdlog::info("Running sample analysis in seeking mode");
            sampleFindings = htsSeekingSampleAnalysis(params, heuristicParams, regionCatalog, bamletWriter);
        }
        else
        {
            spdlog::info("Running sample analysis in streaming mode");
            sampleFindings = htsStreamingSampleAnalysis(params, heuristicParams, regionCatalog, bamletWriter);
        }
        bamletWriter->finish();  // join the writer thread and surface any deferred bamlet write error

        // Filter out hom-ref and/or missing-genotype loci if --skip-hom-ref / --skip-missing-genotypes is enabled
        if (params.skipHomRef() || params.skipMissingGenotypes())
        {
            filterSkippedLoci(regionCatalog, sampleFindings, params.skipHomRef(), params.skipMissingGenotypes());
        }

        spdlog::info("Writing output to disk");
        VcfWriter vcfWriter(sampleParams.id(), reference, regionCatalog, sampleFindings);
        writeToFile(outputPaths.vcf(), vcfWriter);

        JsonWriter jsonWriter(sampleParams, reference.contigInfo(), regionCatalog, sampleFindings, params.copyCatalogFields(), params.genotypeQualityModel().get(), startedEpoch, params.threadCount, params.analysisMode(), commandLine);
        writeToFile(outputPaths.json(), jsonWriter);
    }
    catch (const std::exception& e)
    {
        spdlog::error(e.what());
        return 1;
    }

    return 0;
}
