//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/IterativeJsonWriter.hh"

#include <cctype>
#include <cerrno>
#include <cstring>
#include <exception>
#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "app/Version.hh"
#include "core/Common.hh"
#include "core/ReadSupportCalculator.hh"
#include "genotype_quality/GenotypeQualityModel.hh"

namespace ehunter
{

using std::map;
using std::string;
using Json = nlohmann::json;
using boost::optional;
using std::to_string;
using std::vector;

// Position just past `literal` if `jsonString` has it at `position` (after any spaces), else npos.
static size_t matchAfterSpaces(const std::string& jsonString, size_t position, const char* literal)
{
    if (position >= jsonString.size())
    {
        return std::string::npos; // includes npos from a failed earlier match
    }
    while (position < jsonString.size() && jsonString[position] == ' ')
    {
        ++position;
    }
    const size_t length = std::strlen(literal);
    return jsonString.compare(position, length, literal) == 0 ? position + length : std::string::npos;
}

// Position just past the digits starting at `position`, or npos if there are none.
static size_t matchDigits(const std::string& jsonString, size_t position)
{
    const size_t start = position;
    while (position < jsonString.size() && std::isdigit(static_cast<unsigned char>(jsonString[position])))
    {
        ++position;
    }
    return position > start ? position : std::string::npos;
}

// Puts each MotifComposition count object, which dump(2) spreads over four lines, on one line:
// { "count": 143, "reads": 26 }. Nothing else in a locus record has this shape. A plain scan rather than a
// std::regex, which costs several percent of the run time when applied to every locus record.
static std::string inlineMotifCountObjects(const std::string& jsonString)
{
    if (jsonString.find("\"count\": ") == std::string::npos)
    {
        return jsonString;
    }
    std::string result;
    result.reserve(jsonString.size());
    size_t index = 0;
    while (index < jsonString.size())
    {
        if (jsonString.compare(index, 2, "{\n") == 0)
        {
            const size_t countStart = matchAfterSpaces(jsonString, index + 2, "\"count\": ");
            const size_t countEnd = countStart == std::string::npos ? countStart : matchDigits(jsonString, countStart);
            const size_t readsStart = countEnd == std::string::npos
                ? countEnd
                : matchAfterSpaces(jsonString, matchAfterSpaces(jsonString, countEnd, ",\n"), "\"reads\": ");
            const size_t readsEnd = readsStart == std::string::npos ? readsStart : matchDigits(jsonString, readsStart);
            const size_t objectEnd = readsEnd == std::string::npos
                ? readsEnd
                : matchAfterSpaces(jsonString, matchAfterSpaces(jsonString, readsEnd, "\n"), "}");
            if (objectEnd != std::string::npos)
            {
                result += "{ \"count\": ";
                result.append(jsonString, countStart, countEnd - countStart);
                result += ", \"reads\": ";
                result.append(jsonString, readsStart, readsEnd - readsStart);
                result += " }";
                index = objectEnd;
                continue;
            }
        }
        result += jsonString[index++];
    }
    return result;
}

std::string jsonDocumentHeader(const SampleParameters& sampleParams)
{
    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams.sex());

    const std::string jsonString = std::regex_replace(sampleParametersRecord.dump(2), std::regex("\n"), "\n  ");
    return "{\n  \"SampleParameters\": " + jsonString + ",\n  \"LocusResults\": {";
}

IterativeJsonWriter::IterativeJsonWriter(
    const SampleParameters& sampleParams,
    const ReferenceContigInfo& contigInfo,
    const std::string& outputFilePath,
    bool copyCatalogFields,
    const gq::GenotypeQualityModel* qualityModel,
    std::time_t startedEpoch,
    int threadCount,
    AnalysisMode analysisMode,
    const std::string& commandLine,
    JsonOutputMode outputMode,
    bool hasExistingRecords,
    OptimizedStreamingGenotypingApproach genotypingApproach)
    : contigInfo_(contigInfo)
    , outputFilePath_(outputFilePath)
    , firstRecord_(outputMode == JsonOutputMode::kTruncate || !hasExistingRecords)
    , copyCatalogFields_(copyCatalogFields)
    , qualityModel_(qualityModel)
    , startedEpoch_(startedEpoch)
    , threadCount_(threadCount)
    , analysisMode_(analysisMode)
    , genotypingApproach_(genotypingApproach)
    , commandLine_(commandLine)
{
    const bool append = outputMode == JsonOutputMode::kAppendAfterHeader;
    const std::ios::openmode openMode
        = std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc);
	outFile_.open(outputFilePath, openMode);
	if (!outFile_)
	{
		throw std::runtime_error("Failed to open file: " + outputFilePath);
	}

	// Check if we need to compress with Gzip
	if (outputFilePath.size() > 2 && outputFilePath.substr(outputFilePath.size() - 2) == "gz")
	{
		outStream_.push(boost::iostreams::gzip_compressor());
	}

	outStream_.push(outFile_);

    // In append mode the header (and any already-written records) are whatever the interrupted run left
    // behind, so only a fresh document writes one.
    if (!append)
    {
        outStream_ << jsonDocumentHeader(sampleParams);
    }
}


void IterativeJsonWriter::addRecord(const LocusSpecification& locusSpec, const LocusFindings& locusFindings) {
	const std::string& locusId(locusSpec.locusId());

    Json locusRecord;

    // Copy extra annotation fields from input catalog first (if enabled),
    // so computed values take precedence in case of field name collisions
    if (copyCatalogFields_ && locusSpec.extraFields().has_value())
    {
        const nlohmann::json& extraFields = locusSpec.extraFields().value();
        for (auto it = extraFields.begin(); it != extraFields.end(); ++it)
        {
            locusRecord[it.key()] = it.value();
        }
    }

    locusRecord["LocusId"] = locusId;
    // See JsonWriter::write -- the emitted field and the model's `coverage` feature share one value.
    const double locusCoverage = std::round(locusFindings.stats.depth() * 100) / 100.0;
    locusRecord["Coverage"] = locusCoverage;
    locusRecord["ReadLength"] = locusFindings.stats.meanReadLength();
    locusRecord["FragmentLength"] = locusFindings.stats.meanFragLength();
    locusRecord["AlleleCount"] = static_cast<int>(locusFindings.stats.alleleCount());

    Json variantRecords;
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

        VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec, qualityModel_, locusCoverage);
        variantIdAndFindings.second->accept(&variantWriter);
        variantRecords[variantId] = variantWriter.record();
    }

    if (!variantRecords.empty())
    {
        locusRecord["Variants"] = variantRecords;
    }
    if (locusFindings.repeatAllelePhasing)
    {
        locusRecord["RepeatAllelePhasing"] = encodeRepeatAllelePhasing(*locusFindings.repeatAllelePhasing);
    }

    std::string jsonString
        = std::regex_replace(inlineMotifCountObjects(locusRecord.dump(2)), std::regex("\n"), "\n    ");
    if (!firstRecord_)
        outStream_ << ", ";
    outStream_ << "\n    \"" << locusId << "\": " << jsonString;

    firstRecord_ = false;
}

void IterativeJsonWriter::addSkippedRecord(const std::string& locusId, const std::string& reason) {
    Json locusRecord;
    locusRecord["LocusId"] = locusId;
    locusRecord["Status"] = "skipped";
    locusRecord["Reason"] = reason;

    std::string jsonString = std::regex_replace(locusRecord.dump(2), std::regex("\n"), "\n    ");
    if (!firstRecord_)
        outStream_ << ", ";
    outStream_ << "\n    \"" << locusId << "\": " << jsonString;

    firstRecord_ = false;
}

std::uintmax_t IterativeJsonWriter::flushAndGetFileSize()
{
    outStream_.flush();
    outFile_.flush();
    if (!outStream_ || !outFile_)
    {
        throw std::runtime_error("Failed to write " + outputFilePath_ + " (" + std::strerror(errno) + ")");
    }
    return boost::filesystem::file_size(outputFilePath_);
}

void IterativeJsonWriter::close()
{
    if (closed_)
    {
        return;
    }
    closed_ = true;

    // RunInfo is written here, after all records, because "Completed" is only known once
    // genotyping has actually finished (the SampleParameters header above is written up front,
    // before any genotyping happens).
    //
    // close() also runs from the destructor so an exception unwinding past this writer still leaves
    // a well-formed, parseable JSON file (see the destructor comment) -- but that means this can run
    // mid-crash. std::uncaught_exceptions() tells us which case we're in: only report "Completed"/
    // "Runtime" when they reflect a real, successful finish, so a crashed run doesn't look complete.
    const bool completedNormally = std::uncaught_exceptions() == 0;

    Json runInfoRecord;
    runInfoRecord["Source"] = kSourceUrl;
    runInfoRecord["Version"] = kCommitSha;
    runInfoRecord["AnalysisMode"] = analysisModeToString(analysisMode_);
    if (analysisMode_ == AnalysisMode::kOptimizedStreaming)
    {
        // --genotyping-approach changes what optimized-streaming does at every locus (only-full writes no
        // QuickGenotype fields), so the approach is recorded next to the mode, as buildRunInfoJson does for
        // the multi-threaded run.
        runInfoRecord["GenotypingApproach"] = optimizedStreamingGenotypingApproachToString(genotypingApproach_);
    }
    runInfoRecord["Threads"] = threadCount_;
    runInfoRecord["Started"] = formatLocalTimestamp(startedEpoch_);
    if (completedNormally)
    {
        const std::time_t completedEpoch = currentEpochSeconds();
        runInfoRecord["Completed"] = formatLocalTimestamp(completedEpoch);
        runInfoRecord["Runtime"] = formatRuntime(completedEpoch - startedEpoch_);
    }
    runInfoRecord["PeakRssMemoryMb"] = peakRssMemoryMB();
    runInfoRecord["CommandLine"] = commandLine_;
    if (qualityModel_)
    {
        runInfoRecord["GenotypeQualityModelVersion"] = qualityModel_->version;
    }
    std::string runInfoJson = std::regex_replace(runInfoRecord.dump(2), std::regex("\n"), "\n  ");

    // Close the "LocusResults" object, then append RunInfo and close the outer scope JSON object.
    outStream_ << "\n  },\n  \"RunInfo\": " << runInfoJson << "\n}\n";
    outStream_.flush();
    outStream_.reset();
    outFile_.flush();
    const bool failed = !outFile_;
    if (outFile_.is_open())
    {
        outFile_.close();
    }
    // A write that silently failed (a full disk, most likely) would leave a truncated document that the
    // end-of-run merge would accept as complete, quietly dropping the loci it could not write.
    if (failed || !outFile_)
    {
        throw std::runtime_error("Failed to write " + outputFilePath_ + " (" + std::strerror(errno) + ")");
    }
}

IterativeJsonWriter::~IterativeJsonWriter()
{
    try
    {
        close();
    }
    catch (...)
    {
        // Never throw from a destructor.
    }
}


}
