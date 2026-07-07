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

#include <exception>
#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string/join.hpp>
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

IterativeJsonWriter::IterativeJsonWriter(
    const SampleParameters& sampleParams,
    const ReferenceContigInfo& contigInfo,
    const std::string& outputFilePath,
    bool copyCatalogFields,
    const gq::GenotypeQualityModel* qualityModel,
    std::time_t startedEpoch,
    int threadCount,
    AnalysisMode analysisMode,
    const std::string& commandLine)
    : contigInfo_(contigInfo)
    , firstRecord_(true)
    , copyCatalogFields_(copyCatalogFields)
    , qualityModel_(qualityModel)
    , startedEpoch_(startedEpoch)
    , threadCount_(threadCount)
    , analysisMode_(analysisMode)
    , commandLine_(commandLine)
{
	outFile_.open(outputFilePath, std::ios::out | std::ios::binary);
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

    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams.sex());

	std::string jsonString = std::regex_replace(sampleParametersRecord.dump(2), std::regex("\n"), "\n  ");
    outStream_ << "{\n";
    outStream_ << "  \"SampleParameters\": " << jsonString << ",\n";
    outStream_ << "  \"LocusResults\": {";
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
    locusRecord["Coverage"] = std::round(locusFindings.stats.depth() * 100) / 100.0;
    locusRecord["ReadLength"] = locusFindings.stats.meanReadLength();
    locusRecord["FragmentLength"] = locusFindings.stats.meanFragLength();
    locusRecord["AlleleCount"] = static_cast<int>(locusFindings.stats.alleleCount());

    Json variantRecords;
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

        VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec, qualityModel_);
        variantIdAndFindings.second->accept(&variantWriter);
        variantRecords[variantId] = variantWriter.record();
    }

    if (!variantRecords.empty())
    {
        locusRecord["Variants"] = variantRecords;
    }

	std::string jsonString = std::regex_replace(locusRecord.dump(2), std::regex("\n"), "\n    ");
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
    runInfoRecord["Threads"] = threadCount_;
    runInfoRecord["Started"] = formatLocalTimestamp(startedEpoch_);
    if (completedNormally)
    {
        const std::time_t completedEpoch = currentEpochSeconds();
        runInfoRecord["Completed"] = formatLocalTimestamp(completedEpoch);
        runInfoRecord["Runtime"] = formatRuntime(completedEpoch - startedEpoch_);
    }
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
    if (outFile_.is_open())
    {
        outFile_.close();
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
