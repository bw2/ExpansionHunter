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

#include <iomanip>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "core/Common.hh"
#include "core/ReadSupportCalculator.hh"

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
    bool copyCatalogFields)
    : contigInfo_(contigInfo), firstRecord_(true), copyCatalogFields_(copyCatalogFields)
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

        VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
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

void IterativeJsonWriter::close()
{
    outStream_ << "\n  }\n}\n"; //close the "LocusRecords" and outer scope JSON objects.
    outStream_.flush();
    outStream_.reset();
    if (outFile_.is_open())
    {
        outFile_.close();
    }
}


}
