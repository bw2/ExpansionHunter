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

#pragma once

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <cstdint>
#include <fstream>

#include "core/Parameters.hh"
#include "io/JsonWriter.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

using Json = nlohmann::json;

// The exact bytes an IterativeJsonWriter emits before its first record: the opening brace, the
// SampleParameters object, and the start of the LocusResults object. --resume uses its size to tell
// whether a temp file it cut back still holds any record (see io/ResumeCheckpoint.hh).
std::string jsonDocumentHeader(const SampleParameters& sampleParams);

// How a writer should attach to its output file.
enum class JsonOutputMode
{
    kTruncate,          // create/overwrite the file and write the document header
    kAppendAfterHeader, // open an existing partial document and continue after its last record
};

class IterativeJsonWriter
{
public:
	IterativeJsonWriter(const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo,
		const std::string& outputFilePath, bool copyCatalogFields = false,
		const gq::GenotypeQualityModel* qualityModel = nullptr, std::time_t startedEpoch = 0,
		int threadCount = 1, AnalysisMode analysisMode = AnalysisMode::kSeeking,
		const std::string& commandLine = "", JsonOutputMode outputMode = JsonOutputMode::kTruncate,
		bool hasExistingRecords = false);
	// Ensure the JSON document is closed even when an exception unwinds past the writer; otherwise
	// the output file is left missing its trailing `}}` braces and is unparseable.
	~IterativeJsonWriter();

	void addRecord(const LocusSpecification& locusSpec,  const LocusFindings& locusFindings);
    void addSkippedRecord(const std::string& locusId, const std::string& reason);
    // Flush everything written so far through to the output file and return the file's size in bytes.
    // --resume records this after each locus so an interrupted run's temp file can be cut back to its last
    // finished locus. Only meaningful for an uncompressed file, which the --resume temp files always are.
    std::uintmax_t flushAndGetFileSize();
    // Close the output file (idempotent). Throws if any of the writing failed, so a truncated output is
    // never mistaken for a complete one; the destructor swallows that, an explicit call propagates it.
    void close();

private:
    const ReferenceContigInfo& contigInfo_;
    std::string outputFilePath_;
    std::ofstream outFile_;
    boost::iostreams::filtering_ostream outStream_;
    bool firstRecord_;
    bool copyCatalogFields_;
    const gq::GenotypeQualityModel* qualityModel_;
    std::time_t startedEpoch_;
    int threadCount_;
    AnalysisMode analysisMode_;
    std::string commandLine_;
    bool closed_ = false;
};

}