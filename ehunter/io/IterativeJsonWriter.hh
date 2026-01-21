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
#include <fstream>

#include "core/Parameters.hh"
#include "io/JsonWriter.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

using Json = nlohmann::json;

class IterativeJsonWriter
{
public:
	IterativeJsonWriter(const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo,
		const std::string& outputFilePath, bool copyCatalogFields = false);

	void addRecord(const LocusSpecification& locusSpec,  const LocusFindings& locusFindings);
    void close();  // Close the output file

private:
    const ReferenceContigInfo& contigInfo_;
    std::ofstream outFile_;
    boost::iostreams::filtering_ostream outStream_;
    bool firstRecord_;
    bool copyCatalogFields_;
};

}