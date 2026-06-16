//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/SampleStats.hh"

#include <cstdint>
#include <string>
#include <utility>

#include <boost/filesystem.hpp>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "core/HtsHelpers.hh"

using std::pair;
using std::string;
using std::vector;

namespace ehunter
{

ReferenceContigInfo extractReferenceContigInfo(const std::string& htsFilePath)
{
    std::unique_ptr<samFile, decltype(&hts_close)> htsFilePtr(hts_open(htsFilePath.c_str(), "r"), hts_close);
    if (!htsFilePtr)
    {
        throw std::runtime_error("Failed to read " + htsFilePath);
    }

    std::unique_ptr<bam_hdr_t, decltype(&bam_hdr_destroy)> htsHeaderPtr(
        sam_hdr_read(htsFilePtr.get()), bam_hdr_destroy);
    if (!htsHeaderPtr)
    {
        throw std::runtime_error("Failed to read the header of " + htsFilePath);
    }

    return htshelpers::decodeContigInfo(htsHeaderPtr.get());
}

}
