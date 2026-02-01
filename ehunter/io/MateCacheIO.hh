//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
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

#include <string>

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"
#include "sample/MateExtractor.hh"

namespace ehunter
{

namespace htshelpers
{

class MateCacheIO
{
public:
    // Write mate cache to BAM file
    static void writeCache(
        const std::string& outputPath,
        const ReferenceContigInfo& contigInfo,
        const MateCache& cache,
        size_t catalogSize);

    // Read mate cache from BAM file, returns true on success
    static bool readCache(
        const std::string& inputPath,
        const ReferenceContigInfo& contigInfo,
        MateCache& cache,
        size_t expectedCatalogSize);
};

}

}
