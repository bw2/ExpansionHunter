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

#include <string>
#include <vector>

namespace ehunter
{

// Merge per-region temp VCF files (written by separate region workers, already position-sorted
// within each file) into one final VCF in genomic region order. The header of regionTempPaths[0]
// is kept; the duplicate headers of subsequent files are dropped. The final output is gzip-compressed
// iff finalPath ends in "gz"; region temp files are read as plain (uncompressed) text.
void mergeRegionVcfFiles(const std::string& finalPath, const std::vector<std::string>& regionTempPaths);

// Merge per-region temp JSON files (written by separate region workers) into one final JSON in
// genomic region order. The final output is gzip-compressed iff finalPath ends in "gz"; region temp
// files are read as plain (uncompressed) text. Each region file's own "RunInfo" record (each worker's
// local start/finish time) is discarded; runInfoJson (an already-formatted, 2-space-indented JSON
// object, run-wide and supplied by the caller) is appended instead.
void mergeRegionJsonFiles(
    const std::string& finalPath, const std::vector<std::string>& regionTempPaths, const std::string& runInfoJson);

}
