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

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
}

#include <stdexcept>

#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{

// Thrown when a read fails (see isReadFailure). Code that catches a locus's errors so the run can go on lets this one
// through, since it isn't a problem with the locus: the run stops instead of quietly writing results from fewer reads
class ReadFailureError : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

hts_idx_t* openHtsIndex(htsFile* htsFilePtr, std::string htsFilePath, const std::string& htsIndexPath);
LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr);
bool isPrimaryAlignment(bam1_t* htsAlignPtr);
// Says whether a sam_itr_next or sam_read1 return code means that a read failed (for example because the connection
// to a gs:// file dropped), rather than the normal end of the region or file. For a CRAM read with a thread pool it
// can miss a failed read of a container's header, so whole-file reads use isWholeFileReadFailure instead
bool isReadFailure(htsFile* htsFilePtr, int samReadReturnCode);
// Like isReadFailure, for sam_read1 reading a whole file, which a CRAM ends by reading its EOF container
bool isWholeFileReadFailure(htsFile* htsFilePtr, int samReadReturnCode);
ReadId decodeReadId(bam1_t* htsAlignPtr);
Read decodeRead(bam1_t* htsAlignPtr);
ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr);

} // namespace htshelpers

}
