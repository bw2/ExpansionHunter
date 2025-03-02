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

#include "core/HtsHelpers.hh"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>

#include "spdlog/spdlog.h"

using std::pair;
using std::string;
using std::vector;

namespace spd = spdlog;

namespace ehunter
{
namespace htshelpers
{

/* Open a local or remote .bai or .crai index file. This supports direct access to files in Google Storage or S3.
 *
 * @param htsFilePtr pointer to an already open .bam or .cram file
 * @param htsFilePath path of the .bam or .cram file
 */
hts_idx_t* openHtsIndex(htsFile* htsFilePtr, std::string htsFilePath) {
    std::string indexFilePath = htsFilePath + (htsFilePtr->format.format == cram ? ".crai" : ".bai");
    hts_idx_t* htsIndexPtr_ = sam_index_load3(htsFilePtr, htsFilePath.c_str(), indexFilePath.c_str(), HTS_IDX_SAVE_REMOTE);

    /*
    // if file starts with gs://
    if ( this->file_name.find("gs://") == 0 ) {
        // run gcloud auth command to refresh the GCS_OAUTH_TOKEN environment variable with up to 5 retries
        for(int32_t i=0; i < 5; ++i) {
            std::string oauth_token = exec_cmd("gcloud auth application-default print-access-token");
            //notice("Debug message: refreshed OAUTH TOKEN to %s", oauth_token.c_str());
            if ( setenv("GCS_OAUTH_TOKEN",oauth_token.c_str(),1) != 0 ) {
                error("Failed to update the environment variable GCS_OAUTH_TOKEN to %s", oauth_token.c_str());
            }
            notice("Successfully refreshed OAUTH TOKEN");
            file = hts_open(this->file_name.c_str(), "r");
            if ( file ) break;
        }
    }
    */

    if (!htsIndexPtr_)
    {
        //try the alternative index filename, replacing the .cram or .bam suffix with .crai or .bai
        if (htsFilePtr->format.format == cram) {
            indexFilePath = htsFilePath.substr(0, htsFilePath.length() - 4) + ".crai";
        } else {
            indexFilePath = htsFilePath.substr(0, htsFilePath.length() - 3) + ".bai";
        }
        htsIndexPtr_ = sam_index_load3(htsFilePtr, htsFilePath.c_str(), indexFilePath.c_str(), HTS_IDX_SAVE_REMOTE);
    }

    return htsIndexPtr_;
}

string decodeQuals(bam1_t* htsAlignPtr)
{
    string quals;
    uint8_t* htsQualsPtr = bam_get_qual(htsAlignPtr);
    const int readLength = htsAlignPtr->core.l_qseq;
    quals.resize(readLength);

    for (int index = 0; index < readLength; ++index)
    {
        quals[index] = static_cast<char>(33 + htsQualsPtr[index]);
    }

    return quals;
}

string decodeBases(bam1_t* htsAlignPtr)
{
    string bases;
    uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
    const int32_t readLength = htsAlignPtr->core.l_qseq;
    bases.resize(readLength);

    for (int32_t index = 0; index < readLength; ++index)
    {
        bases[index] = seq_nt16_str[bam_seqi(htsSeqPtr, index)];
    }

    return bases;
}

LinearAlignmentStats decodeAlignmentStats(bam1_t* htsAlignPtr)
{
    LinearAlignmentStats alignmentStats;
    alignmentStats.chromId = htsAlignPtr->core.tid;
    alignmentStats.pos = htsAlignPtr->core.pos; //0-based
    alignmentStats.mapq = htsAlignPtr->core.qual;
    alignmentStats.mateChromId = htsAlignPtr->core.mtid;
    alignmentStats.matePos = htsAlignPtr->core.mpos;

    uint32_t samFlag = htsAlignPtr->core.flag;
    alignmentStats.isPaired = samFlag & BAM_FPAIRED;
    alignmentStats.isMapped = !(samFlag & BAM_FUNMAP);
    alignmentStats.isMateMapped = !(samFlag & BAM_FMUNMAP);
    alignmentStats.isSecondaryAlignment = samFlag & BAM_FSECONDARY;
    alignmentStats.isSupplementaryAlignment = samFlag & BAM_FSUPPLEMENTARY;

    if (alignmentStats.isMapped && htsAlignPtr->core.n_cigar > 0) {
        const uint32_t* cigarPtr = bam_get_cigar(htsAlignPtr);
        if (cigarPtr) {
            alignmentStats.cigar.assign(cigarPtr, cigarPtr + htsAlignPtr->core.n_cigar);
        }
    }

    return alignmentStats;
}

bool isPrimaryAlignment(bam1_t* htsAlignPtr)
{
    return !((htsAlignPtr->core.flag & BAM_FSECONDARY) || (htsAlignPtr->core.flag & BAM_FSUPPLEMENTARY));
}

/// Lower-case version of htslib seq_nt16_str table;
class Cache_seq_nt16_str_lc
{
public:
    Cache_seq_nt16_str_lc()
    {
        for (unsigned i(0); i < 16; ++i)
        {
            data[i] = std::tolower(seq_nt16_str[i]);
        }
    }

    char data[16];
};

static Cache_seq_nt16_str_lc seq_nt16_str_lc;

ReadId decodeReadId(bam1_t* htsAlignPtr) {
    const uint32_t samFlag = htsAlignPtr->core.flag;
    const bool isFirstMate = samFlag & BAM_FREAD1;

    const char* qname(bam_get_qname(htsAlignPtr));
    MateNumber mateNumber = isFirstMate ? MateNumber::kFirstMate : MateNumber::kSecondMate;
    ReadId readId(qname, mateNumber);

    return readId;
}

Read decodeRead(bam1_t* htsAlignPtr)
{
    ReadId readId = decodeReadId(htsAlignPtr);

    const uint32_t samFlag = htsAlignPtr->core.flag;
    const bool isReversed = samFlag & BAM_FREVERSE;

    string bases;

    {
        // Decode bases and convert low-quality bases to lowercase:
        static const uint8_t lowBaseQualityCutoff(20);

        const uint8_t* htsSeqPtr = bam_get_seq(htsAlignPtr);
        const uint8_t* htsQualPtr = bam_get_qual(htsAlignPtr);
        const int32_t readLength = htsAlignPtr->core.l_qseq;
        bases.resize(readLength);

        for (int32_t index = 0; index < readLength; ++index)
        {
            const char* converter((htsQualPtr[index] <= lowBaseQualityCutoff) ? seq_nt16_str_lc.data : seq_nt16_str);
            bases[index] = converter[bam_seqi(htsSeqPtr, index)];
        }
    }

    return { std::move(readId), std::move(bases), isReversed };
}

ReferenceContigInfo decodeContigInfo(bam_hdr_t* htsHeaderPtr)
{
    vector<pair<string, int64_t>> contigNamesAndSizes;
    const int32_t numContigs = htsHeaderPtr->n_targets;
    contigNamesAndSizes.reserve(numContigs);

    for (int32_t contigIndex = 0; contigIndex != numContigs; ++contigIndex)
    {
        const string contig = htsHeaderPtr->target_name[contigIndex];
        int64_t size = htsHeaderPtr->target_len[contigIndex];
        contigNamesAndSizes.push_back(std::make_pair(contig, size));
    }

    return ReferenceContigInfo(contigNamesAndSizes);
}

} // namespace htshelpers
}
