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

#include "io/MateCacheIO.hh"

#include <cstring>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <string>

#include "spdlog/spdlog.h"

// cppcheck-suppress missingInclude
#include "htslib/bgzf.h"
// cppcheck-suppress missingInclude
#include "htslib/hts.h"
// cppcheck-suppress missingInclude
#include "htslib/sam.h"

using std::string;

namespace spd = spdlog;

namespace ehunter
{
namespace htshelpers
{

// Set a base in the bit-backed read sequence representation of a BAM record. (4 bits per base).
// Sets position i to base c in sequence s.
#define bam1_seq_seti(s, i, c) ((s)[(i) >> 1] = ((s)[(i) >> 1] & 0xf << (((i)&1) << 2)) | (c) << ((~(i)&1) << 2))

static const string kCatalogSizeComment = "@CO\tEH:CATALOG-LOCI-COUNT:";

void MateCacheIO::writeCache(
    const string& outputPath,
    const ReferenceContigInfo& contigInfo,
    const MateCache& cache,
    size_t catalogSize)
{
    // Open BAM file for writing
    htsFile* filePtr = hts_open(outputPath.c_str(), "wb");
    if (!filePtr)
    {
        throw std::runtime_error("Failed to open BAM file for writing: " + outputPath);
    }

    // Create and initialize header
    bam_hdr_t* headerPtr = bam_hdr_init();
    if (!headerPtr)
    {
        hts_close(filePtr);
        throw std::runtime_error("Failed to initialize BAM header");
    }

    // Build header text with catalog size comment
    string headerText = "@HD\tVN:1.4\tSO:unknown\n";
    headerText += kCatalogSizeComment + std::to_string(catalogSize) + "\n";

    headerPtr->l_text = headerText.length();
    headerPtr->text = strdup(headerText.c_str());
    headerPtr->n_targets = contigInfo.numContigs();

    // Allocate memory for contig names and lengths
    headerPtr->target_len = (uint32_t*)calloc(contigInfo.numContigs(), sizeof(uint32_t));
    headerPtr->target_name = (char**)calloc(contigInfo.numContigs(), sizeof(char*));

    if (!headerPtr->target_len || !headerPtr->target_name)
    {
        bam_hdr_destroy(headerPtr);
        hts_close(filePtr);
        throw std::runtime_error("Failed to allocate memory for BAM header");
    }

    for (int32_t index = 0; index < contigInfo.numContigs(); ++index)
    {
        const string& contigName = contigInfo.getContigName(index);
        headerPtr->target_name[index] = (char*)malloc(contigName.length() + 1);
        if (!headerPtr->target_name[index])
        {
            bam_hdr_destroy(headerPtr);
            hts_close(filePtr);
            throw std::runtime_error("Failed to allocate memory for contig name");
        }
        memcpy(headerPtr->target_name[index], contigName.c_str(), contigName.length() + 1);
        headerPtr->target_len[index] = contigInfo.getContigSize(index);
    }

    // Write header
    if (bam_hdr_write(filePtr->fp.bgzf, headerPtr) != 0)
    {
        bam_hdr_destroy(headerPtr);
        hts_close(filePtr);
        throw std::runtime_error("Failed to write BAM header");
    }

    // Write each read in the cache
    for (const auto& entry : cache)
    {
        const FullRead& fullRead = entry.second;
        const Read& r = fullRead.r;
        const LinearAlignmentStats& s = fullRead.s;

        bam1_t* alignmentPtr = bam_init1();
        if (!alignmentPtr)
        {
            bam_hdr_destroy(headerPtr);
            hts_close(filePtr);
            throw std::runtime_error("Failed to initialize BAM record");
        }

        // Set core fields
        alignmentPtr->core.tid = s.chromId;
        alignmentPtr->core.pos = s.pos;
        alignmentPtr->core.qual = s.mapq;
        alignmentPtr->core.mtid = s.mateChromId;
        alignmentPtr->core.mpos = s.matePos;

        // Set flags
        uint16_t flag = 0;
        if (r.mateNumber() == MateNumber::kFirstMate)
        {
            flag |= BAM_FREAD1;
        }
        else
        {
            flag |= BAM_FREAD2;
        }
        if (r.isReversed())
        {
            flag |= BAM_FREVERSE;
        }
        if (s.isPaired)
        {
            flag |= BAM_FPAIRED;
        }
        if (!s.isMapped)
        {
            flag |= BAM_FUNMAP;
        }
        if (!s.isMateMapped)
        {
            flag |= BAM_FMUNMAP;
        }
        alignmentPtr->core.flag = flag;

        // Set qname
        const string& qname = r.fragmentId();
        alignmentPtr->core.l_qname = qname.length() + 1;  // +1 for null terminator

        // Set sequence length
        const string& sequence = r.sequence();
        alignmentPtr->core.l_qseq = sequence.length();

        // Set CIGAR
        alignmentPtr->core.n_cigar = s.cigar.size();

        // Calculate data length: qname + cigar + seq + qual
        int cigarLen = s.cigar.size() * sizeof(uint32_t);
        int seqLen = (sequence.length() + 1) / 2;  // 4 bits per base
        int qualLen = sequence.length();

        alignmentPtr->l_data = alignmentPtr->core.l_qname + cigarLen + seqLen + qualLen;
        alignmentPtr->m_data = alignmentPtr->l_data;
        kroundup32(alignmentPtr->m_data);
        alignmentPtr->data = (uint8_t*)realloc(alignmentPtr->data, alignmentPtr->m_data);

        if (!alignmentPtr->data)
        {
            bam_destroy1(alignmentPtr);
            bam_hdr_destroy(headerPtr);
            hts_close(filePtr);
            throw std::runtime_error("Failed to allocate memory for BAM record data");
        }

        // Copy qname
        memcpy(alignmentPtr->data, qname.c_str(), alignmentPtr->core.l_qname);

        // Copy CIGAR
        if (!s.cigar.empty())
        {
            uint32_t* cigarPtr = bam_get_cigar(alignmentPtr);
            memcpy(cigarPtr, s.cigar.data(), cigarLen);
        }

        // Encode sequence
        uint8_t* seqPtr = bam_get_seq(alignmentPtr);
        memset(seqPtr, 0, seqLen);  // Initialize to zero
        for (size_t i = 0; i < sequence.length(); ++i)
        {
            bam1_seq_seti(seqPtr, i, seq_nt16_table[(unsigned char)sequence[i]]);
        }

        // Set quality scores (use default high quality for uppercase, low for lowercase)
        uint8_t* qualPtr = bam_get_qual(alignmentPtr);
        const int kLowQualityScore = 0;
        const int kHighQualityScore = 40;
        for (size_t i = 0; i < sequence.length(); ++i)
        {
            qualPtr[i] = isupper(sequence[i]) ? kHighQualityScore : kLowQualityScore;
        }

        // Write record
        if (bam_write1(filePtr->fp.bgzf, alignmentPtr) < 0)
        {
            bam_destroy1(alignmentPtr);
            bam_hdr_destroy(headerPtr);
            hts_close(filePtr);
            throw std::runtime_error("Failed to write BAM record");
        }

        bam_destroy1(alignmentPtr);
    }

    bam_hdr_destroy(headerPtr);
    if (hts_close(filePtr) != 0)
    {
        throw std::runtime_error("Failed to close BAM file: " + outputPath);
    }
}

bool MateCacheIO::readCache(
    const string& inputPath,
    [[maybe_unused]] const ReferenceContigInfo& contigInfo,
    MateCache& cache,
    size_t expectedCatalogSize)
{
    // Check if file exists
    std::ifstream testFile(inputPath);
    if (!testFile.good())
    {
        return false;
    }
    testFile.close();

    // Open BAM file for reading
    htsFile* filePtr = sam_open(inputPath.c_str(), "r");
    if (!filePtr)
    {
        spd::warn("Failed to open mate cache BAM file: {}", inputPath);
        return false;
    }

    // Read header
    bam_hdr_t* headerPtr = sam_hdr_read(filePtr);
    if (!headerPtr)
    {
        spd::warn("Failed to read header from mate cache BAM file: {}", inputPath);
        hts_close(filePtr);
        return false;
    }

    // Parse catalog size from header comment
    size_t catalogSize = 0;
    bool foundCatalogSize = false;
    if (headerPtr->text)
    {
        string headerText(headerPtr->text, headerPtr->l_text);
        std::regex catalogRegex("@CO\\tEH:CATALOG-LOCI-COUNT:(\\d+)");
        std::smatch match;
        if (std::regex_search(headerText, match, catalogRegex))
        {
            catalogSize = std::stoull(match[1].str());
            foundCatalogSize = true;
        }
    }

    if (!foundCatalogSize)
    {
        spd::warn("Mate cache BAM file missing catalog size comment: {}", inputPath);
        bam_hdr_destroy(headerPtr);
        hts_close(filePtr);
        return false;
    }

    if (catalogSize != expectedCatalogSize)
    {
        spd::warn(
            "Mate cache catalog size mismatch: expected {}, found {} in {}",
            expectedCatalogSize, catalogSize, inputPath);
        bam_hdr_destroy(headerPtr);
        hts_close(filePtr);
        return false;
    }

    // Read records
    bam1_t* alignmentPtr = bam_init1();
    if (!alignmentPtr)
    {
        bam_hdr_destroy(headerPtr);
        hts_close(filePtr);
        return false;
    }

    while (sam_read1(filePtr, headerPtr, alignmentPtr) >= 0)
    {
        // Decode fragment ID (qname)
        const char* qname = bam_get_qname(alignmentPtr);
        string fragmentId(qname);

        // Decode mate number from flags
        uint32_t samFlag = alignmentPtr->core.flag;
        bool isFirstMate = samFlag & BAM_FREAD1;
        MateNumber mateNumber = isFirstMate ? MateNumber::kFirstMate : MateNumber::kSecondMate;

        // Decode isReversed
        bool isReversed = samFlag & BAM_FREVERSE;

        // Decode sequence
        string sequence;
        uint8_t* seqPtr = bam_get_seq(alignmentPtr);
        int32_t seqLength = alignmentPtr->core.l_qseq;
        sequence.resize(seqLength);

        // Use quality scores to determine case (like in HtsHelpers.cpp)
        static const uint8_t lowBaseQualityCutoff = 20;
        uint8_t* qualPtr = bam_get_qual(alignmentPtr);

        for (int32_t i = 0; i < seqLength; ++i)
        {
            char base = seq_nt16_str[bam_seqi(seqPtr, i)];
            if (qualPtr[i] <= lowBaseQualityCutoff)
            {
                base = std::tolower(base);
            }
            sequence[i] = base;
        }

        // Create ReadId and Read
        ReadId readId(fragmentId, mateNumber);
        Read read(readId, sequence, isReversed);

        // Decode alignment stats
        LinearAlignmentStats stats;
        stats.chromId = alignmentPtr->core.tid;
        stats.pos = alignmentPtr->core.pos;
        stats.mapq = alignmentPtr->core.qual;
        stats.mateChromId = alignmentPtr->core.mtid;
        stats.matePos = alignmentPtr->core.mpos;
        stats.isPaired = samFlag & BAM_FPAIRED;
        stats.isMapped = !(samFlag & BAM_FUNMAP);
        stats.isMateMapped = !(samFlag & BAM_FMUNMAP);
        stats.isSecondaryAlignment = samFlag & BAM_FSECONDARY;
        stats.isSupplementaryAlignment = samFlag & BAM_FSUPPLEMENTARY;

        // Decode CIGAR
        if (alignmentPtr->core.n_cigar > 0)
        {
            const uint32_t* cigarPtr = bam_get_cigar(alignmentPtr);
            stats.cigar.assign(cigarPtr, cigarPtr + alignmentPtr->core.n_cigar);
        }

        // Add to cache
        FullRead fullRead(read, stats);
        cache.emplace(readId, fullRead);
    }

    bam_destroy1(alignmentPtr);
    bam_hdr_destroy(headerPtr);
    hts_close(filePtr);

    return true;
}

}  // namespace htshelpers
}  // namespace ehunter
