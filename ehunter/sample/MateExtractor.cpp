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

#include "sample/MateExtractor.hh"

#include <stdexcept>

#include "core/HtsHelpers.hh"
#include "spdlog/spdlog.h"
#include "htslib/hts.h"

namespace ehunter
{

using boost::optional;
using std::string;

namespace htshelpers
{


MateExtractor::MateExtractor(const string& htsFilePath, const std::string& htsReferencePath, const bool cacheMates, const int farAwayMateDistanceThreshold)
    : htsFilePath_(htsFilePath)
    , htsReferencePath_(htsReferencePath)
    , contigInfo_({})
    , cacheMates_(cacheMates)
    , farAwayMateDistanceThreshold_(farAwayMateDistanceThreshold)
{
    openFile();
    loadHeader();
    loadIndex();
    htsAlignmentPtr_ = bam_init1();
}

MateExtractor::~MateExtractor()
{
    bam_destroy1(htsAlignmentPtr_);
    htsAlignmentPtr_ = nullptr;

    hts_idx_destroy(htsIndexPtr_);
    htsIndexPtr_ = nullptr;

    bam_hdr_destroy(htsHeaderPtr_);
    htsHeaderPtr_ = nullptr;

    sam_close(htsFilePtr_);
    htsFilePtr_ = nullptr;
}

void MateExtractor::openFile()
{
    htsFilePtr_ = sam_open(htsFilePath_.c_str(), "r");

    if (!htsFilePtr_)
    {
        throw std::runtime_error("Failed to read BAM file " + htsFilePath_);
    }

    // Required step for parsing of some CRAMs
    if (hts_set_fai_filename(htsFilePtr_, htsReferencePath_.c_str()) != 0)
    {
        throw std::runtime_error("Failed to set index of: " + htsReferencePath_);
    }
}

void MateExtractor::loadHeader()
{
    htsHeaderPtr_ = sam_hdr_read(htsFilePtr_);

    if (!htsHeaderPtr_)
    {
        throw std::runtime_error("Failed to read header of " + htsFilePath_);
    }

    contigInfo_ = htshelpers::decodeContigInfo(htsHeaderPtr_);
}

void MateExtractor::loadIndex()
{
    htsIndexPtr_ = openHtsIndex(htsFilePtr_, htsFilePath_);

    if (!htsIndexPtr_)
    {
        throw std::runtime_error("Failed to read index of " + htsFilePath_);
    }
}


void MateExtractor::addMateToCache(const ReadId& mateReadId, const FullRead& mate) {
    if (mateCache_.count(mateReadId) > 0) {
        throw std::logic_error("Mate cache already contains mate with id " + mateReadId.fragmentId());
    }
    mateCache_.emplace(mateReadId, std::move(mate));
}

boost::optional<FullRead> MateExtractor::extractMate(const ReadId& mateReadId, const GenomicRegion& genomicRegion) {
    extractedTotalCounter_ += 1;
    if (cacheMates_ && mateCache_.count(mateReadId) > 0) {
        return mateCache_.at(mateReadId);
    }

    hts_itr_t* htsRegionPtr_ = sam_itr_queryi(htsIndexPtr_, genomicRegion.contigIndex(), genomicRegion.start(), genomicRegion.end());

    if (!htsRegionPtr_)
    {
        const string& contigName = contigInfo_.getContigName(genomicRegion.contigIndex());
        const string regionEncoding
            = contigName + ":" + std::to_string(genomicRegion.start()) + "-" + std::to_string(genomicRegion.end());

        throw std::logic_error("Unable to jump to " + regionEncoding + " to recover a mate");
    }

    boost::optional<FullRead> mate = boost::none;
    while (sam_itr_next(htsFilePtr_, htsRegionPtr_, htsAlignmentPtr_) >= 0)
    {
        const bool isSecondaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSECONDARY;
        const bool isSupplementaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSUPPLEMENTARY;
        const bool isPrimaryAlignment = !(isSecondaryAlignment || isSupplementaryAlignment);
        if (!isPrimaryAlignment)
        {
            continue;
        }

        ReadId putativeMateReadId = decodeReadId(htsAlignmentPtr_);

        // if this is one of the fragment ids we're looking for, check if forms a proper pair with the 1st read with
        // this fragmentId. If it does, decode this read and add it to the list of mates to return
        const bool foundRequestedMate = putativeMateReadId == mateReadId;
        if (foundRequestedMate || cacheMates_) {
            Read putativeMate = htshelpers::decodeRead(htsAlignmentPtr_);
            LinearAlignmentStats mateStats = decodeAlignmentStats(htsAlignmentPtr_);
            FullRead putativeFullMate = {std::move(putativeMate), std::move(mateStats)};

            if (cacheMates_) {
                // cache the mate if --cache-mates flag was specified and the mate mapped sufficiently far away from the read
                mateCache_.emplace(putativeFullMate.r.readId(), std::move(putativeFullMate));
            }
            if (foundRequestedMate) {
                extractedFromDiskCounter_ += 1;
                mate = putativeFullMate;
                break;
            }
        }
    }

    hts_itr_destroy(htsRegionPtr_);

    if (mate == boost::none) {
        spdlog::warn("Failed to recover mate: {}", mateReadId.fragmentId());
    }

    return mate;
}


std::vector<FullRead> MateExtractor::extractMates(const MateRegionToRecover& mateRegionToRecover)
{
    std::unordered_set<ReadId, boost::hash<ReadId>> mateReadIdsNotFoundYet = mateRegionToRecover.mateReadIds;
    std::vector<FullRead> extractedMates;

    if (cacheMates_) {
        for (ReadId mateReadId : mateRegionToRecover.mateReadIds) {
            if (mateCache_.count(mateReadId) > 0) {
                extractedMates.push_back(mateCache_.at(mateReadId));
                mateReadIdsNotFoundYet.erase(mateReadId);
            }
        }
        if (mateReadIdsNotFoundYet.size() == 0) {
            extractedTotalCounter_ += extractedMates.size();
            return extractedMates;
        }
    }

    const int32_t searchRegionContigIndex = mateRegionToRecover.genomicRegion.contigIndex();
    const int32_t searchRegionStart = mateRegionToRecover.genomicRegion.start();
    const int32_t searchRegionEnd = mateRegionToRecover.genomicRegion.end();

    hts_itr_t* htsRegionPtr_
        = sam_itr_queryi(htsIndexPtr_, searchRegionContigIndex, searchRegionStart, searchRegionEnd);

    if (!htsRegionPtr_)
    {
        const string& contigName = contigInfo_.getContigName(searchRegionContigIndex);
        const string regionEncoding
            = contigName + ":" + std::to_string(searchRegionStart) + "-" + std::to_string(searchRegionEnd);

        throw std::logic_error("Unable to jump to " + regionEncoding + " to recover a mate");
    }

    while (sam_itr_next(htsFilePtr_, htsRegionPtr_, htsAlignmentPtr_) >= 0)
    {
        const bool isSecondaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSECONDARY;
        const bool isSupplementaryAlignment = htsAlignmentPtr_->core.flag & BAM_FSUPPLEMENTARY;
        const bool isPrimaryAlignment = !(isSecondaryAlignment || isSupplementaryAlignment);
        if (!isPrimaryAlignment)
        {
            continue;
        }

        ReadId putativeMateReadId = decodeReadId(htsAlignmentPtr_);

        // if this is one of the fragment ids we're looking for, check if forms a proper pair with the 1st read with
        // this fragmentId. If it does, decode this read and add it to the list of mates to return
        bool foundRequestedMate = mateReadIdsNotFoundYet.count(putativeMateReadId) > 0;
        if (foundRequestedMate || cacheMates_) {
            LinearAlignmentStats mateStats = decodeAlignmentStats(htsAlignmentPtr_);

            // cache the mate if --cache-mates flag was specified and the mate mapped sufficiently far away from the
            // read, suggesting there's a chance that it signals a large expansion at this locus
            bool shouldCacheThisMate = cacheMates_ && (
                mateStats.chromId != mateStats.mateChromId || std::abs(mateStats.pos - mateStats.matePos) >= farAwayMateDistanceThreshold_);

            if (foundRequestedMate || shouldCacheThisMate) {
                Read putativeMate = htshelpers::decodeRead(htsAlignmentPtr_);
                FullRead putativeFullMate = {std::move(putativeMate), std::move(mateStats)};
                if (foundRequestedMate) {
                    extractedFromDiskCounter_ += 1;
                    extractedMates.push_back(putativeFullMate);
                    mateReadIdsNotFoundYet.erase(putativeFullMate.r.readId());
                }
                if (shouldCacheThisMate) {
                    mateCache_.emplace(putativeFullMate.r.readId(), std::move(putativeFullMate));
                }
                if (foundRequestedMate && mateReadIdsNotFoundYet.size() == 0 && !cacheMates_) {
                    break;  // found all requested mates
                }
            }

        }
    }
    hts_itr_destroy(htsRegionPtr_);

    if (mateReadIdsNotFoundYet.size() > 0) {
        std::ostringstream out;
        out <<  "Failed to recover " << mateReadIdsNotFoundYet.size() << " mate(s): ";
        for (const auto& readId: mateReadIdsNotFoundYet) {
            out << readId.fragmentId() << "  ";
        }
        spdlog::warn(out.str());
    }

    extractedTotalCounter_ += extractedMates.size();
    return extractedMates;
}

}

}
