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

#include <boost/optional.hpp>

extern "C"
{
#include "htslib/hts.h"
#include "htslib/sam.h"
}

#include "core/GenomicRegion.hh"
#include "core/Read.hh"
#include "core/ReferenceContigInfo.hh"

namespace ehunter
{

namespace htshelpers
{


// Represents one or more mates that need to be recovered from a specific genomic region.
struct MateRegionToRecover {
    GenomicRegion genomicRegion;

    //stores the ReadId of each mate that needs to be recovered from the genomicRegion
    std::unordered_set<ReadId, boost::hash<ReadId>> mateReadIds;
};

using MateCache = std::unordered_map<ReadId, FullRead, boost::hash<ReadId>>;

class MateExtractor
{
public:
    MateExtractor(const std::string& htsFilePath, const std::string& htsReferencePath, const bool cacheMates,
                  const int farAwayMateDistanceThreshold = 1000);
    ~MateExtractor();

    boost::optional<FullRead> extractMate(const ReadId& mateReadId, const GenomicRegion& genomicRegion);
    std::vector<FullRead> extractMates(const MateRegionToRecover& mateRegionToRecover);
    void addMateToCache(const ReadId& readId, const FullRead& mate);
    size_t mateCacheSize() const { return mateCache_.size(); }
    uint64_t extractedMatesCounter() const { return extractedTotalCounter_; }
    uint64_t extractedMatesFromDiskCounter() const { return extractedFromDiskCounter_; }
    const MateCache& mateCache() const { return mateCache_; }
    MateCache& mutableMateCache() { return mateCache_; }

private:
    void openFile();
    void loadHeader();
    void loadIndex();

    std::string htsFilePath_;
    std::string htsReferencePath_;
    ReferenceContigInfo contigInfo_;
    MateCache mateCache_;
    bool cacheMates_;

    htsFile* htsFilePtr_ = nullptr;
    bam_hdr_t* htsHeaderPtr_ = nullptr;
    hts_idx_t* htsIndexPtr_ = nullptr;
    bam1_t* htsAlignmentPtr_ = nullptr;

    int farAwayMateDistanceThreshold_ = 1000;

    uint64_t extractedTotalCounter_ = 0;
    uint64_t extractedFromDiskCounter_ = 0;

};

}

}
