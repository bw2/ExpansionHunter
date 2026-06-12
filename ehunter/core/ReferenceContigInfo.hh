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

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "core/Common.hh"

namespace ehunter
{

// Thrown when a contig name is not present in the reference. Derives from std::logic_error so existing
// `catch (const std::logic_error&)` handlers keep working, while allowing callers to catch this case
// specifically (e.g. to warn and skip a catalog locus rather than aborting the whole run).
class MissingContigError : public std::logic_error
{
public:
    explicit MissingContigError(const std::string& contigName)
        : std::logic_error("Invalid contig name " + contigName)
    {
    }
};

// Handles translation between contig names and indexes
class ReferenceContigInfo
{
public:
    explicit ReferenceContigInfo(std::vector<std::pair<std::string, int64_t>> namesAndSizes);

    int32_t numContigs() const { return namesAndSizes_.size(); }
    const std::string& getContigName(int32_t contigIndex) const;
    int64_t getContigSize(int32_t contigIndex) const;
    const ChromType& getChromType(int32_t contigIndex) const;
    int32_t getContigId(const std::string& contigName) const;

private:
    void assertValidIndex(int32_t contigIndex) const;

    std::vector<ChromType> chromTypes_;
    std::vector<std::pair<std::string, int64_t>> namesAndSizes_;
    std::unordered_map<std::string, int32_t> nameToIndex_;
};

std::ostream& operator<<(std::ostream& out, const ReferenceContigInfo& contigInfo);

}
