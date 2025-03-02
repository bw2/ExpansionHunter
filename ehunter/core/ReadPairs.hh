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
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "core/Read.hh"

namespace ehunter
{

/**
 * Read pair container class
 */
class FullReadPairs
{
public:
    typedef std::unordered_map<std::string, FullReadPair>::const_iterator const_iterator;
    typedef std::unordered_map<std::string, FullReadPair>::iterator iterator;
    const_iterator begin() const { return fragmentIdToReadPair_.begin(); }
    const_iterator end() const { return fragmentIdToReadPair_.end(); }
    iterator begin() { return fragmentIdToReadPair_.begin(); }
    iterator end() { return fragmentIdToReadPair_.end(); }

    FullReadPairs() = default;
    void Clear();
    void Add(FullRead read);
    void AddMateToExistingRead(FullRead mate);

    const FullReadPair& operator[](const std::string& fragmentId) const;

    int32_t NumReads() const { return numReads_; }
    int32_t NumCompletePairs() const;

    bool operator==(const FullReadPairs& other) const
    {
        return (fragmentIdToReadPair_ == other.fragmentIdToReadPair_ && numReads_ == other.numReads_);
    }

private:
    std::unordered_map<std::string, FullReadPair> fragmentIdToReadPair_;
    int32_t numReads_ = 0;
};

}
