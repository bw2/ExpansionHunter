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

#include "core/ReadPairs.hh"

#include <stdexcept>

using std::string;
using std::vector;

namespace ehunter
{

void FullReadPairs::Add(FullRead read)
{
    FullReadPair& readPair = fragmentIdToReadPair_[read.r.fragmentId()];
    const int originalMateCount = readPair.numMatesSet();

    if (read.r.isFirstMate() && readPair.firstMate == std::nullopt)
    {
        readPair.firstMate = std::move(read);
    }

    if (read.r.isSecondMate() && readPair.secondMate == std::nullopt)
    {
        readPair.secondMate = std::move(read);
    }

    const int mateCountAfterAdd = readPair.numMatesSet();

    numReads_ += mateCountAfterAdd - originalMateCount;
}

void FullReadPairs::AddMateToExistingRead(FullRead mate)
{
    FullReadPair& readPair = fragmentIdToReadPair_.at(mate.r.fragmentId());
    if (mate.r.isFirstMate() && readPair.firstMate == std::nullopt)
    {
        readPair.firstMate = std::move(mate);
        ++numReads_;
    }
    else if (mate.r.isSecondMate() && readPair.secondMate == std::nullopt)
    {
        readPair.secondMate = std::move(mate);
        ++numReads_;
    }
    else
    {
        throw std::logic_error("Unable to find read placement");
    }
}

const FullReadPair& FullReadPairs::operator[](const string& fragment_id) const
{
    if (fragmentIdToReadPair_.find(fragment_id) == fragmentIdToReadPair_.end())
    {
        throw std::logic_error("Fragment " + fragment_id + " does not exist");
    }
    return fragmentIdToReadPair_.at(fragment_id);
}

int32_t FullReadPairs::NumCompletePairs() const
{
    int32_t numCompletePairs = 0;
    for (const auto& fragmentIdAndReads : fragmentIdToReadPair_)
    {
        const FullReadPair& readPair = fragmentIdAndReads.second;
        if (readPair.firstMate && readPair.secondMate)
        {
            ++numCompletePairs;
        }
    }

    return numCompletePairs;
}

void FullReadPairs::Clear()
{
    fragmentIdToReadPair_.clear();
    numReads_ = 0;
}


}
