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

#include "core/RepeatPurity.hh"

#include <cctype>
#include <cstddef>

namespace ehunter
{

MotifPurity motifTilingPurity(const std::string& sequence, const std::string& motif)
{
    MotifPurity result;
    const std::size_t motifLength = motif.length();
    if (motifLength == 0)
    {
        return result;
    }

    result.totalBases = static_cast<long>(sequence.length());
    for (std::size_t index = 0; index < sequence.length(); ++index)
    {
        const char seqBase = static_cast<char>(std::toupper(static_cast<unsigned char>(sequence[index])));
        const char motifBase = static_cast<char>(std::toupper(static_cast<unsigned char>(motif[index % motifLength])));
        if (seqBase == motifBase)
        {
            ++result.matchedBases;
        }
    }

    return result;
}

}
