//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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
// Adapted from REViewer (Copyright 2020 Illumina, Inc., GPL-3.0 license)
//

#include "reviewer/Origin.hh"

#include <cstdlib>

namespace ehunter
{
namespace reviewer
{

using graphtools::Path;
using std::string;
using std::vector;

FragAssignment getBestFragAssignment(const vector<Path>& /* hapPaths */, const FragPathAlignsById& fragPathAlignsById)
{
    vector<string> fragIds;
    fragIds.reserve(fragPathAlignsById.size());
    for (const auto& idAndFragAligns : fragPathAlignsById)
    {
        fragIds.push_back(idAndFragAligns.first);
    }

    vector<int> alignIndexByFrag(fragIds.size(), 0);
    for (size_t fragIndex = 0; fragIndex != fragIds.size(); ++fragIndex)
    {
        const auto& fragId = fragIds[fragIndex];
        int numOrigins = static_cast<int>(fragPathAlignsById.at(fragId).size());
        int originIndex = rand() % numOrigins;
        alignIndexByFrag[fragIndex] = originIndex;
    }

    return { fragIds, alignIndexByFrag };
}

}  // namespace reviewer
}  // namespace ehunter
