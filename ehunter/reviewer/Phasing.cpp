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

#include "reviewer/Phasing.hh"

#include <algorithm>
#include <cassert>

#include "reviewer/Projection.hh"

namespace ehunter
{
namespace reviewer
{

using std::pair;
using std::string;
using std::vector;

static int scorePath(int pathIndex, const PairPathAlignById& pairPathAlignById)
{
    int pathScore = 0;
    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& pairAlign = idAndPairPathAlign.second;
        for (const auto& readAlign : pairAlign.readAligns)
        {
            if (readAlign.pathIndex == pathIndex)
            {
                pathScore += score(*readAlign.align);
                break;
            }
        }

        for (const auto& mateAlign : pairAlign.mateAligns)
        {
            if (mateAlign.pathIndex == pathIndex)
            {
                pathScore += score(*mateAlign.align);
                break;
            }
        }
    }

    return pathScore;
}

ScoredDiplotypes scoreDiplotypes(const FragById& fragById, const vector<Diplotype>& diplotypes)
{
    vector<ScoredDiplotype> scoredDiplotypes;

    for (const auto& diplotype : diplotypes)
    {
        assert(diplotype.size() == 1 || diplotype.size() == 2);

        auto pairPathAlignById = project(diplotype, fragById);

        int genotypeScore = scorePath(0, pairPathAlignById);
        if (diplotype.size() == 2)
        {
            genotypeScore += scorePath(1, pairPathAlignById);
        }

        scoredDiplotypes.emplace_back(diplotype, genotypeScore);
    }

    std::sort(
        scoredDiplotypes.begin(), scoredDiplotypes.end(),
        [](const ScoredDiplotype& gt1, const ScoredDiplotype& gt2) { return gt1.second > gt2.second; });

    assert(!scoredDiplotypes.empty());
    return scoredDiplotypes;
}

}  // namespace reviewer
}  // namespace ehunter
