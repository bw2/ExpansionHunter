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

#include "reviewer/FragLenFilter.hh"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace ehunter
{
namespace reviewer
{

static int calcFragLen(const ReadPathAlign& readAlign, const ReadPathAlign& mateAlign)
{
    if (readAlign.pathIndex != mateAlign.pathIndex)
    {
        return std::numeric_limits<int>::max();
    }

    return std::max(readAlign.end, mateAlign.end) - std::min(readAlign.begin, mateAlign.begin);
}

FragPathAlignsById resolveByFragLen(int meanFragLen, const Diplotype& paths, const PairPathAlignById& pairPathAlignById)
{
    FragPathAlignsById fragPathAlignsById;

    for (const auto& idAndPairPathAlign : pairPathAlignById)
    {
        const auto& fragId = idAndPairPathAlign.first;
        const PairPathAlign& pairPathAlign = idAndPairPathAlign.second;

        int bestFragLen = std::numeric_limits<int>::max();
        for (const auto& readAlign : pairPathAlign.readAligns)
        {
            for (const auto& mateAlign : pairPathAlign.mateAligns)
            {
                if (readAlign.pathIndex != mateAlign.pathIndex)
                {
                    continue;
                }

                // Suppress unused variable warning - paths is kept for API consistency
                (void)paths;
                const int fragLen = calcFragLen(readAlign, mateAlign);

                if (std::abs(fragLen - meanFragLen) < std::abs(bestFragLen - meanFragLen))
                {
                    bestFragLen = fragLen;
                    fragPathAlignsById[fragId].clear();
                }

                if (fragLen == bestFragLen)
                {
                    fragPathAlignsById[fragId].emplace_back(readAlign, mateAlign);
                }

                // Legacy behavior (standalone REViewer): sets bestFragLen to meanFragLen, so the equality check
                // against meanFragLen always passes and all same-path alignments are retained. This is a bug because
                // it defeats the fragment-length filtering step; replaced by the code above.
                // if (std::abs(fragLen - meanFragLen) < std::abs(bestFragLen - meanFragLen))
                // {
                //     bestFragLen = meanFragLen;
                //     fragPathAlignsById[fragId].clear();
                // }
                //
                // if (meanFragLen == bestFragLen)
                // {
                //     fragPathAlignsById[fragId].emplace_back(readAlign, mateAlign);
                // }
            }
        }
    }

    return fragPathAlignsById;
}

int getMeanFragLen(const FragById& fragById)
{
    if (fragById.empty())
    {
        throw std::runtime_error("There are no read alignments in the target region");
    }

    graphtools::NodeId leftFlankId = 0;
    const auto& firstFrag = fragById.begin()->second;
    graphtools::NodeId rightFlankId = firstFrag.read.align.path().graphRawPtr()->numNodes() - 1;

    double fragLenSum = 0;
    int numFlankingReads = 0;
    for (const auto& idAndFrag : fragById)
    {
        const auto& frag = idAndFrag.second;

        const auto readStartNode = frag.read.align.path().getNodeIdByIndex(0);
        const auto mateStartNode = frag.mate.align.path().getNodeIdByIndex(0);

        const bool matesStartOnLeftFlank = readStartNode == leftFlankId && mateStartNode == leftFlankId;
        const bool matesStartOnRightFlank = readStartNode == rightFlankId && mateStartNode == rightFlankId;

        if (!matesStartOnLeftFlank && !matesStartOnRightFlank)
        {
            continue;
        }

        const int readStart = static_cast<int>(frag.read.align.path().startPosition());
        const int readEnd = readStart + static_cast<int>(frag.read.align.queryLength());

        const int mateStart = static_cast<int>(frag.mate.align.path().startPosition());
        const int mateEnd = mateStart + static_cast<int>(frag.mate.align.queryLength());

        fragLenSum += readEnd <= mateEnd ? mateEnd - readStart : readEnd - mateStart;
        ++numFlankingReads;
    }

    if (numFlankingReads == 0)
    {
        throw std::runtime_error("Unable to determine fragment length due to missing flanking read fragments");
    }

    return static_cast<int>(fragLenSum / numFlankingReads);
}

}  // namespace reviewer
}  // namespace ehunter
