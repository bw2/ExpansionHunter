//
// ExpansionHunter
// Copyright 2016-2020 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Konrad Scheffler <kscheffler@illumina.com>,
//         Felix Schlesinger <fschlesinger@illumina.com>
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

#include <memory>
#include <unordered_set>
#include <vector>

#include "genotyping/FragLogliks.hh"
#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{
namespace strgt
{

class TwoAlleleGenotyper
{
public:
    TwoAlleleGenotyper(
        int motifLen, int fragLen, int readLen, std::vector<double> topFragLogliks, FragLogliks* fragLogliksPtr,
        bool isOptimizedStreamingMode = false)
        : motifLen_(motifLen)
        , fragLen_(fragLen)
        , readLen_(readLen)
        , topFragLogliks_(std::move(topFragLogliks))
        , fragLogliks_(*fragLogliksPtr)
        , isOptimizedStreamingMode_(isOptimizedStreamingMode)
    {
    }

    RepeatGenotype genotype(const std::unordered_set<int>& alleleSizeCandidates);

private:
    using ReadIndexAndNumMotifs = std::pair<int, int>;

    RepeatGenotype getMostLikelyGenotype(const std::unordered_set<int>& alleleSizeCandidates);
    double getShortAndLongAlleleLoglik(int shortAlleleSize, int longAlleleSize);
    double getLongAndShortAlleleLoglik(int longAlleleSize, int shortAlleleSize);

    int motifLen_;
    int fragLen_;
    int readLen_;
    std::vector<double> topFragLogliks_;
    FragLogliks& fragLogliks_;
    bool isOptimizedStreamingMode_;
};

}
}
