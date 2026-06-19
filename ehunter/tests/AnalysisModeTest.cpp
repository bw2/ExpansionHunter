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

#include <stdexcept>
#include <string>

#include "gtest/gtest.h"

#include "core/Parameters.hh"

namespace ehunter
{
// decodeAnalysisMode is an internal (not header-declared) helper in ParameterLoading.cpp with external
// linkage; forward-declare it here to exercise the analysis-mode whitelist directly.
AnalysisMode decodeAnalysisMode(const std::string& encoding);
}

using namespace ehunter;

TEST(DecodeAnalysisMode, ValidModes_Decoded)
{
    EXPECT_EQ(decodeAnalysisMode("seeking"), AnalysisMode::kSeeking);
    EXPECT_EQ(decodeAnalysisMode("streaming"), AnalysisMode::kStreaming);
    EXPECT_EQ(decodeAnalysisMode("low-mem-streaming"), AnalysisMode::kLowMemStreaming);
    EXPECT_EQ(decodeAnalysisMode("optimized-streaming"), AnalysisMode::kOptimizedStreaming);
}

TEST(DecodeAnalysisMode, RegionParallelStreaming_Rejected)
{
    // region-parallel-streaming was removed; it must no longer decode to any mode.
    EXPECT_THROW(decodeAnalysisMode("region-parallel-streaming"), std::logic_error);
}

TEST(DecodeAnalysisMode, UnknownMode_Rejected)
{
    EXPECT_THROW(decodeAnalysisMode("not-a-real-mode"), std::logic_error);
}
