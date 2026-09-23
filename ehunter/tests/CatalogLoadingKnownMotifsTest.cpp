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

#include "io/CatalogLoading.hh"

#include <stdexcept>
#include <string>
#include <vector>

#include "gtest/gtest.h"

using namespace ehunter;
using Json = nlohmann::json;
using std::string;
using std::vector;

TEST(DecodeKnownMotifs, ArrayOfStringsIsReturnedAsWritten)
{
    const Json locus = Json::parse(R"({"LocusId": "L", "KnownMotifs": ["CAG", "caa", "CCG"]})");
    EXPECT_EQ(vector<string>({ "CAG", "caa", "CCG" }), decodeKnownMotifs(locus, "L"));
}

TEST(DecodeKnownMotifs, MissingFieldGivesAnEmptyList)
{
    EXPECT_TRUE(decodeKnownMotifs(Json::parse(R"({"LocusId": "L"})"), "L").empty());
}

TEST(DecodeKnownMotifs, ValueThatIsNotAnArrayOfStringsThrows)
{
    EXPECT_THROW(decodeKnownMotifs(Json::parse(R"({"KnownMotifs": "CAG"})"), "L"), std::logic_error);
    EXPECT_THROW(decodeKnownMotifs(Json::parse(R"({"KnownMotifs": ["CAG", 1]})"), "L"), std::logic_error);
}
