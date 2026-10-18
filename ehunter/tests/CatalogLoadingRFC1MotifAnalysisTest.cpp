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

#include <string>

#include "gtest/gtest.h"

using namespace ehunter;
using Json = nlohmann::json;
using std::string;

namespace
{

// The RFC1 record of Illumina's hg38 catalog, with its ReferenceRegion replaced.
Json makeRFC1Record(const string& referenceRegion)
{
    return Json { { "LocusId", "RFC1" }, { "LocusStructure", "(AARRG)*" }, { "ReferenceRegion", referenceRegion },
        { "VariantType", "Repeat" } };
}

}

TEST(DecodeUseRFC1MotifAnalysis, RFC1LocusIsRecognizedInEachReferenceWithAndWithoutChrPrefix)
{
    for (const string& referenceRegion :
         { "chr4:39348424-39348479", "4:39348424-39348479", "chr4:39350044-39350099", "4:39350044-39350099",
           "chr4:39318077-39318132", "4:39318077-39318132" })
    {
        EXPECT_TRUE(decodeUseRFC1MotifAnalysis(makeRFC1Record(referenceRegion))) << referenceRegion;
    }
}

TEST(DecodeUseRFC1MotifAnalysis, OneElementArraysAreRecognizedLikeStrings)
{
    // loadLocusDescription turns ReferenceRegion and VariantType into arrays before deciding.
    Json record = makeRFC1Record("chr4:39348424-39348479");
    record["ReferenceRegion"] = Json::array({ "chr4:39348424-39348479" });
    record["VariantType"] = Json::array({ "Repeat" });
    EXPECT_TRUE(decodeUseRFC1MotifAnalysis(record));
}

TEST(DecodeUseRFC1MotifAnalysis, FalseTurnsTheAnalysisOff)
{
    Json record = makeRFC1Record("chr4:39348424-39348479");
    record["RFC1MotifAnalysis"] = false;
    EXPECT_FALSE(decodeUseRFC1MotifAnalysis(record));
}

TEST(DecodeUseRFC1MotifAnalysis, OtherValuesOfTheFieldAreIgnored)
{
    for (const Json& value : { Json(true), Json::object(), Json("yes"), Json(0) })
    {
        Json rfc1Record = makeRFC1Record("chr4:39348424-39348479");
        rfc1Record["RFC1MotifAnalysis"] = value;
        EXPECT_TRUE(decodeUseRFC1MotifAnalysis(rfc1Record)) << value;

        Json otherRecord = makeRFC1Record("chr4:1000-1055");
        otherRecord["RFC1MotifAnalysis"] = value;
        EXPECT_FALSE(decodeUseRFC1MotifAnalysis(otherRecord)) << value;
    }
}

TEST(DecodeUseRFC1MotifAnalysis, OtherDefinitionsAtTheRFC1CoordinatesAreNotRecognized)
{
    Json referenceMotif = makeRFC1Record("chr4:39348424-39348479");
    referenceMotif["LocusStructure"] = "(AAAAG)*";
    EXPECT_FALSE(decodeUseRFC1MotifAnalysis(referenceMotif));

    Json rareRepeat = makeRFC1Record("chr4:39348424-39348479");
    rareRepeat["VariantType"] = "RareRepeat";
    EXPECT_FALSE(decodeUseRFC1MotifAnalysis(rareRepeat));

    Json shiftedByOneBase = makeRFC1Record("chr4:39348423-39348479");
    EXPECT_FALSE(decodeUseRFC1MotifAnalysis(shiftedByOneBase));

    Json twoRepeats = makeRFC1Record("chr4:39348424-39348479");
    twoRepeats["LocusStructure"] = "(AARRG)*(AAGGG)*";
    twoRepeats["ReferenceRegion"] = Json::array({ "chr4:39348424-39348479", "chr4:39348479-39348483" });
    twoRepeats["VariantType"] = Json::array({ "Repeat", "Repeat" });
    EXPECT_FALSE(decodeUseRFC1MotifAnalysis(twoRepeats));
}
