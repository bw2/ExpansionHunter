//
// ExpansionHunter
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

#include "io/CatalogLoading.hh"

#include <stdexcept>

#include "gtest/gtest.h"

#include "core/LocusStats.hh"

using namespace ehunter;
using Json = nlohmann::json;

namespace
{

Json makeRecord()
{
    return Json { { "LocusId", "PAR_LOCUS" }, { "LocusStructure", "(CAG)*" }, { "ReferenceRegion", "chrX:3000-3036" },
        { "VariantType", "Repeat" } };
}

Json makeRecordWithField(const Json& value)
{
    Json record = makeRecord();
    record["IsInsidePseudoautosomalRegion"] = value;
    return record;
}

}

TEST(DecodeIsInsidePseudoautosomalRegion, AbsentFieldMeansFalse)
{
    EXPECT_FALSE(decodeIsInsidePseudoautosomalRegion(makeRecord(), "PAR_LOCUS", ChromType::kX));
    EXPECT_FALSE(decodeIsInsidePseudoautosomalRegion(makeRecord(), "PAR_LOCUS", ChromType::kAutosome));
}

TEST(DecodeIsInsidePseudoautosomalRegion, TrueIsAcceptedOnChrX)
{
    EXPECT_TRUE(decodeIsInsidePseudoautosomalRegion(makeRecordWithField(true), "PAR_LOCUS", ChromType::kX));
}

TEST(DecodeIsInsidePseudoautosomalRegion, FalseIsAcceptedOnAnyChromosome)
{
    for (ChromType contigChromType : { ChromType::kX, ChromType::kY, ChromType::kAutosome })
    {
        EXPECT_FALSE(decodeIsInsidePseudoautosomalRegion(makeRecordWithField(false), "PAR_LOCUS", contigChromType));
    }
}

TEST(DecodeIsInsidePseudoautosomalRegion, TrueOnChrYOrAnAutosomeIsAnError)
{
    for (ChromType contigChromType : { ChromType::kY, ChromType::kAutosome })
    {
        EXPECT_THROW(
            decodeIsInsidePseudoautosomalRegion(makeRecordWithField(true), "PAR_LOCUS", contigChromType),
            std::runtime_error);
    }
}

TEST(DecodeIsInsidePseudoautosomalRegion, ValuesOtherThanTrueOrFalseAreAnError)
{
    for (const Json& value : { Json("true"), Json(1), Json(nullptr), Json::array({ true }) })
    {
        EXPECT_THROW(
            decodeIsInsidePseudoautosomalRegion(makeRecordWithField(value), "PAR_LOCUS", ChromType::kX),
            std::runtime_error)
            << value;
    }
}

TEST(DetermineExpectedAlleleCount, PseudoautosomalLociHaveTwoAllelesInEverySample)
{
    EXPECT_EQ(AlleleCount::kTwo, determineExpectedAlleleCount(ChromType::kPseudoautosomalRegion, Sex::kMale));
    EXPECT_EQ(AlleleCount::kTwo, determineExpectedAlleleCount(ChromType::kPseudoautosomalRegion, Sex::kFemale));
}
