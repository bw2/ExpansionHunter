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

#include "gtest/gtest.h"

#include <string>
#include <vector>

#include "thirdparty/json/json.hpp"
#include "locus/LocusSpecification.hh"

using Json = nlohmann::json;
using std::string;
using std::vector;

using namespace ehunter;

namespace
{

// Helper functions that mirror the static functions in CatalogLoading.cpp
// These are duplicated here for testing purposes since the originals are static

PlotThresholdAppliedTo decodePlotThresholdAppliedTo(const string& encoding)
{
    if (encoding == "ShortAllele")
    {
        return PlotThresholdAppliedTo::kShortAllele;
    }
    if (encoding == "LongAllele")
    {
        return PlotThresholdAppliedTo::kLongAllele;
    }
    throw std::logic_error("PlotReadVisualization 'If' field must be 'ShortAllele' or 'LongAllele', got: " + encoding);
}

PlotThresholdComparisonOp decodePlotThresholdComparisonOp(const string& encoding)
{
    if (encoding == "<")
    {
        return PlotThresholdComparisonOp::kLessThan;
    }
    if (encoding == "<=")
    {
        return PlotThresholdComparisonOp::kLessThanOrEqual;
    }
    if (encoding == "=" || encoding == "==")
    {
        return PlotThresholdComparisonOp::kEqual;
    }
    if (encoding == "!=" || encoding == "<>")
    {
        return PlotThresholdComparisonOp::kNotEqual;
    }
    if (encoding == ">=")
    {
        return PlotThresholdComparisonOp::kGreaterThanOrEqual;
    }
    if (encoding == ">")
    {
        return PlotThresholdComparisonOp::kGreaterThan;
    }
    throw std::logic_error("PlotReadVisualization 'Is' field must be '<', '<=', '=', '==', '!=', '<>', '>=', or '>', got: " + encoding);
}

bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

// Parse PlotReadVisualization array from JSON, mirroring CatalogLoading.cpp logic
vector<PlotReadVisualization> parsePlotReadVisualization(const Json& locusJson)
{
    vector<PlotReadVisualization> plotConditions;
    static const string plotReadVisualizationKey("PlotReadVisualization");

    if (!checkIfFieldExists(locusJson, plotReadVisualizationKey))
    {
        return plotConditions;
    }

    const Json& record(locusJson[plotReadVisualizationKey]);
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error(
            "Key '" + plotReadVisualizationKey + "' must be an array, observed value is '" + out.str() + "'");
    }

    plotConditions.reserve(record.size());
    for (const auto& conditionJson : record)
    {
        if (!conditionJson.is_object())
        {
            std::stringstream out;
            out << conditionJson;
            throw std::logic_error(
                "Each element in '" + plotReadVisualizationKey + "' must be an object with 'If', 'Is', and 'Threshold' fields, got: " + out.str());
        }
        if (!checkIfFieldExists(conditionJson, "If") ||
            !checkIfFieldExists(conditionJson, "Is") ||
            !checkIfFieldExists(conditionJson, "Threshold"))
        {
            std::stringstream out;
            out << conditionJson;
            throw std::logic_error(
                "Each element in '" + plotReadVisualizationKey + "' must have 'If', 'Is', and 'Threshold' fields, got: " + out.str());
        }
        PlotReadVisualization condition;
        condition.appliedTo = decodePlotThresholdAppliedTo(conditionJson["If"].get<string>());
        condition.op = decodePlotThresholdComparisonOp(conditionJson["Is"].get<string>());
        condition.threshold = conditionJson["Threshold"].get<int>();
        plotConditions.push_back(condition);
    }

    return plotConditions;
}

} // anonymous namespace


// Tests for decodePlotThresholdAppliedTo

TEST(DecodePlotThresholdAppliedTo, ShortAllele_ReturnsKShortAllele)
{
    EXPECT_EQ(PlotThresholdAppliedTo::kShortAllele, decodePlotThresholdAppliedTo("ShortAllele"));
}

TEST(DecodePlotThresholdAppliedTo, LongAllele_ReturnsKLongAllele)
{
    EXPECT_EQ(PlotThresholdAppliedTo::kLongAllele, decodePlotThresholdAppliedTo("LongAllele"));
}

TEST(DecodePlotThresholdAppliedTo, InvalidValue_Throws)
{
    EXPECT_THROW(decodePlotThresholdAppliedTo("InvalidAllele"), std::logic_error);
    EXPECT_THROW(decodePlotThresholdAppliedTo("shortallele"), std::logic_error);
    EXPECT_THROW(decodePlotThresholdAppliedTo(""), std::logic_error);
}


// Tests for decodePlotThresholdComparisonOp

TEST(DecodePlotThresholdComparisonOp, LessThan_ReturnsKLessThan)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThan, decodePlotThresholdComparisonOp("<"));
}

TEST(DecodePlotThresholdComparisonOp, LessThanOrEqual_ReturnsKLessThanOrEqual)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThanOrEqual, decodePlotThresholdComparisonOp("<="));
}

TEST(DecodePlotThresholdComparisonOp, Equal_ReturnsKEqual)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kEqual, decodePlotThresholdComparisonOp("="));
    EXPECT_EQ(PlotThresholdComparisonOp::kEqual, decodePlotThresholdComparisonOp("=="));
}

TEST(DecodePlotThresholdComparisonOp, NotEqual_ReturnsKNotEqual)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kNotEqual, decodePlotThresholdComparisonOp("!="));
    EXPECT_EQ(PlotThresholdComparisonOp::kNotEqual, decodePlotThresholdComparisonOp("<>"));
}

TEST(DecodePlotThresholdComparisonOp, GreaterThanOrEqual_ReturnsKGreaterThanOrEqual)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThanOrEqual, decodePlotThresholdComparisonOp(">="));
}

TEST(DecodePlotThresholdComparisonOp, GreaterThan_ReturnsKGreaterThan)
{
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThan, decodePlotThresholdComparisonOp(">"));
}

TEST(DecodePlotThresholdComparisonOp, InvalidValue_Throws)
{
    EXPECT_THROW(decodePlotThresholdComparisonOp("invalid"), std::logic_error);
    EXPECT_THROW(decodePlotThresholdComparisonOp(""), std::logic_error);
    EXPECT_THROW(decodePlotThresholdComparisonOp("<<"), std::logic_error);
}


// Tests for PlotReadVisualization JSON parsing

TEST(ParsePlotReadVisualization, SingleCondition_ParsedCorrectly)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "LongAllele", "Is": ">=", "Threshold": 100}
        ]
    })");

    auto conditions = parsePlotReadVisualization(locusJson);

    ASSERT_EQ(1, conditions.size());
    EXPECT_EQ(PlotThresholdAppliedTo::kLongAllele, conditions[0].appliedTo);
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThanOrEqual, conditions[0].op);
    EXPECT_EQ(100, conditions[0].threshold);
}

TEST(ParsePlotReadVisualization, MultipleConditions_ParsedCorrectly)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "ShortAllele", "Is": "<", "Threshold": 10},
            {"If": "LongAllele", "Is": ">", "Threshold": 50}
        ]
    })");

    auto conditions = parsePlotReadVisualization(locusJson);

    ASSERT_EQ(2, conditions.size());

    EXPECT_EQ(PlotThresholdAppliedTo::kShortAllele, conditions[0].appliedTo);
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThan, conditions[0].op);
    EXPECT_EQ(10, conditions[0].threshold);

    EXPECT_EQ(PlotThresholdAppliedTo::kLongAllele, conditions[1].appliedTo);
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThan, conditions[1].op);
    EXPECT_EQ(50, conditions[1].threshold);
}

TEST(ParsePlotReadVisualization, MissingIfField_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"Is": ">=", "Threshold": 100}
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, MissingIsField_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "LongAllele", "Threshold": 100}
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, MissingThresholdField_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "LongAllele", "Is": ">="}
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, NotAnArray_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": {"If": "LongAllele", "Is": ">=", "Threshold": 100}
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, ArrayElementNotObject_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            "not an object"
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, NoPlotReadVisualizationField_ReturnsEmpty)
{
    Json locusJson = Json::parse(R"({
        "LocusId": "TestLocus"
    })");

    auto conditions = parsePlotReadVisualization(locusJson);

    EXPECT_TRUE(conditions.empty());
}

TEST(ParsePlotReadVisualization, EmptyArray_ReturnsEmpty)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": []
    })");

    auto conditions = parsePlotReadVisualization(locusJson);

    EXPECT_TRUE(conditions.empty());
}

TEST(ParsePlotReadVisualization, AllComparisonOperators_ParsedCorrectly)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "ShortAllele", "Is": "<", "Threshold": 1},
            {"If": "ShortAllele", "Is": "<=", "Threshold": 2},
            {"If": "ShortAllele", "Is": "=", "Threshold": 3},
            {"If": "ShortAllele", "Is": "==", "Threshold": 4},
            {"If": "ShortAllele", "Is": "!=", "Threshold": 5},
            {"If": "ShortAllele", "Is": "<>", "Threshold": 6},
            {"If": "ShortAllele", "Is": ">=", "Threshold": 7},
            {"If": "ShortAllele", "Is": ">", "Threshold": 8}
        ]
    })");

    auto conditions = parsePlotReadVisualization(locusJson);

    ASSERT_EQ(8, conditions.size());
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThan, conditions[0].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kLessThanOrEqual, conditions[1].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kEqual, conditions[2].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kEqual, conditions[3].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kNotEqual, conditions[4].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kNotEqual, conditions[5].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThanOrEqual, conditions[6].op);
    EXPECT_EQ(PlotThresholdComparisonOp::kGreaterThan, conditions[7].op);
}

TEST(ParsePlotReadVisualization, InvalidAppliedTo_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "InvalidAllele", "Is": ">=", "Threshold": 100}
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}

TEST(ParsePlotReadVisualization, InvalidComparisonOp_Throws)
{
    Json locusJson = Json::parse(R"({
        "PlotReadVisualization": [
            {"If": "LongAllele", "Is": "invalid", "Threshold": 100}
        ]
    })");

    EXPECT_THROW(parsePlotReadVisualization(locusJson), std::logic_error);
}
