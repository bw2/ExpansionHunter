//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
// Concept: Michael Eberle <meberle@illumina.com>
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

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/optional.hpp>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/Reference.hh"
#include "io/LocusSpecDecoding.hh"

using boost::optional;
using graphtools::NodeId;
using std::map;
using std::ostream;
using std::string;
using std::to_string;
using std::vector;

using Json = nlohmann::json;

namespace spd = spdlog;

namespace ehunter
{

static bool checkIfFieldExists(const Json& record, const string& fieldName)
{
    return record.find(fieldName) != record.end();
}

static void assertFieldExists(const Json& record, const string& fieldName)
{
    if (!checkIfFieldExists(record, fieldName))
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Field " + fieldName + " must be present in " + out.str());
    }
}

static void assertRecordIsArray(const Json& record)
{
    if (!record.is_array())
    {
        std::stringstream out;
        out << record;
        throw std::logic_error("Expected array but got this instead " + out.str());
    }
}

static void makeArray(Json& record)
{
    if (record.type() != Json::value_t::array)
    {
        record = Json::array({ record });
    }
}

static VariantTypeFromUser decodeVariantTypeFromUser(const string& encoding)
{
    if (encoding == "RareRepeat")
    {
        return VariantTypeFromUser::kRareRepeat;
    }
    if (encoding == "Repeat")
    {
        return VariantTypeFromUser::kCommonRepeat;
    }
    if (encoding == "SmallVariant")
    {
        return VariantTypeFromUser::kSmallVariant;
    }
    if (encoding == "SMN")
    {
        return VariantTypeFromUser::kSMN;
    }
    else
    {
        throw std::logic_error("Encountered invalid variant type: " + encoding);
    }
}

static vector<string> generateIds(const std::string& locusId, const Json& variantRegionEncodings)
{
    if (variantRegionEncodings.size() == 1)
    {
        return { locusId };
    }

    vector<string> variantIds;
    for (const auto& regionEncoding : variantRegionEncodings)
    {
        variantIds.push_back(locusId + "_" + regionEncoding.get<string>());
    }

    return variantIds;
}

/// \brief Translate a single locus from the catalog file json structure into an intermediate locus configuration
///
static LocusDescriptionFromUser loadUserDescription(Json& locusJson, const ReferenceContigInfo& contigInfo)
{
    LocusDescriptionFromUser userDescription;

    assertFieldExists(locusJson, "LocusId");
    userDescription.locusId = locusJson["LocusId"].get<string>();

    assertFieldExists(locusJson, "ReferenceRegion");
    makeArray(locusJson["ReferenceRegion"]);
    for (const auto& encoding : locusJson["ReferenceRegion"])
    {
        GenomicRegion region = decode(contigInfo, encoding.get<string>());
        userDescription.referenceRegions.push_back(region);
    }

    assertFieldExists(locusJson, "LocusStructure");
    userDescription.locusStructure = locusJson["LocusStructure"].get<string>();

    assertFieldExists(locusJson, "VariantType");
    makeArray(locusJson["VariantType"]);
    for (const auto& encoding : locusJson["VariantType"])
    {
        userDescription.variantTypesFromUser.push_back(decodeVariantTypeFromUser(encoding.get<string>()));
    }

    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        for (const auto& locusEncoding : locusJson["TargetRegion"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.targetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "VariantId"))
    {
        makeArray(locusJson["VariantId"]);
        for (const auto& variantId : locusJson["VariantId"])
        {
            userDescription.variantIds.push_back(variantId.get<string>());
        }
    }
    else
    {
        userDescription.variantIds = generateIds(userDescription.locusId, locusJson["ReferenceRegion"]);
    }

    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const auto& locusEncoding : locusJson["OfftargetRegions"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            userDescription.offtargetRegions.push_back(region);
        }
    }

    if (checkIfFieldExists(locusJson, "ErrorRate"))
    {
        userDescription.errorRate = locusJson["ErrorRate"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "LikelihoodRatioThreshold"))
    {
        userDescription.likelihoodRatioThreshold = locusJson["LikelihoodRatioThreshold"].get<double>();
    }
    if (checkIfFieldExists(locusJson, "MinimalLocusCoverage"))
    {
        userDescription.minLocusCoverage = locusJson["MinimalLocusCoverage"].get<double>();
    }

    static const std::string rfc1MotifAnalysisKey("RFC1MotifAnalysis");
    if (checkIfFieldExists(locusJson, rfc1MotifAnalysisKey))
    {
        const Json& record(locusJson[rfc1MotifAnalysisKey]);
        if (record.type() == Json::value_t::boolean)
        {
            userDescription.useRFC1MotifAnalysis = record.get<bool>();
        }
        else if (record.type() == Json::value_t::object)
        {
            userDescription.useRFC1MotifAnalysis = true;
        }
        else
        {
            std::stringstream out;
            out << record;
            throw std::logic_error(
                "Key '" + rfc1MotifAnalysisKey
                + "' must have either a boolean or object value type, observed value is '" + out.str() + "'");
        }
    }

    return userDescription;
}

void printMessageIfLocusIdArgNotFullyProcessed(vector<std::string>& locusIdsFromArg,
    vector<std::string>& processedLocusIds, const std::string& locusIdArg) {

    if (!locusIdsFromArg.empty() && processedLocusIds.size() < locusIdsFromArg.size()) {
        if (processedLocusIds.size() == 0) {
            if (locusIdsFromArg.size() > 1) {
                spdlog::error("ERROR: --locus arg: none of the locus ids {} were not found in the catalog", locusIdArg);
            } else {
                spdlog::error("ERROR: --locus arg: {} not found in the catalog", locusIdArg);
            }
        } else {
            std::string message = "";
            int i = 0;
            for (const auto& locusId : locusIdsFromArg) {
                if (std::find(processedLocusIds.begin(), processedLocusIds.end(), locusId) == processedLocusIds.end()) {
                    if (i > 0) {
                        message += ", ";
                    }
                    message += locusId;
                    i += 1;
                }
            }
            if (i > 1) {
                spdlog::warn("WARNING: --locus arg: {} locus ids were not found in the catalog: {}", i, message);
            } else {
                spdlog::warn("WARNING: --locus arg: 1 locus id was not found in the catalog: {}", message);
            }

        }
    }
}

RegionCatalog loadLocusCatalogFromDisk(
    const string& catalogPath, const std::string& locusIdArg, const HeuristicParameters& heuristicParams,
    const Reference& reference)
{
    vector<std::string> locusIdsFromArg;
    if (!locusIdArg.empty()) {
        boost::split(locusIdsFromArg, locusIdArg, boost::is_any_of(","));
    }

    std::ifstream inputStream(catalogPath.c_str());

    if (!inputStream.is_open())
    {
        throw std::runtime_error("Failed to open catalog file " + catalogPath);
    }

    Json catalogJson;
    if (boost::algorithm::ends_with(catalogPath, "gz"))
    {
        std::ifstream binaryInputStream(catalogPath.c_str(), std::ios::binary);
        boost::iostreams::filtering_istreambuf bufferedInputStream;
        bufferedInputStream.push(boost::iostreams::gzip_decompressor());
        bufferedInputStream.push(binaryInputStream);
        std::istream(&bufferedInputStream) >> catalogJson;

    } else {
        std::ifstream inputStream(catalogPath.c_str());

        if (!inputStream.is_open())
        {
            throw std::runtime_error("Failed to open catalog file " + catalogPath);
        }
        inputStream >> catalogJson;
    }
    makeArray(catalogJson);

    RegionCatalog catalog;
    catalog.reserve(catalogJson.size());
    vector<std::string> processedLocusIds;
    for (auto& locusJson : catalogJson)
    {
        if (!locusIdsFromArg.empty()) {
            // check if currentLocusId is in the vector of --locus ids
            assertFieldExists(locusJson, "LocusId");
            std::string currentLocusId = locusJson["LocusId"].get<string>();
            if (std::find(locusIdsFromArg.begin(), locusIdsFromArg.end(), currentLocusId) == locusIdsFromArg.end()) {
                continue;
            }
            processedLocusIds.push_back(currentLocusId);
        }

        LocusDescriptionFromUser userDescription = loadUserDescription(locusJson, reference.contigInfo());
        try {
            LocusSpecification locusSpec = decodeLocusSpecification(userDescription, reference, heuristicParams);
            catalog.emplace_back(std::move(locusSpec));
        }
        catch (const std::exception& e)
        {
            std::cout << "Error on locus spec " + userDescription.locusId + ": " + e.what() << std::endl;
        }
    }

    printMessageIfLocusIdArgNotFullyProcessed(locusIdsFromArg, processedLocusIds, locusIdArg);

    return catalog;
}

}
