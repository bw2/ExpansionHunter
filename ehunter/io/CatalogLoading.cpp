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
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>

#include "core/Common.hh"
#include "core/Reference.hh"
#include "graphcore/GraphReferenceMapping.hh"
#include "io/CatalogLoading.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/StringUtils.hh"
#include "locus/LocusSpecification.hh"
#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

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

static PlotThresholdAppliedTo decodePlotThresholdAppliedTo(const string& encoding)
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

static PlotThresholdComparisonOp decodePlotThresholdComparisonOp(const string& encoding)
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
static LocusDescription loadLocusDescription(
    Json& locusJson, const ReferenceContigInfo& contigInfo, const HeuristicParameters& heuristicParams)
{
    assertFieldExists(locusJson, "LocusId");
    std::string locusId = locusJson["LocusId"].get<string>();

    std::vector<GenomicRegion> referenceRegions;
    assertFieldExists(locusJson, "ReferenceRegion");
    makeArray(locusJson["ReferenceRegion"]);
    referenceRegions.reserve(locusJson["ReferenceRegion"].size());
    for (const auto& encoding : locusJson["ReferenceRegion"])
    {
        GenomicRegion region = decode(contigInfo, encoding.get<string>());
        referenceRegions.push_back(region);
    }

    if (referenceRegions.empty()) {
        throw std::runtime_error("No ReferenceRegion specified for locus " + locusId);
    }

    assertFieldExists(locusJson, "LocusStructure");
    std::string locusStructure = locusJson["LocusStructure"].get<string>();

    std::vector<VariantTypeFromUser> variantTypesFromUser;
    assertFieldExists(locusJson, "VariantType");
    makeArray(locusJson["VariantType"]);
    variantTypesFromUser.reserve(locusJson["VariantType"].size());
    for (const auto& encoding : locusJson["VariantType"])
    {
        variantTypesFromUser.push_back(decodeVariantTypeFromUser(encoding.get<string>()));
    }

    std::vector<GenomicRegion> targetRegions;
    if (checkIfFieldExists(locusJson, "TargetRegion"))
    {
        makeArray(locusJson["TargetRegion"]);
        targetRegions.reserve(locusJson["TargetRegion"].size());
        for (const auto& locusEncoding : locusJson["TargetRegion"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            targetRegions.push_back(region);
        }
    }

    std::vector<std::string> variantIds;
    if (checkIfFieldExists(locusJson, "VariantId"))
    {
        makeArray(locusJson["VariantId"]);
        variantIds.reserve(locusJson["VariantId"].size());
        for (const auto& variantId : locusJson["VariantId"])
        {
            variantIds.push_back(variantId.get<string>());
        }
    }
    else
    {
        variantIds = generateIds(locusId, locusJson["ReferenceRegion"]);
    }

    std::vector<GenomicRegion> offtargetRegions;
    if (checkIfFieldExists(locusJson, "OfftargetRegions"))
    {
        offtargetRegions.reserve(locusJson["OfftargetRegions"].size());
        assertRecordIsArray(locusJson["OfftargetRegions"]);
        for (const auto& locusEncoding : locusJson["OfftargetRegions"])
        {
            GenomicRegion region = decode(contigInfo, locusEncoding.get<string>());
            offtargetRegions.push_back(region);
        }
    }

    std::optional<double> errorRate;
    if (checkIfFieldExists(locusJson, "ErrorRate"))
    {
        errorRate = locusJson["ErrorRate"].get<double>();
    }

    std::optional<double> likelihoodRatioThreshold;
    if (checkIfFieldExists(locusJson, "LikelihoodRatioThreshold"))
    {
        likelihoodRatioThreshold = locusJson["LikelihoodRatioThreshold"].get<double>();
    }

    std::optional<double> minimalLocusCoverage;
    if (checkIfFieldExists(locusJson, "MinimalLocusCoverage"))
    {
        minimalLocusCoverage = locusJson["MinimalLocusCoverage"].get<double>();
    }

    bool useRFC1MotifAnalysis = false;
    static const std::string rfc1MotifAnalysisKey("RFC1MotifAnalysis");
    if (checkIfFieldExists(locusJson, rfc1MotifAnalysisKey))
    {
        const Json& record(locusJson[rfc1MotifAnalysisKey]);
        if (record.type() == Json::value_t::boolean)
        {
            useRFC1MotifAnalysis = record.get<bool>();
        }
        else if (record.type() == Json::value_t::object)
        {
            useRFC1MotifAnalysis = true;
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

    std::vector<PlotReadVisualization> plotConditions;
    static const std::string plotReadVisualizationKey("PlotReadVisualization");
    if (checkIfFieldExists(locusJson, plotReadVisualizationKey))
    {
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
    }

    // Define known fields that ExpansionHunter processes
    static const std::set<std::string> knownFields = {
        "LocusId", "ReferenceRegion", "LocusStructure", "VariantType",
        "TargetRegion", "VariantId", "OfftargetRegions", "ErrorRate",
        "LikelihoodRatioThreshold", "MinimalLocusCoverage",
        "RFC1MotifAnalysis", "PlotReadVisualization"
    };

    // Capture any extra fields not in the known set
    Json extraFields;
    for (auto it = locusJson.begin(); it != locusJson.end(); ++it) {
        if (knownFields.find(it.key()) == knownFields.end()) {
            extraFields[it.key()] = it.value();
        }
    }

    const int64_t kExtensionLength = heuristicParams.regionExtensionLength();
    const auto firstReferenceRegion = referenceRegions.front();
    const auto lastReferenceRegion = referenceRegions.back();

    const auto chromType = contigInfo.getChromType(firstReferenceRegion.contigIndex());
    const int32_t locusContigIndex = firstReferenceRegion.contigIndex();
    const int64_t locusAndFlanksStart = firstReferenceRegion.start() - kExtensionLength;
    const int64_t locusAndFlanksEnd = lastReferenceRegion.end() + kExtensionLength;
    const int64_t locusWithoutFlanksStart = firstReferenceRegion.start();
    const int64_t locusWithoutFlanksEnd = lastReferenceRegion.end();


    return LocusDescription(
        locusId,
        chromType,
        locusStructure,
        locusContigIndex,
        locusAndFlanksStart,
        locusAndFlanksEnd,
        locusWithoutFlanksStart,
        locusWithoutFlanksEnd,
        useRFC1MotifAnalysis,
        variantIds,
        referenceRegions,
        targetRegions,
        offtargetRegions,
        variantTypesFromUser,
        errorRate,
        likelihoodRatioThreshold,
        minimalLocusCoverage,
        plotConditions,
        extraFields
    );
}


void printMessageIfLocusArgNotFullyProcessed(
	vector<std::string>& locusIdsFromArg, vector<std::string>& processedLocusIds, const std::string& locusArg) {

    if (!locusIdsFromArg.empty() && processedLocusIds.size() < locusIdsFromArg.size()) {
        if (processedLocusIds.size() == 0) {
            if (locusIdsFromArg.size() > 1) {
                spdlog::error("ERROR: --locus arg: none of the locus ids {} were not found in the catalog", locusArg);
            } else {
                spdlog::error("ERROR: --locus arg: {} not found in the catalog", locusArg);
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


RegionCatalog convertLocusDescriptionsToLocusSpecs(
    LocusDescriptionCatalog& locusDescriptionCatalog, const HeuristicParameters& heuristicParams,
    Reference& reference)
{
    RegionCatalog catalog;
    catalog.reserve(locusDescriptionCatalog.size());
    for (auto& locusDescription : locusDescriptionCatalog)
    {
        try {
            LocusSpecification locusSpec = decodeLocusSpecification(locusDescription, reference, heuristicParams, true);
            catalog.emplace_back(std::move(locusSpec));
        }
        catch (const std::exception& e)
        {
            std::cout << "Error on locus spec " + locusDescription.locusId() + ": " + e.what() << std::endl;
        }
    }

    return catalog;
}


void sortAndFilterCatalog(
	LocusDescriptionCatalog& locusDescriptionCatalog, const ProgramParameters& programParams, const Reference& reference)

{
	switch (programParams.sortCatalogBy()) {
		case SortCatalogBy::kPosition: {
			spdlog::info("Sorting catalog by genomic coordinates");
			auto sortComparator = [](
				const LocusDescription& a, const LocusDescription& b
			) {
				return a.locusContigIndex() < b.locusContigIndex() || (
					a.locusContigIndex() == b.locusContigIndex() && (
						a.locusAndFlanksStart() < b.locusAndFlanksStart() || (
							a.locusAndFlanksStart() == b.locusAndFlanksStart() && a.locusAndFlanksEnd() < b.locusAndFlanksEnd())));
			};
			std::sort(locusDescriptionCatalog.begin(), locusDescriptionCatalog.end(), sortComparator);
			break;
		}
		case SortCatalogBy::kLocusId: {
			spdlog::info("Sorting catalog alphabetically by LocusId");
			auto sortComparator = [](
				const LocusDescription& a, const LocusDescription& b
			) {
				return a.locusId() < b.locusId();
			};
			std::sort(locusDescriptionCatalog.begin(), locusDescriptionCatalog.end(), sortComparator);
			break;
		}
		case SortCatalogBy::kNone:
			// No sorting needed
			break;
	}

    //if specified, apply the program.interval (which is a string like "chr1:1241251-64542152") to the locusDescriptionCatalog
    if (programParams.region().size() > 0) {
        const graphtools::ReferenceInterval interval = graphtools::ReferenceInterval::parseRegion(programParams.region());
        const auto intervalContigIndex = reference.contigInfo().getContigId(interval.contig);

        //apply filter
        // interval.contig //interval.start //interval.end
        locusDescriptionCatalog.erase(
            std::remove_if(
                locusDescriptionCatalog.begin(), locusDescriptionCatalog.end(),
                [&interval, &intervalContigIndex](const LocusDescription& locusDescription) {
                    return locusDescription.locusContigIndex() != intervalContigIndex ||
                        locusDescription.locusAndFlanksEnd() < interval.start ||
                        locusDescription.locusAndFlanksStart() > interval.end;
                }),
            locusDescriptionCatalog.end());
    }

    // discard the first userParams.startWith() loci
    if (programParams.startWith() > 0 ) {
        if (programParams.startWith() >= locusDescriptionCatalog.size()) {
            spdlog::warn("startWith value " + add_commas_at_thousands(programParams.startWith()) +
                " is greater than the number of loci (" + add_commas_at_thousands(locusDescriptionCatalog.size()) +
                "). Exiting...");
            return;
        }
        locusDescriptionCatalog.erase(locusDescriptionCatalog.begin(), locusDescriptionCatalog.begin() + programParams.startWith());
    }

    if (programParams.nLoci() > 0 && programParams.nLoci() < locusDescriptionCatalog.size()) {
        locusDescriptionCatalog.erase(locusDescriptionCatalog.begin() + programParams.nLoci(), locusDescriptionCatalog.end());
    }

    if (programParams.startWith() > 0 || programParams.nLoci() > 0) {
		spdlog::info("Kept {} loci after applying --start-with {} and -n {}",
		add_commas_at_thousands(locusDescriptionCatalog.size()), programParams.startWith(), programParams.nLoci());
	}
}


LocusDescriptionCatalog loadLocusDescriptions(const ProgramParameters& params, const Reference& reference)
{

	const std::string catalogPath = params.inputPaths().catalog();

    vector<std::string> locusIdsFromArg;
    if (!params.locus().empty()) {
        boost::split(locusIdsFromArg, params.locus(), boost::is_any_of(","));
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

    LocusDescriptionCatalog catalog;
    catalog.reserve(catalogJson.size());

    vector<std::string> processedLocusIds;

    const auto& contigInfo = reference.contigInfo();
    for (auto& locusJson : catalogJson)
    {

        if (!locusIdsFromArg.empty()) {
            assertFieldExists(locusJson, "LocusId");
            std::string currentLocusId = locusJson["LocusId"].get<string>();
            // check if currentLocusId is contained within the locusIdsFromArg vector
            if (std::find(locusIdsFromArg.begin(), locusIdsFromArg.end(), currentLocusId) == locusIdsFromArg.end()) {
                continue;
            }
            processedLocusIds.push_back(currentLocusId);
        }

        LocusDescription locusDescription = loadLocusDescription(locusJson, contigInfo, params.heuristics());
        catalog.emplace_back(std::move(locusDescription));
    }

    printMessageIfLocusArgNotFullyProcessed(locusIdsFromArg, processedLocusIds, params.locus());

	sortAndFilterCatalog(catalog, params, reference);

    return catalog;
}

}
