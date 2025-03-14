//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
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

#include "locus/LocusSpecification.hh"

#include <algorithm>
#include <cassert>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "spdlog/spdlog.h"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/Reference.hh"

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
LocusSpecification::LocusSpecification(
    RegionId locusId, ChromType typeOfChromLocusLocatedOn, std::vector<GenomicRegion> targetReadExtractionRegions,
    graphtools::Graph regionGraph, NodeToRegionAssociation referenceRegions, GenotyperParameters genotyperParams,
    const bool useRFC1MotifAnalysis)
    : locusId_(std::move(locusId))
    , typeOfChromLocusLocatedOn_(typeOfChromLocusLocatedOn)
    , targetReadExtractionRegions_(std::move(targetReadExtractionRegions))
    , regionGraph_(std::move(regionGraph))
    , referenceRegions_(std::move(referenceRegions))
    , parameters_(std::move(genotyperParams))
    , useRFC1MotifAnalysis_(useRFC1MotifAnalysis)
{
}

void LocusSpecification::addVariantSpecification(
    std::string id, VariantClassification classification, GenomicRegion referenceLocus, vector<NodeId> nodes,
    optional<NodeId> refNode)
{
    variantSpecs_.emplace_back(std::move(id), classification, std::move(referenceLocus), std::move(nodes), refNode);
}

const VariantSpecification& LocusSpecification::getVariantSpecById(const std::string& variantSpecId) const
{
    for (const auto& variantSpec : variantSpecs_)
    {
        if (variantSpec.id() == variantSpecId)
        {
            return variantSpec;
        }
    }

    throw std::logic_error("There is no variant " + variantSpecId + " in locus " + locusId_);
}

bool LocusSpecification::requiresGenomeWideDepth() const
{
    for (const auto& variantSpec : variantSpecs_)
    {
        if (variantSpec.classification().subtype == VariantSubtype::kSMN)
        {
            return true;
        }
    }

    return false;
}

std::ostream& operator<<(std::ostream& out, const LocusDescription& locusDescription)
{
	out << "[" << locusDescription.locusId() << ": " << (locusDescription.locusContigIndex()+1) << ":" << locusDescription.locusAndFlanksStart() << "-" << locusDescription.locusAndFlanksEnd() << "]";

	return out;
}

LocusDescription::LocusDescription(
    std::string locusId, ChromType chromType, std::string locusStructure,
    int32_t locusContigIndex,
    int64_t locusAndFlanksStart, int64_t locusAndFlanksEnd,
    int64_t locusWithoutFlanksStart, int64_t locusWithoutFlanksEnd,
    bool useRFC1MotifAnalysis,
    std::vector<std::string> variantIds, std::vector<GenomicRegion> referenceRegions,
    std::vector<GenomicRegion> targetRegions, std::vector<GenomicRegion> offtargetRegions,
    std::vector<VariantTypeFromUser> variantTypesFromUser,
    std::optional<double> errorRate, std::optional<double> likelihoodRatioThreshold,
    std::optional<double> minLocusCoverage)
    : locusId_(std::move(locusId))
    , chromType_(chromType)
    , locusStructure_(std::move(locusStructure))
    , locusContigIndex_(locusContigIndex)
    , locusAndFlanksStart_(locusAndFlanksStart)
    , locusAndFlanksEnd_(locusAndFlanksEnd)
    , locusWithoutFlanksStart_(locusWithoutFlanksStart)
    , locusWithoutFlanksEnd_(locusWithoutFlanksEnd)
    , useRFC1MotifAnalysis_(useRFC1MotifAnalysis)
    , variantIds_(std::move(variantIds))
    , referenceRegions_(std::move(referenceRegions))
    , targetRegions_(std::move(targetRegions))
    , offtargetRegions_(std::move(offtargetRegions))
    , variantTypesFromUser_(std::move(variantTypesFromUser))
    , errorRate_(errorRate)
    , likelihoodRatioThreshold_(likelihoodRatioThreshold)
    , minLocusCoverage_(minLocusCoverage)
{
}

}
