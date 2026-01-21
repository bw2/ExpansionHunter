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

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <optional>

#include <boost/optional.hpp>

#include "graphcore/Graph.hh"
#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/GenomicRegion.hh"
#include "core/Reference.hh"
#include "locus/VariantSpecification.hh"

namespace ehunter
{

// Plot policy determines whether read visualizations are generated for a locus
enum class PlotPolicy {
    kAll,           // Generate read visualizations for all loci (--plot-all)
    kNone,          // Don't generate any read visualizations (--disable-all-plots)
    kConditional    // Generate based on plotConditions_ (default)
};

// Enums for plot conditions
enum class PlotThresholdAppliedTo { kShortAllele, kLongAllele };
enum class PlotThresholdComparisonOp {
    kLessThan,
    kLessThanOrEqual,
    kEqual,
    kNotEqual,
    kGreaterThanOrEqual,
    kGreaterThan
};

// Single plot condition for conditional image generation
struct PlotReadVisualization {
    PlotThresholdAppliedTo appliedTo;    // Which allele to compare
    PlotThresholdComparisonOp op;        // Comparison operator
    int threshold;                        // Value to compare against
};

using RegionId = std::string;
using NodeToRegionAssociation = std::unordered_map<graphtools::NodeId, GenomicRegion>;

class LocusSpecification
{
public:
    LocusSpecification(
        RegionId locusId, ChromType typeOfChromLocusLocatedOn, std::vector<GenomicRegion> targetReadExtractionRegions,
        graphtools::Graph regionGraph, NodeToRegionAssociation referenceRegions, GenotyperParameters genotyperParams,
        bool useRFC1MotifAnalysis,
        std::vector<PlotReadVisualization> plotConditions = {},
        std::optional<nlohmann::json> extraFields = std::nullopt);

    const RegionId& locusId() const { return locusId_; }
    ChromType typeOfChromLocusLocatedOn() const { return typeOfChromLocusLocatedOn_; }
    /*
     * List of all regions in the reference this graph describes
     * i.e. where we expect relevant reads to align
     */
    const std::vector<GenomicRegion>& targetReadExtractionRegions() const { return targetReadExtractionRegions_; }
    /*
     * List of regions that additional relevant reads might be found
     * Require filtering or special considerations
     */
    const std::vector<GenomicRegion>& offtargetReadExtractionRegions() const { return offtargetReadExtractionRegions_; }
    void setOfftargetReadExtractionRegions(const std::vector<GenomicRegion>& offtargetReadExtractionRegions)
    {
        offtargetReadExtractionRegions_ = offtargetReadExtractionRegions;
    }

    const graphtools::Graph& regionGraph() const { return regionGraph_; }
    const GenotyperParameters& genotyperParameters() const { return parameters_; }
    const std::vector<VariantSpecification>& variantSpecs() const { return variantSpecs_; }

    void addVariantSpecification(
        std::string id, VariantClassification classification, GenomicRegion referenceLocus,
        std::vector<graphtools::NodeId> nodes, boost::optional<graphtools::NodeId> optionalRefNode);

    const VariantSpecification& getVariantSpecById(const std::string& variantSpecId) const;

    const NodeToRegionAssociation& referenceProjectionOfNodes() const { return referenceRegions_; }

    bool requiresGenomeWideDepth() const;

    bool useRFC1MotifAnalysis() const { return useRFC1MotifAnalysis_; }

    const std::vector<PlotReadVisualization>& plotConditions() const { return plotConditions_; }
    bool hasPlotConditions() const { return !plotConditions_.empty(); }
    const std::optional<nlohmann::json>& extraFields() const { return extraFields_; }

    PlotPolicy plotPolicy() const { return plotPolicy_; }
    void setPlotPolicy(PlotPolicy policy) { plotPolicy_ = policy; }

    /// Returns true if this locus requires an alignment buffer for post-genotyping analysis
    /// (RFC1 motif analysis, plot-all policy, or conditional image generation)
    bool requiresAlignmentBuffer() const {
        return useRFC1MotifAnalysis_ || plotPolicy_ == PlotPolicy::kAll || hasPlotConditions();
    }

private:
    std::string locusId_;
    ChromType typeOfChromLocusLocatedOn_;
    std::vector<GenomicRegion> targetReadExtractionRegions_;
    std::vector<GenomicRegion> offtargetReadExtractionRegions_;
    graphtools::Graph regionGraph_;
    std::vector<VariantSpecification> variantSpecs_;
    NodeToRegionAssociation referenceRegions_;
    GenotyperParameters parameters_;
    bool useRFC1MotifAnalysis_;
    std::vector<PlotReadVisualization> plotConditions_;
    PlotPolicy plotPolicy_ = PlotPolicy::kConditional;
    std::optional<nlohmann::json> extraFields_;
};

enum class VariantTypeFromUser
{
    kRareRepeat,
    kCommonRepeat,
    kSmallVariant
};

class LocusDescription
{
public:
    LocusDescription(
        std::string locusId, ChromType chromType, std::string locusStructure,
        int32_t locusContigIndex,
        int64_t locusAndFlanksStart, int64_t locusAndFlanksEnd,
        int64_t locusWithoutFlanksStart, int64_t locusWithoutFlanksEnd,
        bool useRFC1MotifAnalysis = false,
        std::vector<std::string> variantIds = {},
        std::vector<GenomicRegion> referenceRegions = {},
        std::vector<GenomicRegion> targetRegions = {},
        std::vector<GenomicRegion> offtargetRegions = {},
        std::vector<VariantTypeFromUser> variantTypesFromUser = {},
        std::optional<double> errorRate = std::nullopt,
        std::optional<double> likelihoodRatioThreshold = std::nullopt,
        std::optional<double> minLocusCoverage = std::nullopt,
        std::vector<PlotReadVisualization> plotConditions = {},
        std::optional<nlohmann::json> extraFields = std::nullopt
    );

    // Getters
    const std::string& locusId() const { return locusId_; }
    ChromType chromType() const { return chromType_; }
    const std::string& locusStructure() const { return locusStructure_; }
    int32_t locusContigIndex() const { return locusContigIndex_; }
    int64_t locusAndFlanksStart() const { return locusAndFlanksStart_; }
    int64_t locusAndFlanksEnd() const { return locusAndFlanksEnd_; }
    int64_t locusWithoutFlanksStart() const { return locusWithoutFlanksStart_; }
    int64_t locusWithoutFlanksEnd() const { return locusWithoutFlanksEnd_; }
    bool useRFC1MotifAnalysis() const { return useRFC1MotifAnalysis_; }
    const std::vector<std::string>& variantIds() const { return variantIds_; }
    const std::vector<GenomicRegion>& referenceRegions() const { return referenceRegions_; }
    const std::vector<GenomicRegion>& targetRegions() const { return targetRegions_; }
    const std::vector<GenomicRegion>& offtargetRegions() const { return offtargetRegions_; }
    const std::vector<VariantTypeFromUser>& variantTypesFromUser() const { return variantTypesFromUser_; }
    const std::optional<double>& errorRate() const { return errorRate_; }
    const std::optional<double>& likelihoodRatioThreshold() const { return likelihoodRatioThreshold_; }
    const std::optional<double>& minLocusCoverage() const { return minLocusCoverage_; }
    const std::vector<PlotReadVisualization>& plotConditions() const { return plotConditions_; }
    bool hasPlotConditions() const { return !plotConditions_.empty(); }
    const std::optional<nlohmann::json>& extraFields() const { return extraFields_; }

private:
    std::string locusId_;
    ChromType chromType_;
    std::string locusStructure_;
    int32_t locusContigIndex_;
    int64_t locusAndFlanksStart_;
    int64_t locusAndFlanksEnd_;
    int64_t locusWithoutFlanksStart_;
    int64_t locusWithoutFlanksEnd_;
    bool useRFC1MotifAnalysis_;
    std::vector<std::string> variantIds_;
    std::vector<GenomicRegion> referenceRegions_;
    std::vector<GenomicRegion> targetRegions_;
    std::vector<GenomicRegion> offtargetRegions_;
    std::vector<VariantTypeFromUser> variantTypesFromUser_;
    std::optional<double> errorRate_;
    std::optional<double> likelihoodRatioThreshold_;
    std::optional<double> minLocusCoverage_;
    std::vector<PlotReadVisualization> plotConditions_;
    std::optional<nlohmann::json> extraFields_;
};

std::ostream& operator<<(std::ostream& out, const LocusDescription& locusDescription);

using RegionCatalog = std::vector<LocusSpecification>;

using LocusDescriptionCatalog = std::vector<LocusDescription>;

}
