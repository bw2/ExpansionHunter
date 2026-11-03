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

#pragma once

#include <limits>

#include "core/Parameters.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"

#include "thirdparty/json/json.hpp"

namespace ehunter
{

namespace gq
{
struct GenotypeQualityModel;
}

class VariantJsonWriter : public VariantFindingsVisitor
{
public:
    // When qualityModel is non-null, per-allele genotype-quality fields
    // (PredictedLengthCorrectionFactor / pOk / pTooShort / pTooLong) are added to each
    // AlleleQualityMetrics allele record, except for a RareRepeat at a locus that specifies OfftargetRegions.
    // PredictedLengthCorrectionFactor is written as 1 where the written pOk is 0.5 or more.
    // A null model leaves output unchanged.
    // locusCoverage is the model's `coverage` feature and must be the same value the caller
    // emits as the locus record's Coverage field; NaN when it is not available.
    VariantJsonWriter(
        const ReferenceContigInfo& contigInfo, const LocusSpecification& locusSpec,
        const VariantSpecification& variantSpec, const gq::GenotypeQualityModel* qualityModel = nullptr,
        double locusCoverage = std::numeric_limits<double>::quiet_NaN())
        : contigInfo_(contigInfo)
        , locusSpec_(locusSpec)
        , variantSpec_(variantSpec)
        , qualityModel_(qualityModel)
        , locusCoverage_(locusCoverage)
    {
    }

    ~VariantJsonWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;
    nlohmann::json record() const { return record_; }

private:
    const ReferenceContigInfo& contigInfo_;
    const LocusSpecification& locusSpec_;
    const VariantSpecification& variantSpec_;
    const gq::GenotypeQualityModel* qualityModel_;
    double locusCoverage_;
    nlohmann::json record_;
};

class JsonWriter
{
public:
    JsonWriter(
        const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings, bool copyCatalogFields = false,
        const gq::GenotypeQualityModel* qualityModel = nullptr, std::time_t startedEpoch = 0, int threadCount = 1,
        AnalysisMode analysisMode = AnalysisMode::kSeeking, const std::string& commandLine = "");

    void write(std::ostream& out);

private:
    const SampleParameters& sampleParams_;
    const ReferenceContigInfo& contigInfo_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
    bool copyCatalogFields_;
    const gq::GenotypeQualityModel* qualityModel_;
    std::time_t startedEpoch_;
    int threadCount_;
    AnalysisMode analysisMode_;
    std::string commandLine_;
};

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter);

// Locus-level "RepeatAllelePhasing" record, shared by JsonWriter and IterativeJsonWriter
nlohmann::json encodeRepeatAllelePhasing(const RepeatAllelePhasing& phasing);

// Every output JSON holds these top-level records in this order: SampleParameters, RunInfo, LocusResults
// (left out by JsonWriter when there are no loci) and RunStats. RunInfo, which says how ExpansionHunter was run,
// comes before the loci so a reader can get it without reading the whole file. RunStats holds what is only
// known once the run ends, so it comes last.

// The RunInfo record. GenotypingApproach is written only in optimized-streaming mode.
nlohmann::json buildRunInfoRecord(
    AnalysisMode analysisMode, OptimizedStreamingGenotypingApproach genotypingApproach, int threadCount,
    const std::string& commandLine, const gq::GenotypeQualityModel* qualityModel);

// The RunStats record: when the run started and finished, how long it took, and its peak memory. Completed and
// Runtime are left out when completedNormally is false (a writer closed while an exception unwinds), so a crashed
// run doesn't look complete.
nlohmann::json buildRunStatsRecord(std::time_t startedEpoch, bool completedNormally);

// The text of a record, with two-space indentation and every line after the first indented by two more spaces,
// ready to be the value of a top-level field.
std::string formatAsTopLevelValue(const nlohmann::json& record);

}
