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
    // (PredictedLengthCorrectionFactor / pTooShort / pTooLong) are added to each
    // AlleleQualityMetrics allele record. A null model leaves output unchanged.
    VariantJsonWriter(
        const ReferenceContigInfo& contigInfo, const LocusSpecification& locusSpec,
        const VariantSpecification& variantSpec, const gq::GenotypeQualityModel* qualityModel = nullptr)
        : contigInfo_(contigInfo)
        , locusSpec_(locusSpec)
        , variantSpec_(variantSpec)
        , qualityModel_(qualityModel)
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

}
