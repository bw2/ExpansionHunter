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

#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"

namespace ehunter
{

class VariantVcfWriter : public VariantFindingsVisitor
{
public:
    VariantVcfWriter(
        Reference& reference, const LocusSpecification& locusSpec, double locusDepth,
        const VariantSpecification& variantSpec, std::ostream& out)
        : reference_(reference)
        , locusSpec_(locusSpec)
        , locusDepth_(locusDepth)
        , variantSpec_(variantSpec)
        , out_(out)
    {
    }

    ~VariantVcfWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;

private:
    Reference& reference_;
    const LocusSpecification& locusSpec_;
    double locusDepth_;
    const VariantSpecification& variantSpec_;
    std::ostream& out_;
};

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class VcfWriter
{
public:
    VcfWriter(
        std::string sampleId, Reference& reference, const RegionCatalog& regionCatalog,
        const SampleFindings& sampleFindings);

    friend std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);

private:
    void writeHeader(std::ostream& out);
    void writeBody(std::ostream& out);
    using LocusIndexAndVariantId = std::pair<unsigned, std::string>;
    const std::vector<LocusIndexAndVariantId> getSortedIdPairs();

    std::string sampleId_;
    Reference& reference_;
    const RegionCatalog& regionCatalog_;
    const SampleFindings& sampleFindings_;
};

std::ostream& operator<<(std::ostream& out, VcfWriter& vcfWriter);
std::string computeFilterSymbol(GenotypeFilter filter);
void writeBodyHeader(const std::string& sampleName, std::ostream& out);
std::string createRepeatAlleleSymbol(int repeatSize);
std::string computeAltSymbol(const boost::optional<RepeatGenotype>& optionalGenotype, int referenceSizeInUnits);
std::string computeInfoFields(const VariantSpecification& variantSpec, const std::string& repeatUnit);
std::string computeAlleleFields(const VariantSpecification& variantSpec, const std::string& repeatUnit, const RepeatFindings& repeatFindings);

// Build the tab-separated VCF record fields for a single variant. Shared by the seeking-mode
// VariantVcfWriter and the streaming-mode IterativeVariantVcfWriter so both modes emit identical
// records from one source of truth.
std::vector<std::string> buildRepeatVcfRecordElements(
    Reference& reference, const LocusSpecification& locusSpec, double locusDepth,
    const VariantSpecification& variantSpec, const RepeatFindings& repeatFindings);
std::vector<std::string> buildSmallVariantVcfRecordElements(
    Reference& reference, const LocusSpecification& locusSpec, double locusDepth,
    const VariantSpecification& variantSpec, const SmallVariantFindings& smallVariantFindings);

}
