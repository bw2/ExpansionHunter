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
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <boost/algorithm/string/join.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "io/VcfHeader.hh"

namespace ehunter
{

class IterativeVariantVcfWriter : public VariantFindingsVisitor
{
public:
    IterativeVariantVcfWriter(
        Reference& reference, const LocusSpecification& locusSpec, double locusDepth,
        const VariantSpecification& variantSpec)
        : reference_(reference)
        , locusSpec_(locusSpec)
        , locusDepth_(locusDepth)
        , variantSpec_(variantSpec)
    {
    }

    ~IterativeVariantVcfWriter() = default;
    void visit(const RepeatFindings* repeatFindingsPtr) override;
    void visit(const SmallVariantFindings* smallVariantFindingsPtr) override;
    std::vector<std::string> getVcfLine() const { return vcfLine_; }

private:
    Reference& reference_;
    const LocusSpecification& locusSpec_;
    double locusDepth_;
    const VariantSpecification& variantSpec_;
    std::vector<std::string> vcfLine_;
};

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class IterativeVcfWriter
{
public:
    IterativeVcfWriter(
        std::string sampleId, Reference& reference, const std::string& outputFilePath);

    void addRecord(const std::string& variantId, const LocusSpecification& locusSpec, const LocusFindings& locusFindings);
    void close();  // Close output file

private:
    void writeVcfLine(const std::vector<std::string>& vcfLine);

    std::string sampleId_;
    Reference& reference_;
    FieldDescriptionCatalog fieldDescriptionCatalog_;

    std::ofstream outFile_;
    boost::iostreams::filtering_ostream outStream_;
    bool headerWritten_;
};

}
