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

#include "io/VcfHeader.hh"
#include "io/IterativeVcfWriter.hh"
#include "core/ReadSupportCalculator.hh"
#include "io/VcfWriter.hh"
#include "io/VcfWriterHelpers.hh"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

namespace ehunter
{

namespace
{

std::string streamToString(const SmallVariantGenotype& genotype)
{
    std::ostringstream encoding;
    encoding << genotype;
    return encoding.str();
}

std::string streamToString(int64_t value)
{
    std::ostringstream encoding;
    encoding << value;
    return encoding.str();
}

}  // namespace

void IterativeVariantVcfWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto repeatNodeId = variantSpec_.nodes().front();
    const std::string& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    const int referenceSizeInUnits = referenceLocus.length() / repeatUnit.length();
    const std::string infoFields = computeInfoFields(variantSpec_, repeatUnit);

    const int posPreceedingRepeat1based = referenceLocus.start();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    const std::string leftFlankingBase
        = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

    const std::string altSymbol = computeAltSymbol(repeatFindingsPtr->optionalGenotype(), referenceSizeInUnits);
    const std::string alleleFields = computeAlleleFields(variantSpec_, repeatUnit, *repeatFindingsPtr);
    const std::string sampleFields = alleleFields + ":" + std::to_string(locusDepth_);

    std::string genotypeFilter = computeFilterSymbol(repeatFindingsPtr->genotypeFilter());

    vcfLine_ = {
        contigName,
        std::to_string(posPreceedingRepeat1based),
        ".",
        leftFlankingBase,
        altSymbol,
        ".",
        genotypeFilter,
        infoFields,
        "GT:SO:REPCN:REPCI:ADSP:ADFL:ADIR:LC",
        sampleFields
    };
}

void IterativeVariantVcfWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
    std::string refSequence;
    std::string altSequence;
    int64_t startPosition = -1;

    if (variantSpec_.classification().subtype == VariantSubtype::kSwap)
    {
        assert(variantSpec_.optionalRefNode());
        const auto refNode = *variantSpec_.optionalRefNode();
        const int refNodeIndex = refNode == variantSpec_.nodes().front() ? 0 : 1;
        const int altNodeIndex = refNode == variantSpec_.nodes().front() ? 1 : 0;

        const auto refNodeId = variantSpec_.nodes()[refNodeIndex];
        const auto altNodeId = variantSpec_.nodes()[altNodeIndex];

        refSequence = locusSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = locusSpec_.regionGraph().nodeSeq(altNodeId);
        startPosition = referenceLocus.start() + 1;
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kDeletion)
    {
        const std::string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int refNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase + locusSpec_.regionGraph().nodeSeq(refNodeId);
        altSequence = refFlankingBase;
        startPosition = referenceLocus.start();
    }
    else if (variantSpec_.classification().subtype == VariantSubtype::kInsertion)
    {
        const std::string refFlankingBase
            = reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

        const int altNodeId = variantSpec_.nodes().front();
        refSequence = refFlankingBase;
        altSequence = refFlankingBase + locusSpec_.regionGraph().nodeSeq(altNodeId);
        startPosition = referenceLocus.start();
    }
    else
    {
        std::ostringstream encoding;
        encoding << variantSpec_.classification().type << "/" << variantSpec_.classification().subtype;
        throw std::logic_error("Unable to generate VCF record for " + encoding.str());
    }

    const std::string infoFields = "VARID=" + variantSpec_.id();

    std::vector<std::string> sampleFields;
    std::vector<std::string> sampleValues;

    sampleFields.emplace_back("GT");
    const auto& optionalGenotype = smallVariantFindingsPtr->optionalGenotype();
    if (optionalGenotype)
    {
        sampleValues.push_back(streamToString(*optionalGenotype));
    }
    else
    {
        sampleValues.emplace_back(smallVariantFindingsPtr->alleleCount() == AlleleCount::kOne ? "." : "./.");
    }

    sampleFields.emplace_back("AD");
    std::ostringstream adEncoding;
    adEncoding << smallVariantFindingsPtr->numRefReads() << "," << smallVariantFindingsPtr->numAltReads();
    sampleValues.push_back(adEncoding.str());

    sampleFields.emplace_back("LC");
    sampleValues.push_back(std::to_string(locusDepth_));

    const std::string sampleField = boost::algorithm::join(sampleFields, ":");
    const std::string sampleValue = boost::algorithm::join(sampleValues, ":");

    std::string genotypeFilter = computeFilterSymbol(smallVariantFindingsPtr->genotypeFilter());

    vcfLine_ = {
        contigName,
        streamToString(startPosition),
        ".",
        refSequence,
        altSequence,
        ".",
        genotypeFilter,
        infoFields,
        sampleField,
        sampleValue
    };
}

IterativeVcfWriter::IterativeVcfWriter(
    std::string sampleId, Reference& reference, const std::string& outputFile)
    : sampleId_(std::move(sampleId)), reference_(reference)
{
    outFile_.open(outputFile, std::ios::out | std::ios::binary);
    if (!outFile_)
    {
        throw std::runtime_error("Failed to open VCF file: " + outputFile);
    }

    // Check if we need to compress with Gzip
    if (outputFile.size() > 2 && outputFile.substr(outputFile.size() - 2) == "gz")
    {
        outStream_.push(boost::iostreams::gzip_compressor());
    }

    outStream_.push(outFile_);

    // Write the VCF header upfront. Per-genotype `##ALT=<ID=STR{n}>` lines are intentionally omitted —
    // streaming mode does not know which STR<N> ALT symbols will be emitted in the body until all
    // records are written, so we don't try to predeclare them. Record bodies still emit `<STR{n}>`
    // ALT symbols; downstream consumers tolerate undeclared symbolic ALTs.
    outStream_ << "##fileformat=VCFv4.1\n";
    FieldDescriptionCatalog catalog;
    addCommonFieldDescriptions(catalog);
    addRepeatFieldDescriptions(catalog);
    addSmallVariantFieldDescriptions(catalog);
    for (const auto& fieldIdAndDescription : catalog)
    {
        outStream_ << fieldIdAndDescription.second << "\n";
    }
    writeBodyHeader(sampleId_, outStream_);
}

void IterativeVcfWriter::addRecord(const std::string& variantId, const LocusSpecification& locusSpec, const LocusFindings& locusFindings)
{
    const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);
    VariantFindings* findingsPtr = locusFindings.findingsForEachVariant.at(variantId).get();
    if (findingsPtr == nullptr)
    {
        return;
    }

    const double locusDepth = std::round(locusFindings.stats.depth() * 100) / 100.0;
    IterativeVariantVcfWriter recordWriter(reference_, locusSpec, locusDepth, variantSpec);
    findingsPtr->accept(&recordWriter);

    const std::vector<std::string>& vcfLine = recordWriter.getVcfLine();
    if (!vcfLine.empty())
    {
        outStream_ << boost::algorithm::join(vcfLine, "\t") << "\n";
    }
}

void IterativeVcfWriter::addRecords(const LocusSpecification& locusSpec, const LocusFindings& locusFindings)
{
    std::vector<std::string> variantIds;
    variantIds.reserve(locusFindings.findingsForEachVariant.size());
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        variantIds.push_back(variantIdAndFindings.first);
    }

    std::sort(
        variantIds.begin(), variantIds.end(),
        [&locusSpec](const std::string& a, const std::string& b) {
            const GenomicRegion& ra = locusSpec.getVariantSpecById(a).referenceLocus();
            const GenomicRegion& rb = locusSpec.getVariantSpecById(b).referenceLocus();
            if (ra.start() != rb.start())
            {
                return ra.start() < rb.start();
            }
            return ra.end() < rb.end();
        });

    for (const std::string& variantId : variantIds)
    {
        addRecord(variantId, locusSpec, locusFindings);
    }
}

void IterativeVcfWriter::close()
{
    if (closed_)
    {
        return;
    }
    closed_ = true;
    outStream_.flush();
    outStream_.reset();  // Ensure proper flushing of data (incl. gzip trailer if compressing)
    if (outFile_.is_open())
    {
        outFile_.close();
    }
}

IterativeVcfWriter::~IterativeVcfWriter()
{
    try
    {
        close();
    }
    catch (...)
    {
        // Never throw from a destructor.
    }
}

}
