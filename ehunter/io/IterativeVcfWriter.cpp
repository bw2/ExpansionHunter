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
#include <cerrno>
#include <cstring>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>

namespace ehunter
{

void IterativeVariantVcfWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    vcfLine_ = buildRepeatVcfRecordElements(reference_, locusSpec_, locusDepth_, variantSpec_, *repeatFindingsPtr);
}

void IterativeVariantVcfWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    vcfLine_ = buildSmallVariantVcfRecordElements(
        reference_, locusSpec_, locusDepth_, variantSpec_, *smallVariantFindingsPtr);
}

std::set<int32_t> contigsWithLoci(const LocusDescriptionCatalog& catalog)
{
    // A locus's variants all lie on its contig: decoding a locus fails unless its reference regions merge
    // into one region (mergeRegions in io/LocusSpecDecoding.cpp).
    std::set<int32_t> contigIndices;
    for (const LocusDescription& locusDescription : catalog)
    {
        contigIndices.insert(locusDescription.locusContigIndex());
    }
    return contigIndices;
}

std::string vcfDocumentHeader(
    const std::string& sampleId, const ReferenceContigInfo& contigInfo, const std::set<int32_t>& headerContigs)
{
    // Per-genotype `##ALT=<ID=STR{n}>` lines are intentionally omitted — streaming mode does not know
    // which STR<N> ALT symbols will be emitted in the body until all records are written, so we don't try
    // to predeclare them. Record bodies still emit `<STR{n}>` ALT symbols; downstream consumers tolerate
    // undeclared symbolic ALTs.
    std::ostringstream header;
    header << "##fileformat=VCFv4.1\n";
    FieldDescriptionCatalog catalog;
    addCommonFieldDescriptions(catalog);
    addRepeatFieldDescriptions(catalog);
    addSmallVariantFieldDescriptions(catalog);
    for (const auto& fieldIdAndDescription : catalog)
    {
        header << fieldIdAndDescription.second << "\n";
    }
    outputVcfContigLines(contigInfo, headerContigs, header);
    writeBodyHeader(sampleId, header);
    return header.str();
}

IterativeVcfWriter::IterativeVcfWriter(
    std::string sampleId, Reference& reference, const std::set<int32_t>& headerContigs, const std::string& outputFile,
    VcfOutputMode outputMode)
    : sampleId_(std::move(sampleId)), reference_(reference), outputFilePath_(outputFile)
{
    const bool append = outputMode == VcfOutputMode::kAppendAfterHeader;
    const std::ios::openmode openMode
        = std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc);
    outFile_.open(outputFile, openMode);
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

    // In append mode the header is whatever the interrupted run left behind, so only a fresh file
    // writes one.
    if (!append)
    {
        outStream_ << vcfDocumentHeader(sampleId_, reference_.contigInfo(), headerContigs);
    }
}

void IterativeVcfWriter::addRecord(const std::string& variantId, const LocusSpecification& locusSpec, const LocusFindings& locusFindings)
{
    const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);
    VariantFindings* findingsPtr = locusFindings.findingsForEachVariant.at(variantId).get();
    if (findingsPtr == nullptr)
    {
        return;
    }

    const double locusDepth = locusFindings.stats.depth();
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

std::uintmax_t IterativeVcfWriter::flushAndGetFileSize()
{
    outStream_.flush();
    outFile_.flush();
    if (!outStream_ || !outFile_)
    {
        throw std::runtime_error("Failed to write " + outputFilePath_ + " (" + std::strerror(errno) + ")");
    }
    return boost::filesystem::file_size(outputFilePath_);
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
    outFile_.flush();
    const bool failed = !outFile_;
    if (outFile_.is_open())
    {
        outFile_.close();
    }
    // As in IterativeJsonWriter::close: a silently truncated VCF would be merged as if complete.
    if (failed || !outFile_)
    {
        throw std::runtime_error("Failed to write " + outputFilePath_ + " (" + std::strerror(errno) + ")");
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
