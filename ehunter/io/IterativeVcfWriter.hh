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

#include <cstdint>
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

// The contigs that have at least one locus in `catalog`, which are the ones an IterativeVcfWriter's header
// declares. Every temp file of a run, and every run resuming it, must declare the same contigs, so this has
// to be computed from the whole catalog, before --resume drops the loci that are already done.
std::set<int32_t> contigsWithLoci(const LocusDescriptionCatalog& catalog);

// The exact header bytes an IterativeVcfWriter emits before its first record, declaring the contigs in
// `headerContigs`. --resume uses its size to sanity-check the sizes it cuts a temp file back to (see
// io/ResumeCheckpoint.hh).
std::string vcfDocumentHeader(
    const std::string& sampleId, const ReferenceContigInfo& contigInfo, const std::set<int32_t>& headerContigs);

// How a writer should attach to its output file.
enum class VcfOutputMode
{
    kTruncate,          // create/overwrite the file and write the VCF header
    kAppendAfterHeader, // open an existing partial VCF and continue after its last record line
};

// TODO: Document the code after multi-unit repeat format is finalized (GT-598)
class IterativeVcfWriter
{
public:
    // `headerContigs` are the contigs the header declares (see contigsWithLoci); unused in append mode.
    IterativeVcfWriter(
        std::string sampleId, Reference& reference, const std::set<int32_t>& headerContigs,
        const std::string& outputFilePath, VcfOutputMode outputMode = VcfOutputMode::kTruncate);
    // Ensure the VCF stream is flushed (and the gzip footer written if compressing) even when an
    // exception unwinds past the writer; otherwise a partially-written VCF may be left on disk.
    ~IterativeVcfWriter();

    void addRecord(const std::string& variantId, const LocusSpecification& locusSpec, const LocusFindings& locusFindings);
    // Emit every variant record for a locus in genomic-position order. findingsForEachVariant is an
    // unordered_map, so the variants are sorted by referenceLocus before output to keep the VCF position-sorted.
    void addRecords(const LocusSpecification& locusSpec, const LocusFindings& locusFindings);
    // Flush everything written so far through to the output file and return the file's size in bytes. See
    // IterativeJsonWriter::flushAndGetFileSize.
    std::uintmax_t flushAndGetFileSize();
    // Close the output file (idempotent). Throws if any of the writing failed, so a truncated output is
    // never mistaken for a complete one; the destructor swallows that, an explicit call propagates it.
    void close();

private:
    std::string sampleId_;
    Reference& reference_;
    std::string outputFilePath_;

    std::ofstream outFile_;
    boost::iostreams::filtering_ostream outStream_;
    bool closed_ = false;
};

}
