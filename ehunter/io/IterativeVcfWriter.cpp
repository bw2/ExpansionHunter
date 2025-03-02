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

namespace ehunter
{

IterativeVcfWriter::IterativeVcfWriter(
    std::string sampleId, Reference& reference, const std::string& outputFile)
    : sampleId_(std::move(sampleId)), reference_(reference), headerWritten_(false)
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
}



void IterativeVcfWriter::addRecord(const std::string& variantId, const LocusSpecification& locusSpec, const LocusFindings& locusFindings)
{

    const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

    const auto repeatFindingsPtr = dynamic_cast<RepeatFindings*>(locusFindings.findingsForEachVariant.at(variantId).get());

	//record info for VCF header
	//const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

	if (!headerWritten_)
	{
		outStream_ << "##fileformat=VCFv4.1\n";
		FieldDescriptionWriter descriptionWriter(locusSpec, variantSpec);
		repeatFindingsPtr->accept(&descriptionWriter);
		descriptionWriter.dumpTo(fieldDescriptionCatalog_);
		for (const auto& fieldIdAndDescription : fieldDescriptionCatalog_)
		{
			const auto& description = fieldIdAndDescription.second;
			outStream_ << description << "\n";
		}

		writeBodyHeader(sampleId_, outStream_);
		headerWritten_ = true;
	}

	//generate VCF row
	const double locusDepth = locusFindings.stats.depth();

	const auto& referenceLocus = variantSpec.referenceLocus();
	const auto repeatNodeId = variantSpec.nodes().front();
	const std::string& repeatUnit = locusSpec.regionGraph().nodeSeq(repeatNodeId);
	const int referenceSizeInUnits = referenceLocus.length() / repeatUnit.length();
	const std::string infoFields = computeInfoFields(variantSpec, repeatUnit);

	const int posPreceedingRepeat1based = referenceLocus.start();
	const auto& contigName = reference_.contigInfo().getContigName(referenceLocus.contigIndex());
	const std::string leftFlankingBase
		= reference_.getSequence(contigName, referenceLocus.start() - 1, referenceLocus.start());

	const std::string altSymbol = computeAltSymbol(repeatFindingsPtr->optionalGenotype(), referenceSizeInUnits);
	const std::string alleleFields = computeAlleleFields(variantSpec, repeatUnit, *repeatFindingsPtr);
	const std::string sampleFields = alleleFields + ":" + std::to_string(locusDepth);

	std::string genotypeFilter = computeFilterSymbol(repeatFindingsPtr->genotypeFilter());

	std::vector<std::string> vcfLine_ = {
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

	outStream_ << boost::algorithm::join(vcfLine_, "\t") << "\n";
}


void IterativeVcfWriter::close()
{
    outStream_.flush();
    outStream_.reset();  // Ensure proper flushing of data
    if (outFile_.is_open())
    {
        outFile_.close();
    }
}

}

