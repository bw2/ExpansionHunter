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

#include <memory>

using std::ostream;
using std::string;

namespace ehunter
{

namespace
{

// Insert a description into the catalog if no entry with the same (type, id) is already present.
void tryAddFieldDescription(FieldDescriptionCatalog& catalog, FieldType type, const string& id, const string& number,
    const string& contentType, const string& description)
{
    FieldDescriptionIdentifier key{ type, id };
    if (catalog.find(key) == catalog.end())
    {
        catalog.emplace(key, FieldDescription(type, id, number, contentType, description));
    }
}

}  // namespace

void addCommonFieldDescriptions(FieldDescriptionCatalog& catalog)
{
    tryAddFieldDescription(
        catalog, FieldType::kInfo, "VARID", "1", "String", "Variant identifier as specified in the variant catalog");
    tryAddFieldDescription(catalog, FieldType::kFormat, "GT", "1", "String", "Genotype");
    tryAddFieldDescription(catalog, FieldType::kFormat, "LC", "1", "Float", "Locus coverage");
    tryAddFieldDescription(catalog, FieldType::kFilter, "PASS", "", "", "All filters passed");
    tryAddFieldDescription(catalog, FieldType::kFilter, "LowDepth", "", "",
        "The overall locus depth is below 10x or number of reads spanning one or both breakends is below 5");
}

void addRepeatFieldDescriptions(FieldDescriptionCatalog& catalog)
{
    tryAddFieldDescription(catalog, FieldType::kInfo, "SVTYPE", "1", "String", "Type of structural variant");
    tryAddFieldDescription(catalog, FieldType::kInfo, "END", "1", "Integer", "End position of the variant");
    tryAddFieldDescription(catalog, FieldType::kInfo, "REF", "1", "Integer", "Reference copy number");
    tryAddFieldDescription(catalog, FieldType::kInfo, "RL", "1", "Integer", "Reference length in bp");
    tryAddFieldDescription(catalog, FieldType::kInfo, "RU", "1", "String", "Repeat unit in the reference orientation");
    tryAddFieldDescription(
        catalog, FieldType::kInfo, "REPID", "1", "String", "Repeat identifier as specified in the variant catalog");
    tryAddFieldDescription(catalog, FieldType::kFormat, "SO", "1", "String",
        "Type of reads that support the allele; can be SPANNING, FLANKING, or INREPEAT meaning that the reads span, "
        "flank, or are fully contained in the repeat");
    tryAddFieldDescription(
        catalog, FieldType::kFormat, "REPCN", "1", "String", "Number of repeat units spanned by the allele");
    tryAddFieldDescription(catalog, FieldType::kFormat, "REPCI", "1", "String", "Confidence interval for REPCN");
    tryAddFieldDescription(
        catalog, FieldType::kFormat, "ADFL", "1", "String", "Number of flanking reads consistent with the allele");
    tryAddFieldDescription(
        catalog, FieldType::kFormat, "ADSP", "1", "String", "Number of spanning reads consistent with the allele");
    tryAddFieldDescription(
        catalog, FieldType::kFormat, "ADIR", "1", "String", "Number of in-repeat reads consistent with the allele");
}

void addSmallVariantFieldDescriptions(FieldDescriptionCatalog& catalog)
{
    tryAddFieldDescription(catalog, FieldType::kFormat, "AD", ".", "Integer",
        "Allelic depths for the ref and alt alleles in the order listed");
}

void FieldDescriptionWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    addCommonFieldDescriptions(fieldDescriptions_);
    addRepeatFieldDescriptions(fieldDescriptions_);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const string& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    const auto& referenceLocus = variantSpec_.referenceLocus();
    const int referenceSize = referenceLocus.length() / repeatUnit.length();

    if (repeatFindingsPtr->optionalGenotype()) {
		const RepeatGenotype& genotype = repeatFindingsPtr->optionalGenotype().get();

		if (genotype.shortAlleleSizeInUnits() != referenceSize)
		{
			const string sizeEncoding = std::to_string(genotype.shortAlleleSizeInUnits());
			const string description = "Allele comprised of " + sizeEncoding + " repeat units";
			tryAddingFieldDescription(FieldType::kAlt, "STR" + sizeEncoding, "", "", description);
		}

		if (genotype.longAlleleSizeInUnits() != referenceSize)
		{
			const string sizeEncoding = std::to_string(genotype.longAlleleSizeInUnits());
			const string description = "Allele comprised of " + sizeEncoding + " repeat units";
			tryAddingFieldDescription(FieldType::kAlt, "STR" + sizeEncoding, "", "", description);
		}
	}
}

void FieldDescriptionWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    if (!smallVariantFindingsPtr->optionalGenotype())
    {
        return;
    }
    addCommonFieldDescriptions(fieldDescriptions_);
    addSmallVariantFieldDescriptions(fieldDescriptions_);
}

void FieldDescriptionWriter::tryAddingFieldDescription(
    FieldType fieldType, const string& id, const string& number, const string& contentType, const string& description)
{
    const auto key = std::make_pair(fieldType, id);
    if (fieldDescriptions_.find(key) == fieldDescriptions_.end())
    {
        FieldDescription fieldDescription(fieldType, id, number, contentType, description);
        fieldDescriptions_.emplace(std::make_pair(key, std::move(fieldDescription)));
    }
}

void FieldDescriptionWriter::dumpTo(FieldDescriptionCatalog& descriptionCatalog)
{
    descriptionCatalog.insert(fieldDescriptions_.begin(), fieldDescriptions_.end());
}

FieldDescription::FieldDescription(
    FieldType fieldType, string id, string number, string contentType, string description)
    : fieldType(fieldType)
    , id(std::move(id))
    , number(std::move(number))
    , contentType(std::move(contentType))
    , description(std::move(description))
{
}

void outputVcfHeader(const RegionCatalog& locusCatalog, const SampleFindings& sampleFindings, ostream& out)
{
    out << "##fileformat=VCFv4.1\n";

    FieldDescriptionCatalog fieldDescriptionCatalog;

    const unsigned locusCount(sampleFindings.size());
    for (unsigned locusIndex(0); locusIndex < locusCount; ++locusIndex)
    {
        const LocusSpecification& locusSpec = locusCatalog[locusIndex];
        const LocusFindings& locusFindings = sampleFindings[locusIndex];

        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

            FieldDescriptionWriter descriptionWriter(locusSpec, variantSpec);
            variantIdAndFindings.second->accept(&descriptionWriter);
            descriptionWriter.dumpTo(fieldDescriptionCatalog);
        }
    }

    for (const auto& fieldIdAndDescription : fieldDescriptionCatalog)
    {
        const auto& description = fieldIdAndDescription.second;
        out << description << "\n";
    }
}

std::ostream& operator<<(std::ostream& out, FieldType fieldType)
{
    switch (fieldType)
    {
    case FieldType::kAlt:
        out << "ALT";
        break;
    case FieldType::kFormat:
        out << "FORMAT";
        break;
    case FieldType::kInfo:
        out << "INFO";
        break;
    case FieldType::kFilter:
        out << "FILTER";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, const FieldDescription& fieldDescription)
{
    switch (fieldDescription.fieldType)
    {
    case FieldType::kInfo:
    case FieldType::kFormat:
        out << "##" << fieldDescription.fieldType << "=<ID=" << fieldDescription.id
            << ",Number=" << fieldDescription.number << ",Type=" << fieldDescription.contentType << ",Description=\""
            << fieldDescription.description << "\">";
        break;
    case FieldType::kAlt:
    case FieldType::kFilter:
        out << "##" << fieldDescription.fieldType << "=<ID=" << fieldDescription.id << ",Description=\""
            << fieldDescription.description << "\">";
        break;
    }

    return out;
}

}
