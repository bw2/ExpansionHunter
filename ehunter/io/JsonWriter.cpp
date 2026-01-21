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

#include "io/JsonWriter.hh"

#include <iomanip>
#include <sstream>
#include <vector>

#include <boost/algorithm/string/join.hpp>
#include <boost/optional.hpp>

#include "core/Common.hh"
#include "core/ReadSupportCalculator.hh"

namespace ehunter
{

using std::map;
using std::string;
using Json = nlohmann::json;
using boost::optional;
using std::to_string;
using std::vector;

// Round double to 3 decimal places for cleaner JSON output
static double round3(double value)
{
    return std::round(value * 1000.0) / 1000.0;
}

std::ostream& operator<<(std::ostream& out, JsonWriter& jsonWriter)
{
    jsonWriter.write(out);
    return out;
}

JsonWriter::JsonWriter(
    const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo, const RegionCatalog& regionCatalog,
    const SampleFindings& sampleFindings, bool copyCatalogFields)
    : sampleParams_(sampleParams)
    , contigInfo_(contigInfo)
    , regionCatalog_(regionCatalog)
    , sampleFindings_(sampleFindings)
    , copyCatalogFields_(copyCatalogFields)
{
}

void JsonWriter::write(std::ostream& out)
{
    Json sampleParametersRecord;
    sampleParametersRecord["SampleId"] = sampleParams_.id();
    sampleParametersRecord["Sex"] = streamToString(sampleParams_.sex());

    Json resultsRecord;
    const unsigned locusCount(sampleFindings_.size());
    for (unsigned locusIndex(0); locusIndex < locusCount; ++locusIndex)
    {
        const LocusSpecification& locusSpec = regionCatalog_[locusIndex];
        const LocusFindings& locusFindings = sampleFindings_[locusIndex];
        const std::string& locusId(locusSpec.locusId());

        Json locusRecord;

        // Copy extra annotation fields from input catalog first (if enabled),
        // so computed values take precedence in case of field name collisions
        if (copyCatalogFields_ && locusSpec.extraFields().has_value())
        {
            const nlohmann::json& extraFields = locusSpec.extraFields().value();
            for (auto it = extraFields.begin(); it != extraFields.end(); ++it)
            {
                locusRecord[it.key()] = it.value();
            }
        }

        locusRecord["LocusId"] = locusId;
        locusRecord["Coverage"] = std::round(locusFindings.stats.depth() * 100) / 100.0;
        locusRecord["ReadLength"] = locusFindings.stats.meanReadLength();
        locusRecord["FragmentLength"] = locusFindings.stats.meanFragLength();
        locusRecord["AlleleCount"] = static_cast<int>(locusFindings.stats.alleleCount());

        Json variantRecords;
        for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
        {
            const string& variantId = variantIdAndFindings.first;
            const VariantSpecification& variantSpec = locusSpec.getVariantSpecById(variantId);

            VariantJsonWriter variantWriter(contigInfo_, locusSpec, variantSpec);
            variantIdAndFindings.second->accept(&variantWriter);
            variantRecords[variantId] = variantWriter.record();
        }

        if (!variantRecords.empty())
        {
            locusRecord["Variants"] = variantRecords;
        }
        resultsRecord[locusId] = locusRecord;
    }

    Json sampleRecords;
    if (!resultsRecord.empty())
    {
        sampleRecords["LocusResults"] = resultsRecord;
    }
    sampleRecords["SampleParameters"] = sampleParametersRecord;

    out << std::setw(2) << sampleRecords << std::endl;
}

static string encodeGenotype(const RepeatGenotype& genotype)
{
    string encoding = std::to_string(genotype.shortAlleleSizeInUnits());

    if (genotype.numAlleles() == 2)
    {
        encoding = encoding + "/" + std::to_string(genotype.longAlleleSizeInUnits());
    }

    return encoding;
}

void VariantJsonWriter::visit(const RepeatFindings* repeatFindingsPtr)
{
    assert(variantSpec_.classification().type == VariantType::kRepeat);

    const RepeatFindings& repeatFindings = *repeatFindingsPtr;

    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);

    const auto repeatNodeId = variantSpec_.nodes().front();
    const auto& repeatUnit = locusSpec_.regionGraph().nodeSeq(repeatNodeId);
    record_["RepeatUnit"] = repeatUnit;

    record_["CountsOfSpanningReads"] = streamToString(repeatFindings.countsOfSpanningReads());
    record_["CountsOfFlankingReads"] = streamToString(repeatFindings.countsOfFlankingReads());
    record_["CountsOfInrepeatReads"] = streamToString(repeatFindings.countsOfInrepeatReads());

    if (repeatFindings.optionalGenotype())
    {
        record_["Genotype"] = encodeGenotype(*repeatFindings.optionalGenotype());
        record_["GenotypeConfidenceInterval"] = streamToString(*repeatFindings.optionalGenotype());
    }

    const auto rfc1Status(repeatFindings.getRFC1Status());
    if (rfc1Status)
    {
        nlohmann::json rfc1Results;
        rfc1Results["Call"] = label(rfc1Status->call);
        rfc1Results["Description"] = rfc1Status->description;
        record_["RFC1MotifAnalysis"] = rfc1Results;
    }

    // Output AlleleQualityMetrics if present
    const auto alleleQualityMetrics = repeatFindings.alleleQualityMetrics();
    if (alleleQualityMetrics && alleleQualityMetrics->hasMetrics)
    {
        nlohmann::json metricsRecord;
        metricsRecord["VariantId"] = alleleQualityMetrics->variantId;

        nlohmann::json allelesArray = nlohmann::json::array();
        for (const auto& allele : alleleQualityMetrics->alleles)
        {
            nlohmann::json alleleRecord;
            alleleRecord["AlleleNumber"] = allele.alleleNumber;
            alleleRecord["AlleleSize"] = allele.alleleSize;
            alleleRecord["Depth"] = round3(allele.depth);
            alleleRecord["QD"] = round3(allele.qd);
            alleleRecord["MeanInsertedBasesWithinRepeats"] = round3(allele.meanInsertedBasesWithinRepeats);
            alleleRecord["MeanDeletedBasesWithinRepeats"] = round3(allele.meanDeletedBasesWithinRepeats);
            alleleRecord["StrandBiasBinomialPhred"] = round3(allele.strandBiasBinomialPhred);
            alleleRecord["LeftFlankNormalizedDepth"] = round3(allele.leftFlankNormalizedDepth);
            alleleRecord["RightFlankNormalizedDepth"] = round3(allele.rightFlankNormalizedDepth);
            alleleRecord["HighQualityUnambiguousReads"] = allele.highQualityUnambiguousReads;
            alleleRecord["ConfidenceIntervalDividedByAlleleSize"] = round3(allele.confidenceIntervalDividedByAlleleSize);

            allelesArray.push_back(alleleRecord);
        }
        metricsRecord["Alleles"] = allelesArray;
        record_["AlleleQualityMetrics"] = metricsRecord;
    }

    // Output ConsensusSequences and related statistics if present
    const auto& consensusSequences = repeatFindings.consensusSequences();
    if (!consensusSequences.empty())
    {
        nlohmann::json consensusArray = nlohmann::json::array();
        for (const auto& seq : consensusSequences)
        {
            consensusArray.push_back(seq);
        }
        record_["ConsensusSequences"] = consensusArray;

        // Output ConsensusSequencesReadSupport if available
        const auto& consensusReadSupport = repeatFindings.consensusReadSupport();
        if (!consensusReadSupport.empty())
        {
            nlohmann::json supportArray = nlohmann::json::array();
            for (const auto& support : consensusReadSupport)
            {
                supportArray.push_back(support);
            }
            record_["ConsensusSequencesReadSupport"] = supportArray;
        }
    }
}

void VariantJsonWriter::visit(const SmallVariantFindings* smallVariantFindingsPtr)
{
    const SmallVariantFindings& findings = *smallVariantFindingsPtr;
    record_.clear();
    record_["VariantId"] = variantSpec_.id();
    record_["VariantType"] = streamToString(variantSpec_.classification().type);
    record_["VariantSubtype"] = streamToString(variantSpec_.classification().subtype);
    record_["ReferenceRegion"] = encode(contigInfo_, variantSpec_.referenceLocus());
    record_["CountOfRefReads"] = findings.numRefReads();
    record_["CountOfAltReads"] = findings.numAltReads();
    record_["StatusOfRefAllele"] = streamToString(findings.refAllelePresenceStatus().status);
    record_["LogLikelihoodRefAllelePresent"] = streamToString(findings.refAllelePresenceStatus().logLikelihoodRatio);
    record_["StatusOfAltAllele"] = streamToString(findings.altAllelePresenceStatus().status);
    record_["LogLikelihoodAltAllelePresent"] = streamToString(findings.altAllelePresenceStatus().logLikelihoodRatio);
    if (findings.optionalGenotype())
    {
        record_["Genotype"] = streamToString(*findings.optionalGenotype());
    }
}

}
