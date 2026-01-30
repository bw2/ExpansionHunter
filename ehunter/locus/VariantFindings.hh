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

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "core/CountTable.hh"
#include "genotyping/AlleleChecker.hh"
#include "genotyping/RepeatGenotype.hh"
#include "genotyping/SmallVariantGenotype.hh"
#include "locus/AlleleQualityMetrics.hh"
#include "locus/RFC1Status.hh"

namespace ehunter
{

enum class GenotypeFilter : unsigned
{
    kLowDepth = 1
};

GenotypeFilter operator|(GenotypeFilter left, GenotypeFilter right);
GenotypeFilter operator&(GenotypeFilter left, GenotypeFilter right);

class VariantFindings;
class RepeatFindings;
class SmallVariantFindings;

struct VariantFindingsVisitor
{
    virtual void visit(const RepeatFindings* findingsPtr) = 0;
    virtual void visit(const SmallVariantFindings* findingsPtr) = 0;
};

class VariantFindings
{
public:
    virtual ~VariantFindings() = default;
    virtual void accept(VariantFindingsVisitor* visitorPtr) = 0;
};

class RepeatFindings : public VariantFindings
{
public:
    RepeatFindings(
        CountTable countsOfSpanningReads, CountTable countsOfFlankingReads, CountTable countsOfInrepeatReads,
        AlleleCount alleleCount, boost::optional<RepeatGenotype> optionalGenotype, GenotypeFilter genotypeFilter)
        : countsOfSpanningReads_(std::move(countsOfSpanningReads))
        , countsOfFlankingReads_(std::move(countsOfFlankingReads))
        , countsOfInrepeatReads_(std::move(countsOfInrepeatReads))
        , alleleCount_(alleleCount)
        , optionalGenotype_(std::move(optionalGenotype))
        , genotypeFilter_(genotypeFilter)
    {
    }

    ~RepeatFindings() override = default;
    void accept(VariantFindingsVisitor* visitorPtr) override { visitorPtr->visit(this); }

    const CountTable& countsOfSpanningReads() const { return countsOfSpanningReads_; }
    const CountTable& countsOfFlankingReads() const { return countsOfFlankingReads_; }
    const CountTable& countsOfInrepeatReads() const { return countsOfInrepeatReads_; }

    AlleleCount alleleCount() const { return alleleCount_; }
    const boost::optional<RepeatGenotype>& optionalGenotype() const { return optionalGenotype_; }
    GenotypeFilter genotypeFilter() const { return genotypeFilter_; }

    void setRFC1Status(const RFC1Status& rfc1Status) { rfc1Status_ = rfc1Status; }

    boost::optional<RFC1Status> getRFC1Status() const { return rfc1Status_; }

    void setAlleleQualityMetrics(const RepeatAlleleQualityMetrics& metrics) { alleleQualityMetrics_ = metrics; }

    boost::optional<RepeatAlleleQualityMetrics> alleleQualityMetrics() const { return alleleQualityMetrics_; }

    void setConsensusSequences(const std::vector<std::string>& sequences) { consensusSequences_ = sequences; }

    const std::vector<std::string>& consensusSequences() const { return consensusSequences_; }

    void setConsensusReadSupport(const std::vector<std::string>& support) { consensusReadSupport_ = support; }

    const std::vector<std::string>& consensusReadSupport() const { return consensusReadSupport_; }

    void setCountsOfHighQualityUnambiguousReads(const CountTable& counts) { countsOfHighQualityUnambiguousReads_ = counts; }

    const CountTable& countsOfHighQualityUnambiguousReads() const { return countsOfHighQualityUnambiguousReads_; }

    // Compares core findings only: read counts and genotype.
    // Deliberately excludes:
    // - alleleCount_: derived from sample sex and chromosome type, not from the findings themselves
    // - genotypeFilter_: metadata about genotyping confidence, not the findings
    // - rfc1Status_, alleleQualityMetrics_: post-hoc enrichment fields
    // - consensusSequences_, consensusReadSupport_: derived consensus annotations
    // These excluded fields are either metadata, derived from external context, or
    // post-hoc annotations that don't represent the fundamental repeat findings.
    bool operator==(const RepeatFindings& other) const
    {
        return countsOfSpanningReads_ == other.countsOfSpanningReads_
            && countsOfFlankingReads_ == other.countsOfFlankingReads_
            && countsOfInrepeatReads_ == other.countsOfInrepeatReads_
            && optionalGenotype_ == other.optionalGenotype_;
    }

private:
    CountTable countsOfSpanningReads_;
    CountTable countsOfFlankingReads_;
    CountTable countsOfInrepeatReads_;
    AlleleCount alleleCount_;
    boost::optional<RepeatGenotype> optionalGenotype_;
    GenotypeFilter genotypeFilter_;
    boost::optional<RFC1Status> rfc1Status_;
    boost::optional<RepeatAlleleQualityMetrics> alleleQualityMetrics_;
    std::vector<std::string> consensusSequences_;   // One consensus sequence per allele
    std::vector<std::string> consensusReadSupport_; // Per-position read support as digit string, one per allele
    CountTable countsOfHighQualityUnambiguousReads_; // Counts by allele size (repeat units)
};

class SmallVariantFindings : public VariantFindings
{
public:
    SmallVariantFindings(
        int numRefReads, int numAltReads, AlleleCheckSummary refAlleleStatus, AlleleCheckSummary altAlleleStatus,
        AlleleCount alleleCount, boost::optional<SmallVariantGenotype> optionalGenotype, GenotypeFilter genotypeFilter)
        : numRefReads_(numRefReads)
        , numAltReads_(numAltReads)
        , refAlleleStatus_(refAlleleStatus)
        , altAlleleStatus_(altAlleleStatus)
        , alleleCount_(alleleCount)
        , optionalGenotype_(std::move(optionalGenotype))
        , genotypeFilter_(genotypeFilter)
    {
    }

    ~SmallVariantFindings() override = default;
    void accept(VariantFindingsVisitor* visitorPtr) override { visitorPtr->visit(this); }

    int numRefReads() const { return numRefReads_; }
    int numAltReads() const { return numAltReads_; }
    AlleleCount alleleCount() const { return alleleCount_; }
    const boost::optional<SmallVariantGenotype>& optionalGenotype() const { return optionalGenotype_; }
    GenotypeFilter genotypeFilter() const { return genotypeFilter_; }

    AlleleCheckSummary refAllelePresenceStatus() const { return refAlleleStatus_; }
    AlleleCheckSummary altAllelePresenceStatus() const { return altAlleleStatus_; }

private:
    int numRefReads_;
    int numAltReads_;
    AlleleCheckSummary refAlleleStatus_;
    AlleleCheckSummary altAlleleStatus_;
    AlleleCount alleleCount_;
    boost::optional<SmallVariantGenotype> optionalGenotype_;
    GenotypeFilter genotypeFilter_;
};

std::ostream& operator<<(std::ostream& out, const RepeatFindings& repeatFindings);

}
