//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
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

#include "io/ResumeCheckpoint.hh"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "core/Common.hh"
#include "core/Parameters.hh"
#include "thirdparty/json/json.hpp"
#include "locus/LocusSpecification.hh"

using namespace ehunter;

namespace
{

const char* const kTestPrefix = "ResumeTest_tmp";

ResumeCheckpointPaths testPaths()
{
    return ResumeCheckpointPaths{ std::string(kTestPrefix) + ".json.unfinished",
        std::string(kTestPrefix) + ".vcf.unfinished", std::string(kTestPrefix) + ".processed_loci.unfinished" };
}

ResumeCheckpointPaths testPathsGz()
{
    return ResumeCheckpointPaths{ std::string(kTestPrefix) + ".json.gz.unfinished",
        std::string(kTestPrefix) + ".vcf.gz.unfinished",
        std::string(kTestPrefix) + ".processed_loci.unfinished" };
}

// A locus on `contigIndex` whose one variant is named after it, matching what the catalog loader produces
// for a single-region locus.
LocusDescription makeLocus(const std::string& locusId, std::int32_t contigIndex, std::int64_t start)
{
    return LocusDescription(locusId, ChromType::kAutosome, "(CAG)*", contigIndex, start - 1000, start + 1000,
        start, start + 30, false, { locusId });
}

// The four loci used by most tests: two on contig 0, one each on contigs 1 and 2.
LocusDescriptionCatalog makeCatalog()
{
    return { makeLocus("L1", 0, 1000), makeLocus("L2", 0, 5000), makeLocus("L3", 1, 2000),
        makeLocus("L4", 2, 3000) };
}

// The record text IterativeJsonWriter would have captured for `locusId`.
std::string jsonRecordFor(const std::string& locusId)
{
    return "\n    \"" + locusId + "\": {\n      \"LocusId\": \"" + locusId + "\",\n      \"Coverage\": 30.0\n    }";
}

// The VCF line IterativeVcfWriter would have captured for `locusId`'s single variant.
std::string vcfLineFor(const std::string& locusId)
{
    return "chr\t100\t.\tA\t.\t.\tPASS\tEND=130;VARID=" + locusId + ";REPID=" + locusId + "\tGT\t0/0\n";
}

void removeIfPresent(const std::string& path) { std::remove(path.c_str()); }

std::string readFile(const std::string& path)
{
    std::ifstream file(path, std::ios::in | std::ios::binary);
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

bool fileExists(const std::string& path)
{
    std::ifstream file(path);
    return static_cast<bool>(file);
}

class ResumeTest : public ::testing::Test
{
protected:
    void SetUp() override { cleanUp(); }
    void TearDown() override { cleanUp(); }

    void cleanUp()
    {
        for (const ResumeCheckpointPaths& paths : { testPaths(), testPathsGz() })
        {
            removeIfPresent(paths.json);
            removeIfPresent(paths.vcf);
            removeIfPresent(paths.processedLoci);
        }
        for (std::int32_t contigIndex = 0; contigIndex < 4; ++contigIndex)
        {
            removeIfPresent(contigTempJsonPath(kTestPrefix, contigIndex));
            removeIfPresent(contigTempVcfPath(kTestPrefix, contigIndex));
        }
        removeIfPresent(singleSliceTempJsonPath(kTestPrefix));
        removeIfPresent(singleSliceTempVcfPath(kTestPrefix));
    }

    // Writes a checkpoint holding `locusIds` in order, as an interrupted run would have left it.
    void writeCheckpoint(const ResumeCheckpointPaths& paths, const std::vector<std::string>& locusIds)
    {
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", false, false, 0);
        for (const std::string& locusId : locusIds)
        {
            writer.recordLocus(locusId, jsonRecordFor(locusId), vcfLineFor(locusId));
        }
        writer.close();
    }

    SampleParameters sampleParams_{ "sample", Sex::kFemale };
};

TEST_F(ResumeTest, CheckpointPathsFollowTheOutputPathsCompression)
{
    const OutputPaths plain("out.vcf", "out.json", "out.bam", "out.tsv", "out");
    EXPECT_EQ(resumeCheckpointPaths(plain).json, "out.json.unfinished");
    EXPECT_EQ(resumeCheckpointPaths(plain).vcf, "out.vcf.unfinished");
    EXPECT_EQ(resumeCheckpointPaths(plain).processedLoci, "out.processed_loci.unfinished");

    const OutputPaths compressed("out.vcf.gz", "out.json.gz", "out.bam", "out.tsv", "out");
    EXPECT_EQ(resumeCheckpointPaths(compressed).json, "out.json.gz.unfinished");
    EXPECT_EQ(resumeCheckpointPaths(compressed).vcf, "out.vcf.gz.unfinished");
}

TEST_F(ResumeTest, AllFinishedLociAreRecovered)
{
    const ResumeCheckpointPaths paths = testPaths();
    writeCheckpoint(paths, { "L1", "L2", "L3" });
    ASSERT_TRUE(resumeCheckpointExists(paths));

    const LocusDescriptionCatalog catalog = makeCatalog();
    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, catalog, kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.size(), 3u);
    EXPECT_EQ(loaded.doneLocusIds.count("L4"), 0u);
    EXPECT_TRUE(loaded.hasExistingJsonRecords);
    // Contigs 0 (L1, L2) and 1 (L3) were rebuilt; contig 2 has nothing finished.
    EXPECT_EQ(loaded.rebuiltSlices.size(), 2u);
    EXPECT_EQ(loaded.rebuiltSlices.count(0), 1u);
    EXPECT_EQ(loaded.rebuiltSlices.count(1), 1u);
    EXPECT_EQ(loaded.rebuiltSlices.count(2), 0u);
}

TEST_F(ResumeTest, RebuiltTempsHoldTheirContigsRecordsInOrder)
{
    writeCheckpoint(testPaths(), { "L1", "L2", "L3" });
    loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true);

    const std::string contig0Json = readFile(contigTempJsonPath(kTestPrefix, 0));
    EXPECT_NE(contig0Json.find("\"LocusResults\": {"), std::string::npos);
    EXPECT_LT(contig0Json.find("\"L1\""), contig0Json.find("\"L2\""));
    EXPECT_EQ(contig0Json.find("\"L3\""), std::string::npos);
    // The temp is left open-ended for the genotyping writer to append to and close.
    EXPECT_EQ(contig0Json.find("\"RunInfo\""), std::string::npos);

    const std::string contig1Json = readFile(contigTempJsonPath(kTestPrefix, 1));
    EXPECT_NE(contig1Json.find("\"L3\""), std::string::npos);
    EXPECT_EQ(contig1Json.find("\"L1\""), std::string::npos);

    const std::string contig0Vcf = readFile(contigTempVcfPath(kTestPrefix, 0));
    EXPECT_NE(contig0Vcf.find("#CHROM"), std::string::npos);
    EXPECT_NE(contig0Vcf.find("VARID=L1"), std::string::npos);
    EXPECT_NE(contig0Vcf.find("VARID=L2"), std::string::npos);
    EXPECT_EQ(contig0Vcf.find("VARID=L3"), std::string::npos);

    EXPECT_FALSE(fileExists(contigTempJsonPath(kTestPrefix, 2)));
}

TEST_F(ResumeTest, ARecordWrittenAfterTheLastProcessedLocusIsDropped)
{
    const ResumeCheckpointPaths paths = testPaths();
    writeCheckpoint(paths, { "L1", "L2" });

    // The run got further in the record files than in the processed-loci list, which is what an
    // interruption between the two writes leaves behind.
    std::ofstream(paths.json, std::ios::app) << ", " << jsonRecordFor("L3");
    std::ofstream(paths.vcf, std::ios::app) << vcfLineFor("L3");

    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.count("L3"), 0u);
    EXPECT_EQ(readFile(paths.json).find("\"L3\""), std::string::npos);
    EXPECT_EQ(readFile(paths.vcf).find("VARID=L3"), std::string::npos);
    EXPECT_FALSE(fileExists(contigTempJsonPath(kTestPrefix, 1)));
}

TEST_F(ResumeTest, ATruncatedTrailingRecordIsDropped)
{
    const ResumeCheckpointPaths paths = testPaths();
    writeCheckpoint(paths, { "L1", "L2", "L3" });

    // Cut the JSON file mid-record, as a lost buffered write would. The processed-loci list still names
    // every locus, so the loss has to be detected rather than silently dropping L3 from the output.
    const std::string content = readFile(paths.json);
    const std::size_t lastRecordStart = content.rfind("\"L3\"");
    ASSERT_NE(lastRecordStart, std::string::npos);
    std::ofstream(paths.json, std::ios::trunc | std::ios::binary) << content.substr(0, lastRecordStart + 10);

    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.size(), 2u);
    EXPECT_EQ(loaded.doneLocusIds.count("L3"), 0u);
    EXPECT_EQ(readFile(paths.processedLoci).find("L3"), std::string::npos);
}

TEST_F(ResumeTest, ALocusThatProducedNoRecordIsStillTreatedAsFinished)
{
    const ResumeCheckpointPaths paths = testPaths();
    {
        // L2 is what --skip-hom-ref leaves behind: genotyped, but with no record in either output file.
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.recordLocus("L2", "", "");
        writer.close();
    }

    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.count("L2"), 1u);
    EXPECT_EQ(readFile(paths.json).find("\"L2\""), std::string::npos);
}

TEST_F(ResumeTest, GzippedCheckpointsRoundTripThroughSeveralAppends)
{
    const ResumeCheckpointPaths paths = testPathsGz();
    {
        // Each close() of a batch ends a gzip member, so several appends make this a multi-member file.
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.close();
    }
    {
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", true, true, 0);
        writer.recordLocus("L2", jsonRecordFor("L2"), vcfLineFor("L2"));
        writer.close();
    }

    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.size(), 2u);
    EXPECT_EQ(loaded.doneLocusIds.count("L1"), 1u);
    EXPECT_EQ(loaded.doneLocusIds.count("L2"), 1u);
    EXPECT_NE(readFile(contigTempJsonPath(kTestPrefix, 0)).find("\"L2\""), std::string::npos);
}

TEST_F(ResumeTest, PerContigSlicesKeepRecordsThatFinishedOutOfCatalogOrder)
{
    // A locus is released for genotyping once the reads pass the end of its flank window, so a locus nested
    // inside a wider one finishes first even though it comes second in the position-sorted catalog. The
    // checkpoint records the order the run actually emitted in, and appending to a rebuilt temp only needs
    // that order, so nothing here may be discarded.
    writeCheckpoint(testPaths(), { "L2", "L1" });

    const ResumeLoadResult loaded
        = loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.size(), 2u);
    // The rebuilt temp has to preserve the emission order, not re-sort into catalog order: the genotyping
    // writer appends after these, and an uninterrupted run would have emitted them in this same order.
    const std::string contig0Json = readFile(contigTempJsonPath(kTestPrefix, 0));
    EXPECT_LT(contig0Json.find("\"L2\""), contig0Json.find("\"L1\""));
}

TEST_F(ResumeTest, SingleSliceKeepsRecordsFromAnAscendingContigSweep)
{
    // One slice covers every contig in one coordinate sweep, which is what a --threads 1 checkpoint is:
    // its contigs only ever move forwards, whatever the order within each of them.
    writeCheckpoint(testPaths(), { "L2", "L1", "L3" });

    const ResumeLoadResult loaded
        = loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, false);

    EXPECT_EQ(loaded.doneLocusIds.size(), 3u);
    EXPECT_EQ(loaded.rebuiltSlices.count(kSingleSliceKey), 1u);
}

TEST_F(ResumeTest, SingleSliceRefusesRecordsThatSkipAnUnfinishedContig)
{
    // A --threads > 1 checkpoint finishes contigs independently, so it can hold contig 1's locus while
    // contig 0 is still running. One coordinate sweep never produces that, and cannot append to it either:
    // it would emit contig 0's loci after records already sitting in the temp for contig 1.
    writeCheckpoint(testPaths(), { "L3", "L1", "L2" });  // contigs 1, 0, 0

    EXPECT_THROW(
        loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, false), std::runtime_error);
}

TEST_F(ResumeTest, RefusingASingleSliceResumeLeavesTheCheckpointIntact)
{
    // The refusal exists so the operator can re-run with --threads > 1 and keep everything, which only
    // works if nothing was rewritten on the way out.
    writeCheckpoint(testPaths(), { "L3", "L1", "L2" });  // contig 1 finished before contig 0
    const std::string checkpointBefore = readFile(testPaths().processedLoci);
    const std::string recordsBefore = readFile(testPaths().json);

    EXPECT_THROW(
        loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, false), std::runtime_error);

    EXPECT_EQ(readFile(testPaths().processedLoci), checkpointBefore);
    EXPECT_EQ(readFile(testPaths().json), recordsBefore);

    // The same checkpoint is still fully usable per contig, which is what the error tells the operator.
    const ResumeLoadResult loaded
        = loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true);
    EXPECT_EQ(loaded.doneLocusIds.size(), 3u);
}

TEST_F(ResumeTest, SingleSliceKeepsASweepThatFinishesEachContigBeforeMovingOn)
{
    // L1 and L2 are all of contig 0, so the sweep may move on to contig 1's L3; every entry is reusable.
    writeCheckpoint(testPaths(), { "L1", "L2", "L3" });

    const ResumeLoadResult loaded
        = loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, false);

    EXPECT_EQ(loaded.doneLocusIds.size(), 3u);
}

TEST_F(ResumeTest, ALocusMissingFromTheCatalogIsRejected)
{
    writeCheckpoint(testPaths(), { "L1", "NOT_IN_CATALOG" });

    EXPECT_THROW(
        loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true), std::runtime_error);
}

TEST_F(ResumeTest, AMissingLocusResultsMarkerIsRejected)
{
    writeCheckpoint(testPaths(), { "L1" });
    std::ofstream(testPaths().json, std::ios::trunc) << "not an ExpansionHunter document";

    EXPECT_THROW(
        loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true), std::runtime_error);
}

TEST_F(ResumeTest, RunSignatureMismatchesAreDescribed)
{
    EXPECT_EQ(describeSignatureMismatch(R"({"Catalog":"a","Threads":2})", R"({"Catalog":"a","Threads":2})"), "");

    const std::string changedValue
        = describeSignatureMismatch(R"({"Catalog":"a"})", R"({"Catalog":"b"})");
    EXPECT_NE(changedValue.find("Catalog"), std::string::npos);
    EXPECT_NE(changedValue.find("\"b\""), std::string::npos);

    EXPECT_NE(describeSignatureMismatch(R"({"A":1,"B":2})", R"({"A":1})").find("B"), std::string::npos);
    EXPECT_NE(describeSignatureMismatch(R"({"A":1})", R"({"A":1,"B":2})").find("B"), std::string::npos);
    EXPECT_FALSE(describeSignatureMismatch("not json", R"({"A":1})").empty());
}

TEST_F(ResumeTest, RunSignatureIsStoredInAndReadBackFromTheCheckpoint)
{
    const std::string signature = R"({
  "Catalog": "catalog.json",
  "SampleId": "sample"
})";
    {
        ResumeCheckpointWriter writer(testPaths(), sampleParams_, signature, false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.close();
    }

    EXPECT_EQ(describeSignatureMismatch(signature, readCheckpointRunSignature(testPaths())), "");
}

TEST_F(ResumeTest, ACatalogWithDuplicateLocusIdsIsRejected)
{
    LocusDescriptionCatalog catalog = makeCatalog();
    EXPECT_NO_THROW(assertCatalogIsResumable(catalog));

    // Resume keys on LocusId alone, so finishing one of two entries sharing an id would mark both done.
    catalog.push_back(makeLocus("L2", 3, 9000));
    EXPECT_THROW(assertCatalogIsResumable(catalog), std::runtime_error);
}

TEST_F(ResumeTest, ACatalogWhoseVariantIdsCollideAcrossLociIsRejected)
{
    // VCF lines are matched back to their locus by VariantId, so a shared one would file a locus's lines
    // under a different locus.
    LocusDescriptionCatalog catalog = makeCatalog();
    catalog.push_back(LocusDescription("L5", ChromType::kAutosome, "(CAG)*", 3, 8000, 10000, 9000, 9030,
        false, { "L1" }));

    EXPECT_THROW(assertCatalogIsResumable(catalog), std::runtime_error);
}

TEST_F(ResumeTest, TheOppositeCompressionSettingsPathsAreDerivable)
{
    const OutputPaths plain("out.vcf", "out.json", "out.bam", "out.tsv", "out");
    EXPECT_EQ(otherCompressionOutputPaths(plain).json(), "out.json.gz");
    EXPECT_EQ(otherCompressionOutputPaths(plain).vcf(), "out.vcf.gz");

    const OutputPaths compressed("out.vcf.gz", "out.json.gz", "out.bam", "out.tsv", "out");
    EXPECT_EQ(otherCompressionOutputPaths(compressed).json(), "out.json");
    EXPECT_EQ(otherCompressionOutputPaths(compressed).vcf(), "out.vcf");
}

TEST_F(ResumeTest, ATruncatedGzipTailKeepsTheMembersBeforeIt)
{
    // A gzip checkpoint is a chain of members, one per flush. Losing the tail of the last one must not
    // cost the complete members before it, or a single lost byte would discard the whole run's progress.
    const ResumeCheckpointPaths paths = testPathsGz();
    {
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.close();
    }
    {
        ResumeCheckpointWriter writer(paths, sampleParams_, "{}", true, true, 0);
        writer.recordLocus("L2", jsonRecordFor("L2"), vcfLineFor("L2"));
        writer.close();
    }

    const std::string content = readFile(paths.json);
    std::ofstream(paths.json, std::ios::trunc | std::ios::binary)
        << content.substr(0, content.size() - 1);  // lose one byte of the final member

    const ResumeLoadResult loaded = loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true);

    EXPECT_EQ(loaded.doneLocusIds.count("L1"), 1u);
}

TEST_F(ResumeTest, CorruptedGzipDataIsRejectedRatherThanPartlyTrusted)
{
    // Damaged bytes are different from a cut-short tail: whatever inflate produces before it notices is
    // not trustworthy, so the checkpoint has to be reported unreadable instead of replayed.
    const ResumeCheckpointPaths paths = testPathsGz();
    writeCheckpoint(paths, { "L1", "L2" });

    std::string content = readFile(paths.json);
    ASSERT_GT(content.size(), 60u);
    content[content.size() / 2] ^= 0xFF;  // flip a byte inside a member's compressed data
    std::ofstream(paths.json, std::ios::trunc | std::ios::binary) << content;

    EXPECT_THROW(loadResumeCheckpoint(paths, sampleParams_, makeCatalog(), kTestPrefix, true),
        ResumeCheckpointCorruptError);
}

TEST_F(ResumeTest, AnOversizedRunSignatureIsStillReadable)
{
    // The signature quotes options verbatim, so a long --locus list makes it far larger than any fixed
    // read window. Failing to read it back would discard the whole checkpoint on every resume.
    nlohmann::json signature;
    signature["Catalog"] = "catalog.json";
    signature["Locus"] = std::string(200 * 1024, 'x');
    const std::string signatureText = signature.dump(2);

    {
        ResumeCheckpointWriter writer(testPaths(), sampleParams_, signatureText, false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.close();
    }

    EXPECT_EQ(describeSignatureMismatch(signatureText, readCheckpointRunSignature(testPaths())), "");
}

TEST_F(ResumeTest, RewritingACheckpointLeavesItsSignatureUnchanged)
{
    // loadResumeCheckpoint rewrites the header using the signature it just read, so anything that is not
    // idempotent there (re-indenting it, say) would compound on every successive resume.
    const std::string signatureText = R"({
  "Catalog": "catalog.json",
  "SampleId": "sample"
})";
    {
        ResumeCheckpointWriter writer(testPaths(), sampleParams_, signatureText, false, false, 0);
        writer.recordLocus("L1", jsonRecordFor("L1"), vcfLineFor("L1"));
        writer.close();
    }

    const std::string afterFirstWrite = readCheckpointRunSignature(testPaths());
    for (int resumeCount = 0; resumeCount != 3; ++resumeCount)
    {
        loadResumeCheckpoint(testPaths(), sampleParams_, makeCatalog(), kTestPrefix, true);
        EXPECT_EQ(readCheckpointRunSignature(testPaths()), afterFirstWrite);
    }
}

TEST_F(ResumeTest, CompletedLociAreRemovedFromTheCatalog)
{
    LocusDescriptionCatalog catalog = makeCatalog();
    removeCompletedLoci(catalog, { "L2", "L4" });

    ASSERT_EQ(catalog.size(), 2u);
    EXPECT_EQ(catalog[0].locusId(), "L1");
    EXPECT_EQ(catalog[1].locusId(), "L3");
}

TEST_F(ResumeTest, RemovingCompletedLociFromAnEmptyFinishedSetChangesNothing)
{
    LocusDescriptionCatalog catalog = makeCatalog();
    removeCompletedLoci(catalog, {});
    EXPECT_EQ(catalog.size(), 4u);
}

}  // namespace
