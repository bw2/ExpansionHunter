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

#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "gtest/gtest.h"

#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "core/ReferenceContigInfo.hh"
#include "io/IterativeJsonWriter.hh"
#include "io/IterativeVcfWriter.hh"
#include "locus/LocusSpecification.hh"
#include "thirdparty/json/json.hpp"

using namespace ehunter;

namespace
{

const char* const kTestPrefix = "ResumeTest_tmp";
const char* const kSignature = R"({"Catalog":"catalog.json","SampleId":"sample"})";

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

// A JSON record shaped like the ones IterativeJsonWriter writes for `locusId`.
std::string jsonRecordFor(const std::string& locusId)
{
    return "\n    \"" + locusId + "\": {\n      \"LocusId\": \"" + locusId + "\",\n      \"Coverage\": 30.0\n    }";
}

// A VCF line shaped like the one IterativeVcfWriter writes for `locusId`'s single variant.
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

void appendToFile(const std::string& path, const std::string& text)
{
    std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::app);
    file << text;
}

std::vector<std::string> readLines(const std::string& path)
{
    std::ifstream file(path);
    std::vector<std::string> lines;
    for (std::string line; std::getline(file, line);)
    {
        lines.push_back(line);
    }
    return lines;
}

// The writers need a Reference to format VCF records, which these tests never do, and for the names and
// lengths of the contigs the VCF header declares: one per contig used by makeCatalog().
class StubReference : public Reference
{
public:
    std::string getSequence(const std::string&, int64_t, int64_t) override { return ""; }
    std::string getSequence(const GenomicRegion&) override { return ""; }
    void loadContigIntoCache(const std::string&) override {}
    void clearContigCache() override {}
    const ReferenceContigInfo& contigInfo() const override { return contigInfo_; }

private:
    ReferenceContigInfo contigInfo_{ { { "chr1", 10000 }, { "chr2", 10000 }, { "chr3", 10000 } } };
};

class ResumeTest : public ::testing::Test
{
protected:
    void SetUp() override { cleanUp(); }
    void TearDown() override { cleanUp(); }

    void cleanUp()
    {
        removeIfPresent(processedLociPath(kTestPrefix));
        removeIfPresent(processedLociPath(kTestPrefix) + ".trimmed");
        for (std::int32_t contigIndex = 0; contigIndex < 4; ++contigIndex)
        {
            removeIfPresent(contigTempJsonPath(kTestPrefix, contigIndex));
            removeIfPresent(contigTempVcfPath(kTestPrefix, contigIndex));
        }
    }

    // Appends a locus's output to its contig's temp files the way the genotyping writers do (a new file
    // starts with the document header, and every JSON record but the first is preceded by ", "), then lists
    // the locus with the files' new sizes, as a --resume run does after each locus.
    void finishLocus(
        ResumeCheckpointWriter& writer, const std::string& locusId, std::int32_t contigIndex, bool writesRecord = true)
    {
        const std::string jsonPath = contigTempJsonPath(kTestPrefix, contigIndex);
        const std::string vcfPath = contigTempVcfPath(kTestPrefix, contigIndex);
        if (!boost::filesystem::exists(jsonPath))
        {
            appendToFile(jsonPath, jsonDocumentHeader(sampleParams_));
            appendToFile(vcfPath, vcfHeader());
        }
        if (writesRecord)
        {
            const bool firstRecord = readFile(jsonPath) == jsonDocumentHeader(sampleParams_);
            appendToFile(jsonPath, (firstRecord ? "" : ", ") + jsonRecordFor(locusId));
            appendToFile(vcfPath, vcfLineFor(locusId));
        }
        writer.recordLocus(
            locusId, boost::filesystem::file_size(jsonPath), boost::filesystem::file_size(vcfPath));
    }

    // Lists L1 and L2 (contig 0) and L3 (contig 1) as finished, as an interrupted run would have left them.
    void finishThreeLoci()
    {
        ResumeCheckpointWriter writer(kTestPrefix, kSignature, false, 0);
        finishLocus(writer, "L1", 0);
        finishLocus(writer, "L2", 0);
        finishLocus(writer, "L3", 1);
    }

    ResumeLoadResult loadCheckpoint(const LocusDescriptionCatalog& catalog = makeCatalog())
    {
        return loadResumeCheckpoint(kTestPrefix, sampleParams_, reference_.contigInfo(), catalog);
    }

    std::string jsonHeader() const { return jsonDocumentHeader(sampleParams_); }
    std::string vcfHeader() const
    {
        return vcfDocumentHeader(sampleParams_.id(), reference_.contigInfo(), contigsWithLoci(makeCatalog()));
    }

    SampleParameters sampleParams_{ "sample", Sex::kFemale };
    StubReference reference_;
};

TEST_F(ResumeTest, PathsAreDerivedFromTheOutputPrefix)
{
    EXPECT_EQ(processedLociPath("out"), "out.processed_loci.txt");
    EXPECT_EQ(contigTempJsonPath("out", 3), "out.contig3.json");
    EXPECT_EQ(contigTempVcfPath("out", 3), "out.contig3.vcf");
}

TEST_F(ResumeTest, TheListIsWhatMarksACheckpoint)
{
    EXPECT_FALSE(resumeCheckpointExists(kTestPrefix));
    ResumeCheckpointWriter writer(kTestPrefix, kSignature, false, 0);
    EXPECT_TRUE(resumeCheckpointExists(kTestPrefix));
}

TEST_F(ResumeTest, AllFinishedLociAreRecovered)
{
    finishThreeLoci();
    const std::string contig0Json = readFile(contigTempJsonPath(kTestPrefix, 0));
    const std::string contig0Vcf = readFile(contigTempVcfPath(kTestPrefix, 0));

    const ResumeLoadResult loaded = loadCheckpoint();

    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L1", "L2", "L3" }));
    // Contigs 0 (L1, L2) and 1 (L3) have finished loci; contig 2 (L4) has none.
    ASSERT_EQ(loaded.resumedContigs.size(), 2u);
    EXPECT_TRUE(loaded.resumedContigs.at(0).hasJsonRecords);
    EXPECT_TRUE(loaded.resumedContigs.at(1).hasJsonRecords);
    // Nothing followed the last finished locus, so nothing is cut off.
    EXPECT_EQ(readFile(contigTempJsonPath(kTestPrefix, 0)), contig0Json);
    EXPECT_EQ(readFile(contigTempVcfPath(kTestPrefix, 0)), contig0Vcf);
}

TEST_F(ResumeTest, TempFilesAreCutBackToTheLastFinishedLocus)
{
    finishThreeLoci();
    const std::string expectedJson = jsonHeader() + jsonRecordFor("L1") + ", " + jsonRecordFor("L2");
    const std::string expectedVcf = vcfHeader() + vcfLineFor("L1") + vcfLineFor("L2");

    // What a kill part way through the next locus leaves behind: half a record in each file. A writer that
    // unwound would have added the closing footer instead; either way it has to go.
    appendToFile(contigTempJsonPath(kTestPrefix, 0), ", \n    \"L9\": {\n      \"LocusId\": ");
    appendToFile(contigTempVcfPath(kTestPrefix, 0), "chr\t500\t.\tA");

    loadCheckpoint();

    EXPECT_EQ(readFile(contigTempJsonPath(kTestPrefix, 0)), expectedJson);
    EXPECT_EQ(readFile(contigTempVcfPath(kTestPrefix, 0)), expectedVcf);
}

TEST_F(ResumeTest, ALocusThatWroteNothingIsStillFinished)
{
    {
        ResumeCheckpointWriter writer(kTestPrefix, kSignature, false, 0);
        finishLocus(writer, "L1", 0, false);  // e.g. dropped by --skip-hom-ref
        finishLocus(writer, "L3", 1);
    }

    const ResumeLoadResult loaded = loadCheckpoint();

    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L1", "L3" }));
    // Contig 0 is resumed, but its JSON temp holds only the header, so the next record needs no separator.
    ASSERT_EQ(loaded.resumedContigs.count(0), 1u);
    EXPECT_FALSE(loaded.resumedContigs.at(0).hasJsonRecords);
    EXPECT_EQ(readFile(contigTempJsonPath(kTestPrefix, 0)), jsonHeader());
}

TEST_F(ResumeTest, LociWhoseOutputIsNotOnDiskAreRedoneFromThereOnTheirContigOnly)
{
    finishThreeLoci();
    // A machine crash can lose file bytes the list already counts. Here contig 0 loses L2's record.
    const std::string contig0JsonAfterL1 = jsonHeader() + jsonRecordFor("L1");
    boost::filesystem::resize_file(contigTempJsonPath(kTestPrefix, 0), contig0JsonAfterL1.size());

    const ResumeLoadResult loaded = loadCheckpoint();

    // L2 is redone; L3 is on another contig and unaffected.
    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L1", "L3" }));
    EXPECT_EQ(readFile(contigTempJsonPath(kTestPrefix, 0)), contig0JsonAfterL1);
    EXPECT_EQ(readFile(contigTempVcfPath(kTestPrefix, 0)), vcfHeader() + vcfLineFor("L1"));

    // The list no longer names L2, so a later resume cannot count it as finished either.
    const std::vector<std::string> lines = readLines(processedLociPath(kTestPrefix));
    ASSERT_EQ(lines.size(), 3u);
    EXPECT_EQ(lines[1].substr(0, 3), "L1\t");
    EXPECT_EQ(lines[2].substr(0, 3), "L3\t");
}

TEST_F(ResumeTest, AMissingTempFileMeansItsContigIsRedone)
{
    finishThreeLoci();
    removeIfPresent(contigTempJsonPath(kTestPrefix, 0));

    const ResumeLoadResult loaded = loadCheckpoint();

    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L3" }));
    EXPECT_EQ(loaded.resumedContigs.count(0), 0u);
}

TEST_F(ResumeTest, APartlyWrittenLastListLineIsIgnored)
{
    finishThreeLoci();
    appendToFile(processedLociPath(kTestPrefix), "L4\t12");  // killed before the line was complete

    const ResumeLoadResult loaded = loadCheckpoint();

    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L1", "L2", "L3" }));
    // Rewritten without it, so the resumed run's first line does not get glued onto the fragment.
    const std::string list = readFile(processedLociPath(kTestPrefix));
    EXPECT_EQ(list.find("L4"), std::string::npos);
    EXPECT_EQ(list.back(), '\n');
}

TEST_F(ResumeTest, AMalformedListLineIsRejectedAsCorrupt)
{
    finishThreeLoci();
    appendToFile(processedLociPath(kTestPrefix), "L4\tnot-a-size\t12\n");

    EXPECT_THROW(loadCheckpoint(), ResumeCheckpointCorruptError);
}

TEST_F(ResumeTest, AListWithoutASignatureLineIsRejectedAsCorrupt)
{
    appendToFile(processedLociPath(kTestPrefix), "L1\t100\t200\n");

    EXPECT_THROW(readCheckpointRunSignature(kTestPrefix), ResumeCheckpointCorruptError);
    EXPECT_THROW(loadCheckpoint(), ResumeCheckpointCorruptError);
}

TEST_F(ResumeTest, ALocusMissingFromTheCatalogIsRejected)
{
    finishThreeLoci();
    LocusDescriptionCatalog catalog = makeCatalog();
    catalog.erase(catalog.begin() + 1);  // L2

    try
    {
        loadCheckpoint(catalog);
        FAIL() << "expected the checkpoint to be rejected";
    }
    catch (const ResumeCheckpointCorruptError&)
    {
        FAIL() << "a checkpoint from a different run must be reported, not discarded as corrupt";
    }
    catch (const std::runtime_error&)
    {
    }
}

TEST_F(ResumeTest, AResumedListCanBeAppendedToAndResumedAgain)
{
    finishThreeLoci();
    loadCheckpoint();
    {
        ResumeCheckpointWriter writer(kTestPrefix, kSignature, true, 0);
        finishLocus(writer, "L4", 2);
    }

    const ResumeLoadResult loaded = loadCheckpoint();

    EXPECT_EQ(loaded.doneLocusIds, (std::unordered_set<std::string>{ "L1", "L2", "L3", "L4" }));
    EXPECT_EQ(readLines(processedLociPath(kTestPrefix)).size(), 5u);
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

TEST_F(ResumeTest, TheRunSignatureSurvivesRepeatedResumes)
{
    finishThreeLoci();
    for (int resumeCount = 0; resumeCount != 3; ++resumeCount)
    {
        loadCheckpoint();
        EXPECT_EQ(readCheckpointRunSignature(kTestPrefix), kSignature);
    }
}

TEST_F(ResumeTest, AnOversizedRunSignatureIsStillReadable)
{
    // The signature quotes options verbatim, so a long --locus list makes it very long. Failing to read it
    // back would discard the whole checkpoint on every resume.
    nlohmann::json signature;
    signature["Catalog"] = "catalog.json";
    signature["Locus"] = std::string(200 * 1024, 'x');
    const std::string signatureText = signature.dump();
    {
        ResumeCheckpointWriter writer(kTestPrefix, signatureText, false, 0);
        finishLocus(writer, "L1", 0);
    }

    EXPECT_EQ(describeSignatureMismatch(signatureText, readCheckpointRunSignature(kTestPrefix)), "");
}

TEST_F(ResumeTest, AMultiLineRunSignatureIsRefused)
{
    EXPECT_THROW(ResumeCheckpointWriter(kTestPrefix, "{\n}", false, 0), std::logic_error);
}

TEST_F(ResumeTest, ACatalogWithDuplicateLocusIdsIsRejected)
{
    LocusDescriptionCatalog catalog = makeCatalog();
    EXPECT_NO_THROW(assertCatalogIsResumable(catalog));

    // Resume keys on LocusId alone, so finishing one of two entries sharing an id would mark both done.
    catalog.push_back(makeLocus("L2", 3, 9000));
    EXPECT_THROW(assertCatalogIsResumable(catalog), std::runtime_error);
}

TEST_F(ResumeTest, WritersReportTheirFileSizeAfterFlushing)
{
    const std::string jsonPath = contigTempJsonPath(kTestPrefix, 0);
    const std::string vcfPath = contigTempVcfPath(kTestPrefix, 0);
    const std::set<int32_t> headerContigs = contigsWithLoci(makeCatalog());
    std::uintmax_t jsonSizeAfterL1 = 0;
    {
        IterativeJsonWriter jsonWriter(sampleParams_, reference_.contigInfo(), jsonPath);
        IterativeVcfWriter vcfWriter(sampleParams_.id(), reference_, headerContigs, vcfPath);
        EXPECT_EQ(jsonWriter.flushAndGetFileSize(), jsonHeader().size());
        EXPECT_EQ(vcfWriter.flushAndGetFileSize(), vcfHeader().size());

        jsonWriter.addSkippedRecord("L1", "error");
        jsonSizeAfterL1 = jsonWriter.flushAndGetFileSize();
        EXPECT_EQ(jsonSizeAfterL1, readFile(jsonPath).size());
        EXPECT_GT(jsonSizeAfterL1, jsonHeader().size());
    }

    // Reopened in append mode, as a resumed run reopens a temp file it cut back. Before anything new is
    // written, the size still has to count what the file already holds.
    boost::filesystem::resize_file(jsonPath, jsonSizeAfterL1);
    IterativeJsonWriter jsonWriter(sampleParams_, reference_.contigInfo(), jsonPath, false, nullptr, 0, 1,
        AnalysisMode::kOptimizedStreaming, "", JsonOutputMode::kAppendAfterHeader, true);
    IterativeVcfWriter vcfWriter(
        sampleParams_.id(), reference_, headerContigs, vcfPath, VcfOutputMode::kAppendAfterHeader);
    EXPECT_EQ(jsonWriter.flushAndGetFileSize(), jsonSizeAfterL1);
    EXPECT_EQ(vcfWriter.flushAndGetFileSize(), vcfHeader().size());

    jsonWriter.addSkippedRecord("L2", "error");
    EXPECT_EQ(jsonWriter.flushAndGetFileSize(), readFile(jsonPath).size());
    EXPECT_NE(readFile(jsonPath).find("}, \n    \"L2\""), std::string::npos);
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
