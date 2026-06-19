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

#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "core/GenomicRegion.hh"
#include "core/ReferenceContigInfo.hh"
#include "io/StreamingOutputMerge.hh"

namespace ehunter
{
// Free helpers defined (with external linkage, not header-declared) in HtsLowMemStreamingSampleAnalysis.cpp.
std::vector<GenomicRegion> wholeContigRegionsForStride(const ReferenceContigInfo& contigInfo, int firstContig, int strideT);
GenomicRegion clampedContigRegion(int32_t contigIndex, int64_t minStart, int64_t maxEnd, int64_t contigSize, int halo);
}

using namespace ehunter;

namespace
{
ReferenceContigInfo makeContigInfo()
{
    // 5 contigs of distinct sizes, ascending index order.
    return ReferenceContigInfo({ { "chr1", 1000 }, { "chr2", 2000 }, { "chr3", 3000 }, { "chr4", 4000 }, { "chr5", 5000 } });
}

std::string readFile(const std::string& path)
{
    std::ifstream in(path, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

void writeFile(const std::string& path, const std::string& content)
{
    std::ofstream out(path, std::ios::binary);
    out << content;
}

// IterativeJsonWriter wrapper format that mergeRegionJsonFiles slices on.
std::string jsonRegionFile(const std::string& body)
{
    return "{\n  \"SampleParameters\": {\n    \"SampleId\": \"s\"\n  },\n  \"LocusResults\": {" + body + "\n  }\n}\n";
}
}  // namespace

TEST(WholeContigRegionsForStride, StrideAssignmentAscendingWholeContig)
{
    const ReferenceContigInfo contigInfo = makeContigInfo();

    // T=2, worker 0 owns contigs 0,2,4; worker 1 owns 1,3.
    const std::vector<GenomicRegion> worker0 = wholeContigRegionsForStride(contigInfo, 0, 2);
    ASSERT_EQ(worker0.size(), 3u);
    EXPECT_EQ(worker0[0].contigIndex(), 0);
    EXPECT_EQ(worker0[0].start(), 0);
    EXPECT_EQ(worker0[0].end(), 1000);
    EXPECT_EQ(worker0[1].contigIndex(), 2);
    EXPECT_EQ(worker0[1].end(), 3000);
    EXPECT_EQ(worker0[2].contigIndex(), 4);
    EXPECT_EQ(worker0[2].end(), 5000);

    const std::vector<GenomicRegion> worker1 = wholeContigRegionsForStride(contigInfo, 1, 2);
    ASSERT_EQ(worker1.size(), 2u);
    EXPECT_EQ(worker1[0].contigIndex(), 1);
    EXPECT_EQ(worker1[1].contigIndex(), 3);

    // T > numContigs: a worker whose firstContig is past the end owns nothing.
    EXPECT_TRUE(wholeContigRegionsForStride(contigInfo, 7, 8).empty());
    // T=1: one worker owns every contig.
    EXPECT_EQ(wholeContigRegionsForStride(contigInfo, 0, 1).size(), 5u);
}

TEST(ClampedContigRegion, HaloClampedToContigBounds)
{
    // Interior locus: full [minStart-halo, maxEnd+halo].
    const GenomicRegion interior = clampedContigRegion(2, 1500, 1800, 3000, 1000);
    EXPECT_EQ(interior.contigIndex(), 2);
    EXPECT_EQ(interior.start(), 500);
    EXPECT_EQ(interior.end(), 2800);

    // Near contig start: clamped to 0.
    const GenomicRegion atStart = clampedContigRegion(0, 200, 400, 1000, 1000);
    EXPECT_EQ(atStart.start(), 0);
    EXPECT_EQ(atStart.end(), 1000);

    // Near contig end: clamped to contigSize.
    const GenomicRegion atEnd = clampedContigRegion(0, 600, 900, 1000, 1000);
    EXPECT_EQ(atEnd.start(), 0);
    EXPECT_EQ(atEnd.end(), 1000);
}

TEST(MergeRegionJsonFiles, AscendingBodiesJoinedSkippingEmpty)
{
    const std::string f0 = "/tmp/eh_chrstride_json_0.json";
    const std::string f1 = "/tmp/eh_chrstride_json_1.json";
    const std::string f2 = "/tmp/eh_chrstride_json_2.json";
    const std::string merged = "/tmp/eh_chrstride_json_merged.json";

    writeFile(f0, jsonRegionFile("\n    \"L0\": {\"a\": 0}"));
    writeFile(f1, jsonRegionFile(""));  // header-only / zero-record interior file
    writeFile(f2, jsonRegionFile("\n    \"L2\": {\"a\": 2}"));

    mergeRegionJsonFiles(merged, { f0, f1, f2 });
    const std::string out = readFile(merged);

    // SampleParameters prefix from file[0], bodies joined with ", ", empty body skipped (no double comma).
    const std::string expected =
        "{\n  \"SampleParameters\": {\n    \"SampleId\": \"s\"\n  },\n  \"LocusResults\": {"
        "\n    \"L0\": {\"a\": 0}, \n    \"L2\": {\"a\": 2}"
        "\n  }\n}\n";
    EXPECT_EQ(out, expected);

    // Path-vector order is load-bearing: a non-ascending order yields different bytes.
    const std::string mergedSwapped = "/tmp/eh_chrstride_json_merged_swapped.json";
    mergeRegionJsonFiles(mergedSwapped, { f2, f1, f0 });
    EXPECT_NE(readFile(mergedSwapped), out);

    for (const std::string& p : { f0, f1, f2, merged, mergedSwapped }) { std::remove(p.c_str()); }
}

TEST(MergeRegionJsonFiles, LeadingEmptyFileStillSuppliesSampleParameters)
{
    const std::string f0 = "/tmp/eh_chrstride_json_lead0.json";
    const std::string f1 = "/tmp/eh_chrstride_json_lead1.json";
    const std::string merged = "/tmp/eh_chrstride_json_lead_merged.json";

    writeFile(f0, jsonRegionFile(""));  // empty body first
    writeFile(f1, jsonRegionFile("\n    \"L1\": {\"a\": 1}"));

    mergeRegionJsonFiles(merged, { f0, f1 });
    const std::string out = readFile(merged);

    // No stray leading ", " before the first non-empty body; SampleParameters still emitted.
    const std::string expected =
        "{\n  \"SampleParameters\": {\n    \"SampleId\": \"s\"\n  },\n  \"LocusResults\": {"
        "\n    \"L1\": {\"a\": 1}"
        "\n  }\n}\n";
    EXPECT_EQ(out, expected);

    for (const std::string& p : { f0, f1, merged }) { std::remove(p.c_str()); }
}

TEST(MergeRegionJsonFiles, SingleFileRoundTrips)
{
    const std::string f0 = "/tmp/eh_chrstride_json_single.json";
    const std::string merged = "/tmp/eh_chrstride_json_single_merged.json";
    const std::string content = jsonRegionFile("\n    \"L0\": {\"a\": 0}");
    writeFile(f0, content);
    mergeRegionJsonFiles(merged, { f0 });
    EXPECT_EQ(readFile(merged), content);
    for (const std::string& p : { f0, merged }) { std::remove(p.c_str()); }
}

TEST(MergeRegionVcfFiles, FirstFileVerbatimRestRecordsOnly)
{
    const std::string v0 = "/tmp/eh_chrstride_vcf_0.vcf";
    const std::string v1 = "/tmp/eh_chrstride_vcf_1.vcf";
    const std::string v2 = "/tmp/eh_chrstride_vcf_2.vcf";
    const std::string merged = "/tmp/eh_chrstride_vcf_merged.vcf";

    writeFile(v0, "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\nchr1\t10\tL0\n");
    writeFile(v1, "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\n");  // header-only (zero records)
    writeFile(v2, "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\nchr2\t20\tL2\n");

    mergeRegionVcfFiles(merged, { v0, v1, v2 });
    const std::string out = readFile(merged);

    const std::string expected =
        "##fileformat=VCFv4.1\n#CHROM\tPOS\tID\nchr1\t10\tL0\n"  // file[0] verbatim (header + record)
        "chr2\t20\tL2\n";                                        // file[2] record only; file[1] contributes nothing
    EXPECT_EQ(out, expected);

    for (const std::string& p : { v0, v1, v2, merged }) { std::remove(p.c_str()); }
}
