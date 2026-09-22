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

#include "io/VcfHeader.hh"

#include <set>
#include <sstream>
#include <string>

#include "gtest/gtest.h"

#include "core/ReferenceContigInfo.hh"
#include "io/IterativeVcfWriter.hh"
#include "locus/LocusSpecification.hh"

using namespace ehunter;

namespace
{

LocusDescription makeLocus(const std::string& locusId, std::int32_t contigIndex, std::int64_t start)
{
    return LocusDescription(locusId, ChromType::kAutosome, "(CAG)*", contigIndex, start - 1000, start + 1000,
        start, start + 30, false, { locusId });
}

const ReferenceContigInfo kContigInfo({ { "chr1", 248956422 }, { "chr2", 242193529 }, { "chrX", 156040895 } });

}

TEST(VcfContigLines, AreWrittenInBamHeaderOrderWithTheirLengths)
{
    std::ostringstream out;
    outputVcfContigLines(kContigInfo, { 2, 0 }, out);
    EXPECT_EQ(out.str(), "##contig=<ID=chr1,length=248956422>\n##contig=<ID=chrX,length=156040895>\n");
}

TEST(IterativeVcfHeader, DeclaresOnlyTheContigsThatHaveLociInTheCatalog)
{
    const LocusDescriptionCatalog catalog{ makeLocus("A", 2, 5000), makeLocus("B", 0, 1000),
        makeLocus("C", 2, 9000) };
    EXPECT_EQ(contigsWithLoci(catalog), (std::set<int32_t>{ 0, 2 }));

    const std::string header = vcfDocumentHeader("sample", kContigInfo, contigsWithLoci(catalog));
    EXPECT_NE(
        header.find("##contig=<ID=chr1,length=248956422>\n##contig=<ID=chrX,length=156040895>\n#CHROM\t"),
        std::string::npos);
    EXPECT_EQ(header.find("chr2"), std::string::npos);
}
