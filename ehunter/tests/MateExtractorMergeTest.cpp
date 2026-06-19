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

#include <stdexcept>
#include <string>

#include "gtest/gtest.h"

#include "core/Read.hh"
#include "sample/MateExtractor.hh"

using namespace ehunter;
using ehunter::htshelpers::MateExtractor;

namespace
{
// A committed multi-contig fixture; MateExtractor's ctor opens the BAM (its file handle is unused by the
// addMateToCache / mergeAndFreeze paths exercised here, but the ctor still requires a readable file).
const std::string kDataDir = std::string(EH_TEST_DATA_DIR) + "/parallel_processing_fixtures";

MateExtractor makeShard()
{
    return MateExtractor(kDataDir + "/reads.bam", kDataDir + "/reads.bam.bai", kDataDir + "/reference.fa", true);
}

FullRead makeFullRead(const std::string& fragmentId, MateNumber mateNumber, const std::string& sequence)
{
    return FullRead(Read(ReadId(fragmentId, mateNumber), sequence, false), LinearAlignmentStats{});
}
}  // namespace

// Scheme (D): each chromosome-stride prepass worker fills a PRIVATE shard; mergeAndFreeze splices the
// key-disjoint shards into one immutable cache. Read1 vs read2 of the same fragment are distinct keys.
TEST(MergeAndFreeze, DisjointShardsSpliceIntoUnion)
{
    MateExtractor shardA = makeShard();
    MateExtractor shardB = makeShard();
    shardA.addMateToCache(ReadId("f1", MateNumber::kFirstMate), makeFullRead("f1", MateNumber::kFirstMate, "ACGTACGT"));
    shardA.addMateToCache(ReadId("f2", MateNumber::kSecondMate), makeFullRead("f2", MateNumber::kSecondMate, "TTTTTTTT"));
    shardB.addMateToCache(ReadId("f3", MateNumber::kFirstMate), makeFullRead("f3", MateNumber::kFirstMate, "GGGGGGGG"));
    // Same fragment id as f1 but the OTHER mate number — a distinct key, must not collide on merge.
    shardB.addMateToCache(ReadId("f1", MateNumber::kSecondMate), makeFullRead("f1", MateNumber::kSecondMate, "CCCCCCCC"));

    std::shared_ptr<const htshelpers::MateCache> frozen = MateExtractor::mergeAndFreeze({ &shardA, &shardB });

    EXPECT_EQ(frozen->size(), 4u);
    EXPECT_EQ(frozen->count(ReadId("f1", MateNumber::kFirstMate)), 1u);
    EXPECT_EQ(frozen->count(ReadId("f1", MateNumber::kSecondMate)), 1u);
    EXPECT_EQ(frozen->count(ReadId("f2", MateNumber::kSecondMate)), 1u);
    EXPECT_EQ(frozen->count(ReadId("f3", MateNumber::kFirstMate)), 1u);
    // Value preserved through the node splice.
    EXPECT_EQ(frozen->at(ReadId("f1", MateNumber::kFirstMate)).r.sequence(), "ACGTACGT");
    EXPECT_EQ(frozen->at(ReadId("f3", MateNumber::kFirstMate)).r.sequence(), "GGGGGGGG");

    // Shards are emptied after splicing (storage released).
    EXPECT_EQ(shardA.mateCacheSize(), 0u);
    EXPECT_EQ(shardB.mateCacheSize(), 0u);
}

// Defense-in-depth: invariant 3 says a key can never appear in two shards (one primary alignment per read),
// but a buggy aligner could emit duplicate primaries — mergeAndFreeze must throw, not silently drop one.
TEST(MergeAndFreeze, DuplicateKeyAcrossShardsThrows)
{
    MateExtractor shardA = makeShard();
    MateExtractor shardB = makeShard();
    shardA.addMateToCache(ReadId("dup", MateNumber::kFirstMate), makeFullRead("dup", MateNumber::kFirstMate, "ACGT"));
    shardB.addMateToCache(ReadId("dup", MateNumber::kFirstMate), makeFullRead("dup", MateNumber::kFirstMate, "TGCA"));

    EXPECT_THROW(MateExtractor::mergeAndFreeze({ &shardA, &shardB }), std::logic_error);
}
