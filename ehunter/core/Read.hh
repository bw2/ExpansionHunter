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
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include <boost/functional/hash.hpp>

#include "alignment/AlignmentClassifier.hh"
#include "core/GenomicRegion.hh"
#include "graphalign/GraphAlignment.hh"
#include "graphutils/SequenceOperations.hh"

namespace ehunter
{

using FragmentId = std::string;

enum class MateNumber
{
    kFirstMate = 1,
    kSecondMate = 2
};

class ReadId
{
public:
    ReadId(FragmentId fragmentId, MateNumber mateNumber)
        : fragmentId_(std::move(fragmentId))
        , mateNumber_(mateNumber)
    {
        if (fragmentId_.empty())
        {
            throw std::logic_error("Encountered an empty fragment id");
        }
    }

    const FragmentId& fragmentId() const { return fragmentId_; }
    MateNumber mateNumber() const { return mateNumber_; }

    bool operator==(const ReadId& other) const
    {
        return fragmentId_ == other.fragmentId_ && mateNumber_ == other.mateNumber_;
    }

    friend std::size_t hash_value(const ReadId& readId)
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, readId.fragmentId_);
        boost::hash_combine(seed, static_cast<int>(readId.mateNumber_));

        return seed;
    }

    std::string toString() const {
        return fragmentId_ + "/" + std::to_string(static_cast<int>(mateNumber_));
    }

private:
    FragmentId fragmentId_;
    MateNumber mateNumber_;
};

std::ostream& operator<<(std::ostream& out, const ReadId& readId);

class Read
{
public:
    Read(ReadId readId, std::string sequence, bool isReversed)
        : readId_(std::move(readId))
        , sequence_(std::move(sequence))
        , isReversed_(isReversed)
    {
        if (sequence_.empty())
        {
            std::ostringstream encoding;
            encoding << readId_;
            throw std::logic_error("Encountered empty query for " + encoding.str());
        }
    }

    const ReadId& readId() const { return readId_; }
    const FragmentId& fragmentId() const { return readId_.fragmentId(); }
    MateNumber mateNumber() const { return readId_.mateNumber(); }
    const std::string& sequence() const { return sequence_; }

    bool isFirstMate() const { return mateNumber() == MateNumber::kFirstMate; }
    bool isSecondMate() const { return mateNumber() == MateNumber::kSecondMate; }
    // Return whether the read is reverse complemented relative to its
    //  original direction during sequencing
    bool isReversed() const { return isReversed_; }

    void reverseComplement()
    {
        sequence_ = graphtools::reverseComplement(sequence_);
        isReversed_ = !isReversed_;
    }

private:
    ReadId readId_;
    std::string sequence_;
    bool isReversed_;
};

struct LinearAlignmentStats
{
    int32_t chromId = -1;
    int32_t pos = -1;
    int32_t mapq = -1;
    int32_t mateChromId = -1;
    int32_t matePos = -1;
    bool isPaired = false;
    bool isMapped = false;
    bool isMateMapped = false;
    bool isSupplementaryAlignment = false;
    bool isSecondaryAlignment = false;
    std::vector<uint32_t> cigar;
};

std::ostream& operator<<(std::ostream& out, const LinearAlignmentStats& alignmentStats);
bool operator==(const Read& read, const Read& mate);
bool operator==(const LinearAlignmentStats& statsA, const LinearAlignmentStats& statsB);

class RepeatAlignmentStats
{
public:
    RepeatAlignmentStats(
        const GraphAlignment& canonical_alignment, AlignmentType canonical_alignment_type,
        int32_t num_repeat_units_spanned)
        : canonical_alignment_(canonical_alignment)
        , canonical_alignment_type_(canonical_alignment_type)
        , num_repeat_units_spanned_(num_repeat_units_spanned)
    {
    }

    const GraphAlignment& canonicalAlignment() const { return canonical_alignment_; }
    AlignmentType canonicalAlignmentType() const { return canonical_alignment_type_; }
    int32_t numRepeatUnitsSpanned() const { return num_repeat_units_spanned_; }

private:
    GraphAlignment canonical_alignment_;
    AlignmentType canonical_alignment_type_;
    int32_t num_repeat_units_spanned_;
};

using ReadIdToRepeatAlignmentStats = std::unordered_map<std::string, RepeatAlignmentStats>;

std::ostream& operator<<(std::ostream& out, const Read& read);

struct FullRead {
    FullRead(const Read& read, const LinearAlignmentStats& alignmentStats) :
        r(std::move(read)), s(std::move(alignmentStats))
    {}

    Read r;
    LinearAlignmentStats s;
};

struct FullReadPair
{
    FullReadPair() = default;
    FullReadPair(FullRead read, FullRead mate) : firstMate(std::move(read)), secondMate(std::move(mate)) {}

    int numMatesSet() const
    {
        return static_cast<int>(firstMate != std::nullopt) + static_cast<int>(secondMate != std::nullopt);
    }

    std::optional<FullRead> firstMate = std::nullopt;
    std::optional<FullRead> secondMate = std::nullopt;
};

bool operator==(const FullRead& read_a, const FullRead& read_b);
bool operator==(const FullReadPair& readPair_a, const FullReadPair& readPair_b);

}
