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

#include "locus/MotifComposition.hh"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <unordered_map>
#include <unordered_set>

#include <htslib/sam.h>

#include "core/Read.hh"
#include "graphutils/BaseMatching.hh"
#include "graphutils/SequenceOperations.hh"
#include "spdlog/spdlog.h"

using std::string;
using std::vector;

namespace ehunter
{

bool isEligibleForMotifComposition(int motifLength, int typicalReadLength)
{
    return motifLength >= 2 && 3 * motifLength <= typicalReadLength;
}

namespace motifcomposition
{

namespace
{

// When no window matches exactly, the shift search and the start-offset search fall back to counting mismatches
// against the catalog motif and at most this many of the most common known motifs. Exact matches are always checked
// against every known motif.
const size_t kMaxMotifsInShiftSearch = 8;

// Per-base error rate of high-quality bases, used to predict how often a sequencing error turns a
// common motif into a one-base variant of it.
const double kMotifCompositionBaseErrorRate = 0.001;
// Target chance that sequencing errors alone add a false motif at a locus.
const double kMotifCompositionLocusFalsePositiveRate = 1e-4;
// A new motif must be seen in at least this many read pairs, and in at least this fraction of the read pairs that
// have any motif-sized substring at the locus.
const int kMotifCompositionMinReadPairs = 2;
const double kMotifCompositionMinReadPairFraction = 0.005;
// A read counts as spanning only with at least this many aligned bases on each side of the repeat.
const int kMotifCompositionMinFlankBases = 10;
// Minimum fraction of bases that equal the base one motif length earlier, for a read to count as an in-repeat read,
// and for a soft clip to be kept as repeat sequence.
const double kMotifCompositionMinInrepeatReadPeriodScore = 0.75;
const double kMotifCompositionMinSoftClipPeriodScore = 0.75;

// Bases with quality <= 20 are stored in lower case (see htshelpers::decodeRead). A no-call is never high quality,
// whatever its case: reverse-complementing a read turns a low-quality 'n' into 'N'.
inline bool isHighQuality(char base) { return base < 'a' && base != 'N'; }

inline char toUpperBase(char base) { return (base >= 'a' && base <= 'z') ? static_cast<char>(base - 'a' + 'A') : base; }

inline bool basesMatch(char left, char right)
{
    left = toUpperBase(left);
    return left == toUpperBase(right) && left != 'N';
}

string toUpperSequence(const string& sequence)
{
    string upper(sequence);
    std::transform(upper.begin(), upper.end(), upper.begin(), toUpperBase);
    return upper;
}

// Mismatches between the window starting at `window` and `motif`, which may contain IUPAC codes. Counting stops
// once it exceeds `limit`.
int countMismatches(const char* window, const string& motif, int limit)
{
    int mismatches = 0;
    for (size_t index = 0; index != motif.size(); ++index)
    {
        if (!graphtools::checkIfReferenceBaseMatchesQueryBase(motif[index], window[index]) && ++mismatches > limit)
        {
            break;
        }
    }
    return mismatches;
}

int countFewestMismatches(const char* window, const string& catalogMotif, const vector<string>& knownMotifs, int limit)
{
    int fewest = countMismatches(window, catalogMotif, limit);
    const size_t motifCount = std::min(knownMotifs.size(), kMaxMotifsInShiftSearch);
    for (size_t index = 0; index != motifCount && fewest > 0; ++index)
    {
        fewest = std::min(fewest, countMismatches(window, knownMotifs[index], fewest - 1));
    }
    return fewest;
}

// True if the upper-case window matches the catalog motif (IUPAC codes allowed) or any known motif exactly.
bool matchesMotifExactly(const char* window, const string& catalogMotif, const vector<string>& knownMotifs)
{
    const size_t motifLength = catalogMotif.size();
    return countMismatches(window, catalogMotif, 0) == 0
        || std::any_of(
               knownMotifs.begin(), knownMotifs.end(),
               [&](const string& motif) { return motif.compare(0, motifLength, window, motifLength) == 0; });
}

// Offset in [0, k) at which to tile motifs over the upper-case sequence: the offset with the most windows that match
// a motif exactly, then the fewest mismatches over the other windows (against the catalog motif and the most common
// known motifs), then the smallest offset. Every offset is scored over the same number of whole windows.
int computeFrameOffsetAgainst(const string& sequence, const string& catalogMotif, const vector<string>& knownMotifs)
{
    const int motifLength = catalogMotif.size();
    const int length = sequence.size();
    if (length < motifLength)
    {
        return 0;
    }
    const int windowCount = std::max(1, (length - motifLength + 1) / motifLength);
    const int lastOffset = std::min(motifLength - 1, length - windowCount * motifLength);
    int bestOffset = 0;
    int bestExactWindows = -1;
    int bestMismatches = INT_MAX;
    for (int offset = 0; offset <= lastOffset; ++offset)
    {
        int exactWindows = 0;
        int mismatches = 0;
        for (int window = 0; window != windowCount; ++window)
        {
            const char* windowStart = sequence.data() + offset + window * motifLength;
            if (matchesMotifExactly(windowStart, catalogMotif, knownMotifs))
            {
                ++exactWindows;
            }
            else
            {
                mismatches += countFewestMismatches(windowStart, catalogMotif, knownMotifs, motifLength);
            }
        }
        if (exactWindows > bestExactWindows || (exactWindows == bestExactWindows && mismatches < bestMismatches))
        {
            bestExactWindows = exactWindows;
            bestMismatches = mismatches;
            bestOffset = offset;
        }
    }
    return bestOffset;
}

// Checks whether a partial motif at the beginning or end of a sequence matches the motif(s) before it. 
// - A partial motif at the end of a sequence (isUnitStart true) is compared with the first `length` bases of a
//   whole motif; a partial motif at the beginning of a sequence (isUnitStart false) is compared with the last `length` bases of a whole motif.
// - It passes if it matches the catalog motif (IUPAC codes allowed) or any of the known motifs.
// - It must have at least min(4, motifLength - 1) bases and fewer than motifLength. Shorter partial units say
//   too little to count. For example, at a 3 bp motif only a 2-base partial unit is judged.
// - It may have one mismatch if it has 4 or more bases; a shorter one must match exactly.
// Example: at a CAG locus, "CA" after the last CAG passes, and "TG" does not.
bool partialUnitMatchesMotif(
    const char* partial, int length, bool isUnitStart, const string& catalogMotif, const vector<string>& knownMotifs)
{
    const int motifLength = catalogMotif.size();
    if (length <= 0 || length >= motifLength || length < std::min(4, motifLength - 1))
    {
        return false;
    }
    const int allowedMismatches = length >= 4 ? 1 : 0;
    const int motifOffset = isUnitStart ? 0 : motifLength - length;
    auto matches = [&](const string& motif)
    {
        int mismatches = 0;
        for (int index = 0; index != length && mismatches <= allowedMismatches; ++index)
        {
            mismatches += !graphtools::checkIfReferenceBaseMatchesQueryBase(motif[motifOffset + index], partial[index]);
        }
        return mismatches <= allowedMismatches;
    };
    if (matches(catalogMotif))
    {
        return true;
    }
    return std::any_of(knownMotifs.begin(), knownMotifs.end(), matches);
}

// Motif list with fast membership test; the order (most common first) drives the shift search.
class MotifList
{
public:
    explicit MotifList(vector<string> motifs)
        : motifs_(std::move(motifs))
    {
        for (const string& motif : motifs_)
        {
            members_.insert(motif);
        }
    }

    const vector<string>& motifs() const { return motifs_; }
    bool contains(const string& motif) const { return members_.count(motif) != 0; }

private:
    vector<string> motifs_;
    std::unordered_set<string> members_;
};

// True if sequence[start, start + motifLength) is a single base repeated.
bool isHomopolymer(const string& upperSequence, int start, int motifLength)
{
    for (int index = 1; index != motifLength; ++index)
    {
        if (upperSequence[start + index] != upperSequence[start])
        {
            return false;
        }
    }
    return true;
}

void mergeAdjacentGaps(vector<SequenceSubstring>& sequenceSubstrings)
{
    vector<SequenceSubstring> merged;
    merged.reserve(sequenceSubstrings.size());
    for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
    {
        if (sequenceSubstring.type == SequenceSubstringType::kGap && !merged.empty()
            && merged.back().type == SequenceSubstringType::kGap)
        {
            merged.back().length += sequenceSubstring.length;
        }
        else
        {
            merged.push_back(sequenceSubstring);
        }
    }
    sequenceSubstrings.swap(merged);
}

// See splitIntoMotifs for the arguments. `upperTract` is `tract` in upper case; `tract` keeps the
// case that marks base quality.
void splitTract(
    const string& tract, const string& upperTract, const string& catalogMotif, const MotifList& knownMotifs,
    int startOffset, bool startsAtRepeatEdge, bool endsAtRepeatEdge, const boost::optional<string>& referenceEndPartial,
    vector<SequenceSubstring>& sequenceSubstrings)
{
    sequenceSubstrings.clear();
    const int motifLength = catalogMotif.size();
    const int length = tract.size();
    const vector<string>& orderedMotifs = knownMotifs.motifs();

    const int leadingPartialLength = std::min(std::max(startOffset, 0), length);
    int position = leadingPartialLength;
    if (position > 0)
    {
        sequenceSubstrings.push_back({ 0, position, SequenceSubstringType::kGap });
    }

    string window;
    while (position + motifLength <= length)
    {
        window.assign(upperTract, position, motifLength);
        if (knownMotifs.contains(window))
        {
            sequenceSubstrings.push_back({ position, motifLength, SequenceSubstringType::kKnownMotif });
            position += motifLength;
            continue;
        }
        const char* windowStart = tract.data() + position;
        if (countMismatches(windowStart, catalogMotif, 0) == 0)
        {
            sequenceSubstrings.push_back({ position, motifLength, SequenceSubstringType::kCatalogMotifMatch });
            position += motifLength;
            continue;
        }

        // Look up to one motif length ahead for where the frame resumes, so an indel or a partial unit does not throw
        // off every later substring: the nearest shift whose window matches a known motif or the catalog motif
        // exactly, or, when there is none, the shift whose window has the fewest mismatches. The fallback keeps the
        // frame at repeats whose units often differ from every known motif by a base or two, such as VNTRs.
        const int maxShift = std::min(motifLength - 1, length - position - motifLength);
        int bestShift = 0;
        for (int shift = 1; shift <= maxShift && bestShift == 0; ++shift)
        {
            window.assign(upperTract, position + shift, motifLength);
            if (knownMotifs.contains(window) || countMismatches(windowStart + shift, catalogMotif, 0) == 0)
            {
                bestShift = shift;
            }
        }
        if (bestShift == 0)
        {
            int bestScore = countFewestMismatches(windowStart, catalogMotif, orderedMotifs, motifLength);
            for (int shift = 1; shift <= maxShift && bestScore > 0; ++shift)
            {
                const int score
                    = countFewestMismatches(windowStart + shift, catalogMotif, orderedMotifs, bestScore - 1);
                if (score < bestScore)
                {
                    bestScore = score;
                    bestShift = shift;
                }
            }
        }

        if (bestShift == 0)
        {
            // A substring made of a single base (such as GGG at a CAG locus) is never a new motif: on NovaSeq data,
            // reads with no signal read as runs of G. A substring that repeats a longer unit (such as ATAT at an ATCT
            // locus) can be one.
            const SequenceSubstringType type = isHomopolymer(upperTract, position, motifLength)
                ? SequenceSubstringType::kGap
                : SequenceSubstringType::kNewMotifCandidate;
            sequenceSubstrings.push_back({ position, motifLength, type });
            position += motifLength;
        }
        else
        {
            sequenceSubstrings.push_back({ position, bestShift, SequenceSubstringType::kGap });
            position += bestShift;
        }
    }
    const int trailingPartialStart = position;
    const int trailingPartialLength = length - position;
    if (trailingPartialLength > 0)
    {
        sequenceSubstrings.push_back({ position, trailingPartialLength, SequenceSubstringType::kGap });
    }

    // Is the tract start (with the partial unit before the first substring) a trusted edge?
    const bool startIsTrusted = startsAtRepeatEdge
        || partialUnitMatchesMotif(tract.data(), leadingPartialLength, false, catalogMotif, orderedMotifs);

    // Is the tract end (with the partial unit after the last substring) a trusted edge? At an end the alignment placed
    // at the repeat's edge, the partial unit must be the reference's own partial last unit: at most loci the
    // repeat length is not a multiple of the motif length, so a partial unit there is normal, while one of any
    // other length shows that an indel shifted the frame.
    bool endIsTrusted;
    if (!referenceEndPartial)
    {
        endIsTrusted = true;
    }
    else if (endsAtRepeatEdge)
    {
        endIsTrusted = trailingPartialLength == static_cast<int>(referenceEndPartial->size());
        if (endIsTrusted && trailingPartialLength > 0)
        {
            const char* partial = tract.data() + trailingPartialStart;
            const int mismatchesVsCatalogMotif = [&]
            {
                int mismatches = 0;
                for (int index = 0; index != trailingPartialLength; ++index)
                {
                    mismatches
                        += !graphtools::checkIfReferenceBaseMatchesQueryBase(catalogMotif[index], partial[index]);
                }
                return mismatches;
            }();
            int mismatchesVsReference = 0;
            for (int index = 0; index != trailingPartialLength; ++index)
            {
                mismatchesVsReference += !basesMatch((*referenceEndPartial)[index], partial[index]);
            }
            endIsTrusted = std::min(mismatchesVsCatalogMotif, mismatchesVsReference) <= 1;
        }
    }
    else
    {
        endIsTrusted = partialUnitMatchesMotif(
            tract.data() + trailingPartialStart, trailingPartialLength, true, catalogMotif, orderedMotifs);
    }

    // A possible new motif is kept only if, on each side, it touches another substring or a trusted edge. Repeat until
    // nothing changes, so a run of candidates survives only if both of its ends are anchored.
    auto isSequenceSubstring = [](SequenceSubstringType type) { return type != SequenceSubstringType::kGap; };
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (size_t index = 0; index != sequenceSubstrings.size(); ++index)
        {
            SequenceSubstring& sequenceSubstring = sequenceSubstrings[index];
            if (sequenceSubstring.type != SequenceSubstringType::kNewMotifCandidate)
            {
                continue;
            }
            const bool leftSideIsAnchored = (index > 0 && isSequenceSubstring(sequenceSubstrings[index - 1].type))
                || (sequenceSubstring.offsetWithinRepeatTract == leadingPartialLength && startIsTrusted);
            const bool rightSideIsAnchored
                = (index + 1 < sequenceSubstrings.size() && isSequenceSubstring(sequenceSubstrings[index + 1].type))
                || (sequenceSubstring.offsetWithinRepeatTract + sequenceSubstring.length == trailingPartialStart
                    && endIsTrusted);
            if (!leftSideIsAnchored || !rightSideIsAnchored)
            {
                sequenceSubstring.type = SequenceSubstringType::kGap;
                changed = true;
            }
        }
    }

    mergeAdjacentGaps(sequenceSubstrings);
}

string computeCanonicalRotation(const string& motif)
{
    string best = motif;
    for (size_t shift = 1; shift < motif.size(); ++shift)
    {
        string rotation = motif.substr(shift) + motif.substr(0, shift);
        if (rotation < best)
        {
            best = std::move(rotation);
        }
    }
    return best;
}

int countDifferences(const string& left, const string& right)
{
    int differences = 0;
    for (size_t index = 0; index != left.size(); ++index)
    {
        differences += left[index] != right[index];
    }
    return differences;
}

// ------------------------------------------------------------------------------------------------------------
// Reads and their repeat tracts
// ------------------------------------------------------------------------------------------------------------

enum class ReadKind
{
    kSpanning,
    kFlanking,
    kInsideRepeat, // mapped inside the repeat, reaching neither edge
    kInrepeat // in-repeat read (IRR): made of repeat, placed elsewhere, with a mate anchored next to the repeat
};

struct RepeatTract
{
    string bases; // reference orientation; lower case marks low-quality bases
    int readPairIndex;
    ReadKind kind;
    int lengthForAssignment; // repeat length in bp used to assign the read to an allele
    bool startsAtRepeatEdge;
    bool endsAtRepeatEdge;
    int startOffset; // frame offset of the first substring; -1 when it has to be found from the sequence
};

// What the read's CIGAR string says about where it sits relative to the repeat [S, E).
struct RepeatEdgeAlignment
{
    int64_t referenceStart = 0;
    int64_t referenceEnd = 0;
    int leftClipLength = 0;
    int rightClipLength = 0;
    // Read position of the tract start placed by the alignment at S, and the reference position its first base
    // stands for (S, the first position after a deletion that contains S, or S minus a whole-motif insertion at
    // S). -1 if the alignment does not reach S.
    int readPositionAtStart = -1;
    int64_t referencePositionAtStart = 0;
    // Read position just past the tract end placed by the alignment at E. -1 if the alignment does not reach E.
    int readPositionAtEnd = -1;
    // True when the tract starts with, or ends with, a whole-motif insertion the aligner put exactly on the edge.
    // The aligner's choice of that spot does not show whether the inserted bases belong to the repeat or to the
    // flank next to it, so such an end is not trusted to anchor a new motif.
    bool startsWithEdgeInsertion = false;
    bool endsWithEdgeInsertion = false;
    int flankBasesBefore = 0;
    int flankBasesAfter = 0;
    // Aligned bases inside [S, E) plus whole-motif insertions near the repeat, computed the way the fast path
    // computes its size vote (processRead in sample/HtsLowMemStreamingHelpers.cpp).
    int lengthForAssignment = 0;
};

int64_t countOverlap(int64_t start, int64_t end, int64_t otherStart, int64_t otherEnd)
{
    return std::max<int64_t>(0, std::min(end, otherEnd) - std::max(start, otherStart));
}

RepeatEdgeAlignment summarizeAlignment(const FullRead& read, int64_t repeatStart, int64_t repeatEnd, int motifLength)
{
    RepeatEdgeAlignment summary;
    int64_t referencePosition = read.s.pos;
    int readPosition = 0;
    bool reachedAlignedPart = false;
    for (const uint32_t cigarOperation : read.s.cigar)
    {
        const int operation = cigarOperation & BAM_CIGAR_MASK;
        const int length = cigarOperation >> BAM_CIGAR_SHIFT;
        switch (operation)
        {
        case BAM_CSOFT_CLIP:
            (reachedAlignedPart ? summary.rightClipLength : summary.leftClipLength) += length;
            readPosition += length;
            break;
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
        {
            reachedAlignedPart = true;
            const int64_t operationEnd = referencePosition + length;
            summary.flankBasesBefore += countOverlap(referencePosition, operationEnd, INT64_MIN, repeatStart);
            summary.flankBasesAfter += countOverlap(referencePosition, operationEnd, repeatEnd, INT64_MAX);
            summary.lengthForAssignment += countOverlap(referencePosition, operationEnd, repeatStart, repeatEnd);
            if (summary.readPositionAtStart == -1 && referencePosition <= repeatStart && repeatStart < operationEnd)
            {
                summary.readPositionAtStart = readPosition + (repeatStart - referencePosition);
                summary.referencePositionAtStart = repeatStart;
            }
            if (summary.readPositionAtEnd == -1 && referencePosition <= repeatEnd && repeatEnd <= operationEnd)
            {
                summary.readPositionAtEnd = readPosition + (repeatEnd - referencePosition);
            }
            referencePosition = operationEnd;
            readPosition += length;
            break;
        }
        case BAM_CINS:
            // A whole-motif insertion is a change in repeat length. As in the fast path, it counts toward the
            // assignment length when it lies in [S - k - 1, E + k], and it belongs to the tract when it sits
            // exactly on an edge. An insertion the aligner placed a few bases further out is separated from the
            // tract by aligned flank bases; taking it in as well is left to a later milestone.
            if (length % motifLength == 0)
            {
                if (referencePosition >= repeatStart - motifLength - 1 && referencePosition <= repeatEnd + motifLength)
                {
                    summary.lengthForAssignment += length;
                }
                if (summary.readPositionAtStart == -1 && referencePosition == repeatStart)
                {
                    summary.readPositionAtStart = readPosition;
                    summary.referencePositionAtStart = repeatStart - length;
                    summary.startsWithEdgeInsertion = true;
                }
                if (referencePosition == repeatEnd && summary.readPositionAtEnd == readPosition)
                {
                    summary.readPositionAtEnd += length;
                    summary.endsWithEdgeInsertion = true;
                }
            }
            readPosition += length;
            break;
        case BAM_CDEL:
        case BAM_CREF_SKIP:
        {
            reachedAlignedPart = true;
            const int64_t operationEnd = referencePosition + length;
            if (summary.readPositionAtStart == -1 && referencePosition <= repeatStart && repeatStart < operationEnd)
            {
                summary.readPositionAtStart = readPosition;
                summary.referencePositionAtStart = operationEnd;
            }
            if (summary.readPositionAtEnd == -1 && referencePosition <= repeatEnd && repeatEnd < operationEnd)
            {
                summary.readPositionAtEnd = readPosition;
            }
            referencePosition = operationEnd;
            break;
        }
        default:
            break;
        }
    }
    summary.referenceStart = read.s.pos;
    summary.referenceEnd = referencePosition;
    return summary;
}

int64_t computeAlignedEnd(const FullRead& read)
{
    int64_t end = read.s.pos;
    for (const uint32_t cigarOperation : read.s.cigar)
    {
        const int operation = cigarOperation & BAM_CIGAR_MASK;
        if (operation == BAM_CMATCH || operation == BAM_CEQUAL || operation == BAM_CDIFF || operation == BAM_CDEL
            || operation == BAM_CREF_SKIP)
        {
            end += cigarOperation >> BAM_CIGAR_SHIFT;
        }
    }
    return end;
}

bool isUsableAlignment(const FullRead& read)
{
    return !read.s.isSecondaryAlignment && !read.s.isSupplementaryAlignment;
}

// The fast path's rule: ignore reads with MAPQ <= 3 at loci whose reads usually map well.
bool hasUnusuallyLowMapq(const FullRead& read, double averageMapq) { return read.s.mapq <= 3 && averageMapq >= 20; }

// Frame offset, in the tract, of the first substring of a tract whose first base stands for `referencePosition`,
// given that units start at S + f + j * k in the reference.
int computeOffsetFromReferenceFrame(int64_t referencePosition, const MotifCompositionLocus& locus, int frameOffset)
{
    const int64_t motifLength = locus.catalogMotif.size();
    const int64_t offset = (locus.referenceRepeatStart + frameOffset - referencePosition) % motifLength;
    return static_cast<int>(offset < 0 ? offset + motifLength : offset);
}

// A 12-mer of reference flank close to the repeat that differs from the repeat, used to find where a read's repeat
// sequence ends and the flank begins.
struct FlankAnchor
{
    string sequence;
    int distanceFromRepeat; // flank bases between the repeat's edge and the anchor
};

const int kFlankAnchorLength = 12;

// True if the anchor matches the sequence at `position` with at most one mismatch.
bool matchesAnchorAt(const string& sequence, int position, const FlankAnchor& anchor)
{
    int mismatches = 0;
    for (int index = 0; index != kFlankAnchorLength && mismatches <= 1; ++index)
    {
        mismatches += !basesMatch(sequence[position + index], anchor.sequence[index]);
    }
    return mismatches <= 1;
}

// Picks the anchor from the reference flank next to the repeat (`flank` is read away from the repeat when
// isRightFlank, toward it otherwise): the 12-mer closest to the repeat, within its first 30 bases, that has at least
// 3 mismatches against every rotation of the catalog motif and of the given known motifs. Requiring the anchor to
// differ from the repeat keeps it from matching inside real repeat sequence.
boost::optional<FlankAnchor> findFlankAnchor(
    const string& flank, bool isRightFlank, const string& catalogMotif, const vector<string>& knownMotifs,
    const string& referenceRepeat)
{
    const int kMaxAnchorDistance = 30 - kFlankAnchorLength;
    const int kMinMismatchesVsRepeat = 3;
    vector<const string*> motifs = { &catalogMotif };
    for (const string& knownMotif : knownMotifs)
    {
        motifs.push_back(&knownMotif);
    }
    const int flankLength = flank.size();
    for (int distance = 0; distance <= std::min(kMaxAnchorDistance, flankLength - kFlankAnchorLength); ++distance)
    {
        const int start = isRightFlank ? distance : flankLength - kFlankAnchorLength - distance;
        const string window = flank.substr(start, kFlankAnchorLength);
        if (window.find('N') != string::npos)
        {
            continue;
        }
        bool differsFromRepeat = true;
        for (const string* motif : motifs)
        {
            const int motifLength = motif->size();
            for (int rotation = 0; rotation != motifLength && differsFromRepeat; ++rotation)
            {
                int mismatches = 0;
                for (int index = 0; index != kFlankAnchorLength; ++index)
                {
                    mismatches += !graphtools::checkIfReferenceBaseMatchesQueryBase(
                        (*motif)[(rotation + index) % motifLength], window[index]);
                }
                differsFromRepeat = mismatches >= kMinMismatchesVsRepeat;
            }
            if (!differsFromRepeat)
            {
                break;
            }
        }
        if (!differsFromRepeat)
        {
            continue;
        }
        // Nor may it match the reference repeat itself, which would cut reads of the reference allele.
        const FlankAnchor anchor { window, distance };
        bool matchesReferenceRepeat = false;
        for (int position = 0; position + kFlankAnchorLength <= static_cast<int>(referenceRepeat.size()); ++position)
        {
            if (matchesAnchorAt(referenceRepeat, position, anchor))
            {
                matchesReferenceRepeat = true;
                break;
            }
        }
        if (!matchesReferenceRepeat)
        {
            return anchor;
        }
    }
    return boost::none;
}

// New end of the tract [start, end) of `sequence`: where the right flank begins, if its anchor is found.
int cutAtRightAnchor(const string& sequence, int start, int end, const FlankAnchor& anchor)
{
    for (int position = start; position + kFlankAnchorLength <= end; ++position)
    {
        if (matchesAnchorAt(sequence, position, anchor))
        {
            return std::max(start, position - anchor.distanceFromRepeat);
        }
    }
    return end;
}

// New start of the tract [start, end) of `sequence`: just past the left flank, if its anchor is found.
int cutAtLeftAnchor(const string& sequence, int start, int end, const FlankAnchor& anchor)
{
    for (int position = end - kFlankAnchorLength; position >= start; --position)
    {
        if (matchesAnchorAt(sequence, position, anchor))
        {
            return std::min(end, position + kFlankAnchorLength + anchor.distanceFromRepeat);
        }
    }
    return start;
}

// ------------------------------------------------------------------------------------------------------------
// Acceptance of new motifs
// ------------------------------------------------------------------------------------------------------------

struct MotifStats
{
    int occurrences = 0;
    vector<int> highQualityBaseCounts; // for each position in the motif, observations with a high-quality base there
};

struct SequenceSubstringLocation
{
    int tractIndex;
    int offsetWithinRepeatTract;
};

// Known motifs ordered most common first (ties: motif sequence), which is the order the shift search uses.
vector<string> orderByOccurrences(const std::unordered_map<string, MotifStats>& stats)
{
    vector<string> motifs;
    motifs.reserve(stats.size());
    for (const auto& motifAndStats : stats)
    {
        motifs.push_back(motifAndStats.first);
    }
    std::sort(
        motifs.begin(), motifs.end(),
        [&](const string& left, const string& right)
        {
            const int leftCount = stats.at(left).occurrences;
            const int rightCount = stats.at(right).occurrences;
            return leftCount != rightCount ? leftCount > rightCount : left < right;
        });
    return motifs;
}

bool isMoreCommon(int occurrences, const string& motif, int otherOccurrences, const string& otherMotif)
{
    return otherOccurrences != occurrences ? otherOccurrences > occurrences : otherMotif < motif;
}

// Positions where `motif` differs from its nearest more common motifs: every one that differs from it at exactly
// one position, or, if there are none, those with the fewest differences. A base there has to be high quality
// for an observation of the motif to count, since a low-quality error at exactly that position would turn the
// more common motif into this one.
vector<int> findDistinguishingPositions(
    const string& motif, int occurrences, const std::unordered_map<string, MotifStats>& knownStats)
{
    int fewestDifferences = INT_MAX;
    vector<const string*> nearest;
    for (const auto& otherAndStats : knownStats)
    {
        const string& other = otherAndStats.first;
        if (other == motif || !isMoreCommon(occurrences, motif, otherAndStats.second.occurrences, other))
        {
            continue;
        }
        const int differences = countDifferences(motif, other);
        if (differences < fewestDifferences)
        {
            fewestDifferences = differences;
            nearest.clear();
        }
        if (differences == fewestDifferences)
        {
            nearest.push_back(&other);
        }
    }
    vector<char> isDistinguishing(motif.size(), false);
    for (const string* other : nearest)
    {
        for (size_t index = 0; index != motif.size(); ++index)
        {
            isDistinguishing[index] = isDistinguishing[index] || motif[index] != (*other)[index];
        }
    }
    vector<int> positions;
    for (size_t index = 0; index != motif.size(); ++index)
    {
        if (isDistinguishing[index])
        {
            positions.push_back(index);
        }
    }
    return positions;
}

bool hasHighQualityBasesAt(const string& bases, int start, const vector<int>& positions)
{
    return std::all_of(
        positions.begin(), positions.end(), [&](int position) { return isHighQuality(bases[start + position]); });
}

struct AcceptanceTest
{
    double errorRate;
    double pValueThreshold; // alpha_locus / n
    int minReadPairs; // max(minReadPairs, ceil(minReadPairFraction * N))
};

// The occurrences of a motif seen at `locations`, and for each of its positions how many of them have a high-quality
// base there.
MotifStats tallyLocations(
    const vector<SequenceSubstringLocation>& locations, const vector<RepeatTract>& tracts, int motifLength)
{
    MotifStats stats;
    stats.occurrences = locations.size();
    stats.highQualityBaseCounts.assign(motifLength, 0);
    for (const SequenceSubstringLocation& location : locations)
    {
        const string& bases = tracts[location.tractIndex].bases;
        for (int index = 0; index != motifLength; ++index)
        {
            stats.highQualityBaseCounts[index] += isHighQuality(bases[location.offsetWithinRepeatTract + index]);
        }
    }
    return stats;
}

// The number of distinct read pairs that show `motif` at `locations` with high-quality bases at every position where it
// differs from its nearest more common known motifs.
int countSupportingReadPairs(
    const string& motif, const vector<SequenceSubstringLocation>& locations,
    const std::unordered_map<string, MotifStats>& knownStats, const vector<RepeatTract>& tracts)
{
    const vector<int> distinguishingPositions = findDistinguishingPositions(motif, locations.size(), knownStats);
    std::unordered_set<int> readPairs;
    for (const SequenceSubstringLocation& location : locations)
    {
        const RepeatTract& tract = tracts[location.tractIndex];
        if (hasHighQualityBasesAt(tract.bases, location.offsetWithinRepeatTract, distinguishingPositions))
        {
            readPairs.insert(tract.readPairIndex);
        }
    }
    return readPairs.size();
}

// Tests candidates in decreasing order of occurrences and adds the accepted ones to knownStats.
vector<string> acceptNewMotifs(
    const std::unordered_map<string, vector<SequenceSubstringLocation>>& candidates,
    std::unordered_map<string, MotifStats>& knownStats, const vector<RepeatTract>& tracts, const AcceptanceTest& test)
{
    vector<const string*> order;
    for (const auto& motifAndLocations : candidates)
    {
        if (static_cast<int>(motifAndLocations.second.size()) >= test.minReadPairs
            && knownStats.count(motifAndLocations.first) == 0)
        {
            order.push_back(&motifAndLocations.first);
        }
    }
    std::sort(
        order.begin(), order.end(),
        [&](const string* left, const string* right)
        {
            const size_t leftCount = candidates.at(*left).size();
            const size_t rightCount = candidates.at(*right).size();
            return leftCount != rightCount ? leftCount > rightCount : *left < *right;
        });

    vector<string> accepted;
    for (const string* candidate : order)
    {
        const vector<SequenceSubstringLocation>& locations = candidates.at(*candidate);
        const int count = countSupportingReadPairs(*candidate, locations, knownStats, tracts);
        if (count < test.minReadPairs)
        {
            continue;
        }

        // Expected number of read pairs showing this candidate through a single high-quality sequencing error of a
        // known motif one base away from it.
        double expectedErrors = 0;
        for (const auto& knownAndStats : knownStats)
        {
            const string& known = knownAndStats.first;
            if (countDifferences(*candidate, known) != 1)
            {
                continue;
            }
            for (size_t index = 0; index != known.size(); ++index)
            {
                if (known[index] != (*candidate)[index])
                {
                    expectedErrors += test.errorRate / 3 * knownAndStats.second.highQualityBaseCounts[index];
                }
            }
        }
        if (computePoissonUpperTail(expectedErrors, count) >= test.pValueThreshold)
        {
            continue;
        }

        knownStats.emplace(*candidate, tallyLocations(locations, tracts, candidate->size()));
        accepted.push_back(*candidate);
    }
    return accepted;
}

// ------------------------------------------------------------------------------------------------------------
// Counting
// ------------------------------------------------------------------------------------------------------------

struct Tally
{
    int occurrences = 0;
    int reads = 0;
    int lastTractIndex = -1;

    void add(int tractIndex)
    {
        ++occurrences;
        if (lastTractIndex != tractIndex)
        {
            ++reads;
            lastTractIndex = tractIndex;
        }
    }
};

struct GroupTallies
{
    explicit GroupTallies(size_t motifCount)
        : motifs(motifCount)
    {
    }

    vector<Tally> motifs; // indexed by position in the final motif list
    std::map<std::pair<int, int>, Tally> motifPairs;
};

MotifCompositionCounts encodeCounts(const GroupTallies& tallies, const vector<int>& motifIds)
{
    MotifCompositionCounts counts;
    for (size_t index = 0; index != tallies.motifs.size(); ++index)
    {
        const Tally& tally = tallies.motifs[index];
        if (tally.occurrences > 0)
        {
            counts.motifs[motifIds[index]] = { tally.occurrences, tally.reads };
        }
    }
    for (const auto& pairAndTally : tallies.motifPairs)
    {
        const auto& pair = pairAndTally.first;
        counts.motifPairs[{ motifIds[pair.first], motifIds[pair.second] }]
            = { pairAndTally.second.occurrences, pairAndTally.second.reads };
    }
    return counts;
}

// 0 = not assigned, 1 = shorter allele, 2 = longer allele.
int assignToAllele(const RepeatTract& tract, int shortAllele, int longAllele, int motifLength, int readLength)
{
    switch (tract.kind)
    {
    case ReadKind::kSpanning:
    {
        const int units = tract.lengthForAssignment / motifLength;
        const int distanceToShort = std::abs(units - shortAllele);
        const int distanceToLong = std::abs(units - longAllele);
        if (distanceToShort < distanceToLong && distanceToShort <= std::max(1.0, 0.1 * shortAllele))
        {
            return 1;
        }
        if (distanceToLong < distanceToShort && distanceToLong <= std::max(1.0, 0.1 * longAllele))
        {
            return 2;
        }
        return 0;
    }
    case ReadKind::kFlanking:
        return tract.lengthForAssignment / motifLength > shortAllele + std::max(1.0, 0.1 * shortAllele) ? 2 : 0;
    case ReadKind::kInsideRepeat:
    case ReadKind::kInrepeat:
        return (longAllele * motifLength >= readLength && shortAllele * motifLength < readLength) ? 2 : 0;
    }
    return 0;
}

} // namespace

double computePeriodScore(const string& sequence, int lag)
{
    if (lag <= 0 || static_cast<int>(sequence.size()) <= lag)
    {
        return 0.0;
    }
    int matches = 0;
    for (size_t index = lag; index != sequence.size(); ++index)
    {
        matches += basesMatch(sequence[index], sequence[index - lag]);
    }
    return static_cast<double>(matches) / (sequence.size() - lag);
}

bool passesPeriodTest(const string& sequence, int motifLength, double minScore)
{
    const double score = computePeriodScore(sequence, motifLength);
    if (score < minScore)
    {
        return false;
    }
    // The small tolerance keeps a margin of exactly 0.1 from failing on rounding.
    const double kMinMarginOverDivisors = 0.1 - 1e-9;
    for (int divisor = 1; divisor < motifLength; ++divisor)
    {
        if (motifLength % divisor == 0 && score - computePeriodScore(sequence, divisor) < kMinMarginOverDivisors)
        {
            return false;
        }
    }
    return true;
}

int computeFrameOffset(const string& sequence, const string& motif)
{
    return computeFrameOffsetAgainst(sequence, motif, {});
}

int computeSoftClipBasesToKeep(
    const string& sequence, int tractStart, int clipStart, int clipEnd, int tractEnd, bool clipOnRight, int motifLength,
    double minScore)
{
    const int clipLength = clipEnd - clipStart;
    if (clipLength <= 0)
    {
        return 0;
    }
    const int windowSize = std::max(2 * motifLength, 12);
    vector<char> window(windowSize, 0);
    int windowHead = 0;
    int comparisons = 0;
    int matchesInWindow = 0;
    int basesThroughLastMatch = 0;
    for (int step = 0; step != clipLength; ++step)
    {
        const int position = clipOnRight ? clipStart + step : clipEnd - 1 - step;
        const int earlierPosition = clipOnRight ? position - motifLength : position + motifLength;
        const bool isComparable = clipOnRight ? earlierPosition >= tractStart : earlierPosition < tractEnd;
        if (!isComparable)
        {
            continue;
        }
        const bool isMatch = basesMatch(sequence[position], sequence[earlierPosition]);
        if (comparisons == windowSize)
        {
            matchesInWindow -= window[windowHead];
        }
        else
        {
            ++comparisons;
        }
        window[windowHead] = isMatch;
        matchesInWindow += isMatch;
        windowHead = (windowHead + 1) % windowSize;
        if (isMatch)
        {
            basesThroughLastMatch = step + 1;
        }
        if (comparisons == windowSize && matchesInWindow < minScore * windowSize)
        {
            return basesThroughLastMatch;
        }
    }
    // A clip too short to fill the window: drop trailing bases that do not continue the repeat.
    if (comparisons > 0 && matchesInWindow < minScore * comparisons)
    {
        return basesThroughLastMatch;
    }
    return clipLength;
}

vector<SequenceSubstring> splitIntoMotifs(
    const string& tract, const string& catalogMotif, const vector<string>& knownMotifs, int startOffset,
    bool startsAtRepeatEdge, bool endsAtRepeatEdge, const boost::optional<string>& referenceEndPartial)
{
    vector<SequenceSubstring> sequenceSubstrings;
    splitTract(
        tract, toUpperSequence(tract), catalogMotif, MotifList(knownMotifs), startOffset, startsAtRepeatEdge,
        endsAtRepeatEdge, referenceEndPartial, sequenceSubstrings);
    return sequenceSubstrings;
}

double computePoissonUpperTail(double mean, int count)
{
    if (count <= 0)
    {
        return 1.0;
    }
    if (mean <= 0.0)
    {
        return 0.0;
    }
    // Sum the tail from `count` up; terms shrink geometrically once past the mean.
    double term = std::exp(-mean + count * std::log(mean) - std::lgamma(count + 1.0));
    double sum = 0.0;
    for (int value = count; value < count + 100000; ++value)
    {
        sum += term;
        term *= mean / (value + 1);
        if (value > mean && term < sum * 1e-15)
        {
            break;
        }
    }
    return std::min(1.0, sum);
}

} // namespace motifcomposition

using namespace motifcomposition;

vector<string> selectCatalogKnownMotifs(const vector<string>& knownMotifs, int motifLength, const string& locusId)
{
    vector<string> selected;
    std::unordered_map<string, string> selectedRotations; // canonical rotation -> the entry as written in the catalog
    for (const string& knownMotif : knownMotifs)
    {
        string motif(knownMotif);
        std::transform(
            motif.begin(), motif.end(), motif.begin(),
            [](unsigned char base) { return static_cast<char>(std::toupper(base)); });
        string problem;
        if (static_cast<int>(motif.size()) != motifLength)
        {
            problem = "its length differs from the LocusStructure motif's";
        }
        else if (motif.find_first_not_of("ACGT") != string::npos)
        {
            problem = "it has bases other than A, C, G and T";
        }
        else
        {
            const auto inserted = selectedRotations.emplace(computeCanonicalRotation(motif), knownMotif);
            if (!inserted.second)
            {
                problem = "it's already listed as '" + inserted.first->second + "'";
            }
        }
        if (!problem.empty())
        {
            spdlog::warn("Skipping KnownMotifs entry '{}' of locus {} because {}", knownMotif, locusId, problem);
            continue;
        }
        selected.push_back(std::move(motif));
    }
    return selected;
}

boost::optional<MotifComposition> computeMotifComposition(
    const MotifCompositionLocus& locus, int typicalReadLength, int flankLength,
    const vector<const FullReadPair*>& readPairs, const boost::optional<RepeatGenotype>& genotype,
    bool onlyLociWithNonRefMotifs)
{
    const string& catalogMotif = locus.catalogMotif;
    const int motifLength = catalogMotif.size();
    const int64_t repeatStart = locus.referenceRepeatStart;
    const int64_t repeatEnd = locus.referenceRepeatEnd;
    const bool catalogMotifIsConcrete = graphtools::checkIfNucleotideReferenceSequence(catalogMotif);

    // --- Reference motifs: parse the reference repeat, both of whose ends are trusted.
    const int frameOffset = computeFrameOffset(locus.referenceRepeatSequence, catalogMotif);
    std::unordered_map<string, MotifStats> knownStats;
    {
        const vector<string> seedMotifs = catalogMotifIsConcrete ? vector<string> { catalogMotif } : vector<string> {};
        vector<SequenceSubstring> sequenceSubstrings;
        splitTract(
            locus.referenceRepeatSequence, locus.referenceRepeatSequence, catalogMotif, MotifList(seedMotifs),
            frameOffset, true, true, boost::none, sequenceSubstrings);
        for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
        {
            if (sequenceSubstring.type != SequenceSubstringType::kGap)
            {
                ++knownStats[locus.referenceRepeatSequence.substr(
                                 sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length)]
                      .occurrences;
            }
        }
        if (catalogMotifIsConcrete)
        {
            knownStats[catalogMotif];
        }
    }
    // The reference motif set R'. Seeds start with the reference's counts only for ordering; the counts used by the
    // acceptance test come from the reads.
    std::unordered_set<string> referenceMotifs;
    const vector<string> seedOrder = orderByOccurrences(knownStats);
    for (auto& motifAndStats : knownStats)
    {
        referenceMotifs.insert(motifAndStats.first);
        motifAndStats.second.occurrences = 0;
        motifAndStats.second.highQualityBaseCounts.assign(motifLength, 0);
    }

    const boost::optional<FlankAnchor> leftAnchor
        = findFlankAnchor(locus.leftFlankSequence, false, catalogMotif, seedOrder, locus.referenceRepeatSequence);
    const boost::optional<FlankAnchor> rightAnchor
        = findFlankAnchor(locus.rightFlankSequence, true, catalogMotif, seedOrder, locus.referenceRepeatSequence);

    // The partial unit the reference has before E, which a read's tract ending at E should also have.
    const int repeatLength = locus.referenceRepeatSequence.size();
    const boost::optional<string> referenceEndPartial = locus.referenceRepeatSequence.substr(
        repeatLength - std::max(0, repeatLength - frameOffset) % motifLength);

    // --- Collect the tracts: reads mapped to the repeat first, then in-repeat reads.
    double averageMapq = 0;
    int mappedReadCount = 0;
    for (const FullReadPair* readPair : readPairs)
    {
        for (const std::optional<FullRead>* mate : { &readPair->firstMate, &readPair->secondMate })
        {
            if (*mate && (*mate)->s.isMapped)
            {
                averageMapq += (*mate)->s.mapq;
                ++mappedReadCount;
            }
        }
    }
    if (mappedReadCount > 0)
    {
        averageMapq /= mappedReadCount;
    }

    const int64_t anchorDistance = locus.meanFragmentLength > 0 ? locus.meanFragmentLength : flankLength;
    // A mate anchored in a flank of the repeat and facing it, so the fragment's other read lies toward the repeat.
    auto isAnchoredFacingRepeat = [&](const FullRead& mate)
    {
        if (!mate.s.isMapped || mate.s.chromId != locus.contigIndex || !isUsableAlignment(mate))
        {
            return false;
        }
        const int64_t mateStart = mate.s.pos;
        const int64_t mateEnd = computeAlignedEnd(mate);
        if (!mate.r.isReversed())
        {
            return mateStart < repeatStart && mateEnd <= repeatEnd && repeatStart - mateStart <= anchorDistance;
        }
        return mateEnd > repeatEnd && mateStart >= repeatStart && mateEnd - repeatEnd <= anchorDistance;
    };

    vector<RepeatTract> tracts;
    vector<RepeatTract> inrepeatTracts;
    for (size_t pairIndex = 0; pairIndex != readPairs.size(); ++pairIndex)
    {
        const FullReadPair& readPair = *readPairs[pairIndex];
        const FullRead* mates[2] = { readPair.firstMate ? &*readPair.firstMate : nullptr,
                                     readPair.secondMate ? &*readPair.secondMate : nullptr };
        for (int mateIndex = 0; mateIndex != 2; ++mateIndex)
        {
            const FullRead* read = mates[mateIndex];
            const FullRead* mate = mates[1 - mateIndex];
            if (read == nullptr || !isUsableAlignment(*read))
            {
                continue;
            }
            const string& bases = read->r.sequence();
            const int readLength = bases.size();
            const bool mateIsAnchored = mate != nullptr && isAnchoredFacingRepeat(*mate);

            RepeatTract tract;
            tract.readPairIndex = static_cast<int>(pairIndex);
            tract.lengthForAssignment = 0;
            tract.startsAtRepeatEdge = false;
            tract.endsAtRepeatEdge = false;
            tract.startOffset = -1;
            int tractStart = 0;
            int tractEnd = 0;
            // True for a read that overlaps the repeat or is soft-clipped toward it; such a read is never an
            // in-repeat read, even if it holds no repeat bases.
            bool isPlacedAtRepeat = false;

            if (read->s.isMapped && read->s.chromId == locus.contigIndex && !read->s.cigar.empty())
            {
                const RepeatEdgeAlignment alignment = summarizeAlignment(*read, repeatStart, repeatEnd, motifLength);
                const bool overlapsRepeat
                    = alignment.referenceStart < repeatEnd && alignment.referenceEnd > repeatStart;
                if (overlapsRepeat)
                {
                    isPlacedAtRepeat = true;
                    const bool leftClipFacesRepeat
                        = alignment.leftClipLength > 0 && alignment.referenceStart >= repeatStart;
                    const bool rightClipFacesRepeat
                        = alignment.rightClipLength > 0 && alignment.referenceEnd <= repeatEnd;
                    if (!leftClipFacesRepeat && alignment.readPositionAtStart >= 0)
                    {
                        tractStart = alignment.readPositionAtStart;
                        tract.startsAtRepeatEdge = !alignment.startsWithEdgeInsertion;
                        tract.startOffset
                            = computeOffsetFromReferenceFrame(alignment.referencePositionAtStart, locus, frameOffset);
                    }
                    else
                    {
                        tractStart = alignment.leftClipLength;
                    }
                    if (!rightClipFacesRepeat && alignment.readPositionAtEnd >= 0)
                    {
                        tractEnd = alignment.readPositionAtEnd;
                        tract.endsAtRepeatEdge = !alignment.endsWithEdgeInsertion;
                    }
                    else
                    {
                        tractEnd = readLength - alignment.rightClipLength;
                    }
                    if (leftClipFacesRepeat)
                    {
                        tractStart -= computeSoftClipBasesToKeep(
                            bases, 0, 0, alignment.leftClipLength, tractEnd, false, motifLength,
                            kMotifCompositionMinSoftClipPeriodScore);
                    }
                    if (rightClipFacesRepeat)
                    {
                        tractEnd += computeSoftClipBasesToKeep(
                            bases, tractStart, tractEnd, readLength, 0, true, motifLength,
                            kMotifCompositionMinSoftClipPeriodScore);
                    }

                    const bool isSpanning = alignment.flankBasesBefore >= kMotifCompositionMinFlankBases
                        && alignment.flankBasesAfter >= kMotifCompositionMinFlankBases;
                    const bool isInside
                        = alignment.referenceStart >= repeatStart && alignment.referenceEnd <= repeatEnd;
                    tract.kind
                        = isSpanning ? ReadKind::kSpanning : (isInside ? ReadKind::kInsideRepeat : ReadKind::kFlanking);
                    tract.lengthForAssignment
                        = isSpanning ? alignment.lengthForAssignment : std::max(0, tractEnd - tractStart);
                }
                else if (
                    alignment.rightClipLength > 0 && alignment.referenceEnd <= repeatStart
                    && alignment.referenceEnd > repeatStart - motifLength)
                {
                    // Aligned part ends in the left flank within k bases of the repeat, soft clip facing it: skip
                    // the clip bases that stand for the rest of the flank, so the tract starts at S.
                    isPlacedAtRepeat = true;
                    tractStart = readLength - alignment.rightClipLength
                        + static_cast<int>(repeatStart - alignment.referenceEnd);
                    if (tractStart < readLength)
                    {
                        tractEnd = tractStart
                            + computeSoftClipBasesToKeep(
                                       bases, tractStart, tractStart, readLength, 0, true, motifLength,
                                       kMotifCompositionMinSoftClipPeriodScore);
                        tract.startsAtRepeatEdge = true;
                        tract.startOffset = frameOffset;
                        tract.kind = ReadKind::kFlanking;
                        tract.lengthForAssignment = tractEnd - tractStart;
                    }
                }
                else if (
                    alignment.leftClipLength > 0 && alignment.referenceStart >= repeatEnd
                    && alignment.referenceStart < repeatEnd + motifLength)
                {
                    // The mirror image: aligned part starts in the right flank within k bases of E.
                    isPlacedAtRepeat = true;
                    tractEnd = alignment.leftClipLength - static_cast<int>(alignment.referenceStart - repeatEnd);
                    if (tractEnd > 0)
                    {
                        tractStart = tractEnd
                            - computeSoftClipBasesToKeep(
                                         bases, 0, 0, tractEnd, tractEnd, false, motifLength,
                                         kMotifCompositionMinSoftClipPeriodScore);
                        tract.endsAtRepeatEdge = true;
                        tract.kind = ReadKind::kFlanking;
                        tract.lengthForAssignment = tractEnd - tractStart;
                    }
                }
            }

            if (isPlacedAtRepeat)
            {
                // A tract can run past the repeat into the flank: a soft clip that continues the repeat's period,
                // or flank bases the aligner placed on the repeat (as in reads from an allele shorter than the
                // reference). Cut it where the flank's anchor appears; the cut is not a trusted edge.
                if (rightAnchor && tractEnd > tractStart)
                {
                    const int cutEnd = cutAtRightAnchor(bases, tractStart, tractEnd, *rightAnchor);
                    if (cutEnd != tractEnd)
                    {
                        tractEnd = cutEnd;
                        tract.endsAtRepeatEdge = false;
                    }
                }
                if (leftAnchor && tractEnd > tractStart)
                {
                    const int cutStart = cutAtLeftAnchor(bases, tractStart, tractEnd, *leftAnchor);
                    if (cutStart != tractStart)
                    {
                        tractStart = cutStart;
                        tract.startsAtRepeatEdge = false;
                        tract.startOffset = -1;
                    }
                }
                if (tract.kind == ReadKind::kFlanking)
                {
                    tract.lengthForAssignment = std::max(0, tractEnd - tractStart);
                }
                if (tractEnd <= tractStart)
                {
                    continue; // no repeat bases, for example an alignment that deletes the whole repeat
                }
                // Low MAPQ is judged on the read itself, except for reads mapped inside the repeat, whose placement
                // is ambiguous by nature: for those it is judged on an anchored mate, if there is one.
                const FullRead* mapqRead
                    = tract.kind == ReadKind::kInsideRepeat ? (mateIsAnchored ? mate : nullptr) : read;
                if (mapqRead != nullptr && hasUnusuallyLowMapq(*mapqRead, averageMapq))
                {
                    continue;
                }
                tract.bases = bases.substr(tractStart, tractEnd - tractStart);
                tracts.push_back(std::move(tract));
                continue;
            }

            // A read that BWA placed confidently near the repeat, but not in it, came from where it was placed, for
            // example a neighbouring repeat of the same period or a repetitive flank; it is not an in-repeat read of
            // this repeat. (Telling apart the in-repeat reads of an expansion that BWA placed in a neighbouring repeat
            // of the same motif needs a comparison between neighbouring repeats, not yet implemented.)
            if (read->s.isMapped && read->s.chromId == locus.contigIndex && read->s.mapq > 3
                && read->s.pos < repeatEnd + flankLength && computeAlignedEnd(*read) > repeatStart - flankLength)
            {
                continue;
            }

            // In-repeat read: made of repeat and placed anywhere else. Its orientation comes from the
            // anchored mate: in a forward-reverse pair it comes from the strand opposite the mate.
            if (!mateIsAnchored || hasUnusuallyLowMapq(*mate, averageMapq)
                || !passesPeriodTest(bases, motifLength, kMotifCompositionMinInrepeatReadPeriodScore))
            {
                continue;
            }
            tract.kind = ReadKind::kInrepeat;
            tract.bases = mate->r.isReversed() == read->r.isReversed() ? graphtools::reverseComplement(bases) : bases;
            {
                // An in-repeat read can still hold some flank at either end.
                int start = 0;
                int end = tract.bases.size();
                if (rightAnchor)
                {
                    end = cutAtRightAnchor(tract.bases, start, end, *rightAnchor);
                }
                if (leftAnchor && end > start)
                {
                    start = cutAtLeftAnchor(tract.bases, start, end, *leftAnchor);
                }
                if (end <= start)
                {
                    continue;
                }
                tract.bases = tract.bases.substr(start, end - start);
            }
            inrepeatTracts.push_back(std::move(tract));
        }
    }
    const size_t mappedTractCount = tracts.size();
    for (RepeatTract& tract : inrepeatTracts)
    {
        tracts.push_back(std::move(tract));
    }
    inrepeatTracts.clear();

    const bool isHeterozygousWithDistinctAlleles = genotype && genotype->numAlleles() == 2
        && genotype->longAlleleSizeInUnits() - genotype->shortAlleleSizeInUnits() >= 2;

    // --- Parse with a given motif list and collect candidates.
    vector<string> upperTracts(tracts.size());
    for (size_t index = 0; index != tracts.size(); ++index)
    {
        upperTracts[index] = toUpperSequence(tracts[index].bases);
    }
    auto startOffsetOf = [&](size_t tractIndex, const vector<string>& motifs)
    {
        const RepeatTract& tract = tracts[tractIndex];
        return tract.startOffset >= 0 ? tract.startOffset
                                      : computeFrameOffsetAgainst(upperTracts[tractIndex], catalogMotif, motifs);
    };

    std::unordered_map<string, vector<SequenceSubstringLocation>> candidates;
    vector<char> readPairHasSequenceSubstrings(readPairs.size(), false);
    int readPairsWithSequenceSubstrings = 0;
    vector<SequenceSubstring> sequenceSubstrings;
    auto parseForDiscovery = [&](size_t firstTract, size_t lastTract, const MotifList& motifList)
    {
        for (size_t tractIndex = firstTract; tractIndex != lastTract; ++tractIndex)
        {
            const RepeatTract& tract = tracts[tractIndex];
            splitTract(
                tract.bases, upperTracts[tractIndex], catalogMotif, motifList,
                startOffsetOf(tractIndex, motifList.motifs()), tract.startsAtRepeatEdge, tract.endsAtRepeatEdge,
                referenceEndPartial, sequenceSubstrings);
            for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
            {
                if (sequenceSubstring.type == SequenceSubstringType::kGap)
                {
                    continue;
                }
                if (!readPairHasSequenceSubstrings[tract.readPairIndex])
                {
                    readPairHasSequenceSubstrings[tract.readPairIndex] = true;
                    ++readPairsWithSequenceSubstrings;
                }
                const string motif = upperTracts[tractIndex].substr(
                    sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length);
                if (sequenceSubstring.type == SequenceSubstringType::kKnownMotif)
                {
                    MotifStats& stats = knownStats[motif];
                    ++stats.occurrences;
                    for (int index = 0; index != motifLength; ++index)
                    {
                        stats.highQualityBaseCounts[index]
                            += isHighQuality(tract.bases[sequenceSubstring.offsetWithinRepeatTract + index]);
                    }
                }
                else
                {
                    candidates[motif].push_back(
                        { static_cast<int>(tractIndex), sequenceSubstring.offsetWithinRepeatTract });
                }
            }
        }
    };

    // The catalog's KnownMotifs are matched by rotation, since the reads, cut in line with the catalog motif, show
    // which rotation actually occurs. A listed motif is found once one of its rotations is known: a reference motif,
    // or a candidate trusted here. For a listed motif not yet found, the candidate among its rotations seen most often
    // (ties: alphabetically first) is trusted without the error test of acceptNewMotifs, provided it has the same
    // minimum read-pair support a new motif needs (so a single sequencing error does not report it); its other
    // rotations are frame shifts and stay ordinary candidates.
    std::unordered_set<string> catalogRotationsToFind;
    for (const string& motif : locus.catalogKnownMotifs)
    {
        catalogRotationsToFind.insert(computeCanonicalRotation(motif));
    }
    for (const string& motif : referenceMotifs)
    {
        catalogRotationsToFind.erase(computeCanonicalRotation(motif));
    }

    // Bonferroni correction over a family fixed before looking at the reads: the one-base variants of the reference
    // motifs and of the listed motifs that are not rotations of one.
    const int trustedMotifCount = std::max<int>(1, referenceMotifs.size() + catalogRotationsToFind.size());
    auto makeAcceptanceTest = [&]()
    {
        AcceptanceTest test;
        test.errorRate = kMotifCompositionBaseErrorRate;
        test.pValueThreshold = kMotifCompositionLocusFalsePositiveRate / (3.0 * motifLength * trustedMotifCount);
        test.minReadPairs = std::max(
            kMotifCompositionMinReadPairs,
            static_cast<int>(std::ceil(kMotifCompositionMinReadPairFraction * readPairsWithSequenceSubstrings)));
        return test;
    };

    std::unordered_set<string> catalogMotifsInReads;
    auto trustCatalogKnownCandidates = [&]()
    {
        std::map<string, string> best; // canonical rotation -> candidate
        for (const auto& motifAndLocations : candidates)
        {
            const string& motif = motifAndLocations.first;
            const string canonical = computeCanonicalRotation(motif);
            if (catalogRotationsToFind.count(canonical) == 0)
            {
                continue;
            }
            auto it = best.find(canonical);
            if (it == best.end())
            {
                best.emplace(canonical, motif);
                continue;
            }
            const size_t count = motifAndLocations.second.size();
            const size_t bestCount = candidates.at(it->second).size();
            if (count > bestCount || (count == bestCount && motif < it->second))
            {
                it->second = motif;
            }
        }
        const int minReadPairs = makeAcceptanceTest().minReadPairs;
        for (const auto& canonicalAndMotif : best)
        {
            const string& motif = canonicalAndMotif.second;
            if (countSupportingReadPairs(motif, candidates.at(motif), knownStats, tracts) < minReadPairs)
            {
                continue;
            }
            knownStats.emplace(motif, tallyLocations(candidates.at(motif), tracts, motifLength));
            catalogMotifsInReads.insert(motif);
            catalogRotationsToFind.erase(canonicalAndMotif.first);
            candidates.erase(motif);
        }
    };

    // Pass 1: reads mapped to the repeat, with the reference motifs as the known list.
    parseForDiscovery(0, mappedTractCount, MotifList(seedOrder));
    trustCatalogKnownCandidates();
    if (onlyLociWithNonRefMotifs && candidates.empty() && catalogMotifsInReads.empty()
        && mappedTractCount == tracts.size())
    {
        return boost::none; // every substring is a reference motif and there are no in-repeat reads
    }
    vector<string> newMotifs = acceptNewMotifs(candidates, knownStats, tracts, makeAcceptanceTest());

    // In-repeat reads, parsed with the motifs accepted so far; their candidates are pooled with pass 1's.
    if (mappedTractCount != tracts.size())
    {
        parseForDiscovery(mappedTractCount, tracts.size(), MotifList(orderByOccurrences(knownStats)));
        trustCatalogKnownCandidates();
        const vector<string> moreNewMotifs = acceptNewMotifs(candidates, knownStats, tracts, makeAcceptanceTest());
        newMotifs.insert(newMotifs.end(), moreNewMotifs.begin(), moreNewMotifs.end());
    }
    candidates.clear();

    // --- Final motif list: merge the rotations of each new motif into the one closest to the catalog
    // motif, and drop new motifs that are rotations of a reference motif or of a catalog known motif found in the
    // reads.
    {
        std::unordered_set<string> trustedRotations;
        for (const std::unordered_set<string>* motifs : { &referenceMotifs, &catalogMotifsInReads })
        {
            for (const string& motif : *motifs)
            {
                trustedRotations.insert(computeCanonicalRotation(motif));
            }
        }
        std::map<string, string> keptRotation; // canonical rotation -> kept new motif
        for (const string& motif : newMotifs)
        {
            const string canonical = computeCanonicalRotation(motif);
            if (trustedRotations.count(canonical) != 0)
            {
                continue;
            }
            auto kept = keptRotation.find(canonical);
            if (kept == keptRotation.end())
            {
                keptRotation.emplace(canonical, motif);
                continue;
            }
            const int mismatches = countMismatches(motif.data(), catalogMotif, motifLength);
            const int keptMismatches = countMismatches(kept->second.data(), catalogMotif, motifLength);
            if (mismatches < keptMismatches || (mismatches == keptMismatches && motif < kept->second))
            {
                kept->second = motif;
            }
        }
        std::unordered_set<string> keptNewMotifs;
        for (const auto& canonicalAndMotif : keptRotation)
        {
            keptNewMotifs.insert(canonicalAndMotif.second);
        }
        for (auto it = knownStats.begin(); it != knownStats.end();)
        {
            const bool keep = referenceMotifs.count(it->first) != 0 || catalogMotifsInReads.count(it->first) != 0
                || keptNewMotifs.count(it->first) != 0;
            it = keep ? std::next(it) : knownStats.erase(it);
        }
        if (onlyLociWithNonRefMotifs && keptNewMotifs.empty() && catalogMotifsInReads.empty())
        {
            return boost::none;
        }
    }

    const vector<string> finalMotifs = orderByOccurrences(knownStats);
    const MotifList finalMotifList(finalMotifs);
    std::unordered_map<string, int> finalMotifIndex;
    vector<vector<int>> distinguishingPositions(finalMotifs.size());
    for (size_t index = 0; index != finalMotifs.size(); ++index)
    {
        finalMotifIndex.emplace(finalMotifs[index], index);
        distinguishingPositions[index] = findDistinguishingPositions(
            finalMotifs[index], knownStats.at(finalMotifs[index]).occurrences, knownStats);
    }

    // --- Final parse and counting. A substring counts only if it is on the final list and its bases
    // are high quality where it differs from its nearest more common motifs; anything else becomes a gap, which
    // breaks pairs.
    GroupTallies locusTallies(finalMotifs.size());
    GroupTallies allele1Tallies(finalMotifs.size());
    GroupTallies allele2Tallies(finalMotifs.size());
    bool anyReadAssigned = false;
    for (size_t tractIndex = 0; tractIndex != tracts.size(); ++tractIndex)
    {
        const RepeatTract& tract = tracts[tractIndex];
        splitTract(
            tract.bases, upperTracts[tractIndex], catalogMotif, finalMotifList, startOffsetOf(tractIndex, finalMotifs),
            tract.startsAtRepeatEdge, tract.endsAtRepeatEdge, referenceEndPartial, sequenceSubstrings);

        int allele = 0;
        if (isHeterozygousWithDistinctAlleles)
        {
            allele = assignToAllele(
                tract, genotype->shortAlleleSizeInUnits(), genotype->longAlleleSizeInUnits(), motifLength,
                typicalReadLength);
            anyReadAssigned = anyReadAssigned || allele != 0;
        }
        GroupTallies* alleleTallies = allele == 1 ? &allele1Tallies : (allele == 2 ? &allele2Tallies : nullptr);
        if (spdlog::should_log(spdlog::level::debug))
        {
            static const char* const kReadKindNames[] = { "spanning", "flanking", "inside", "inrepeat" };
            string encoding;
            for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
            {
                encoding
                    += (sequenceSubstring.type == SequenceSubstringType::kKnownMotif
                            ? tract.bases.substr(sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length)
                            : "("
                                + tract.bases.substr(
                                    sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length)
                                + ")")
                    + " ";
            }
            spdlog::debug(
                "MotifComposition {}:{}-{} read pair {} {} allele {}: {}", locus.contigIndex, repeatStart, repeatEnd,
                tract.readPairIndex, kReadKindNames[static_cast<int>(tract.kind)], allele, encoding);
        }

        int previousMotif = -1;
        for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
        {
            int motif = -1;
            if (sequenceSubstring.type == SequenceSubstringType::kKnownMotif)
            {
                motif = finalMotifIndex.at(
                    upperTracts[tractIndex].substr(
                        sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length));
                if (!hasHighQualityBasesAt(
                        tract.bases, sequenceSubstring.offsetWithinRepeatTract, distinguishingPositions[motif]))
                {
                    motif = -1;
                }
            }
            if (motif == -1)
            {
                previousMotif = -1;
                continue;
            }
            for (GroupTallies* tallies : { &locusTallies, alleleTallies })
            {
                if (tallies == nullptr)
                {
                    continue;
                }
                tallies->motifs[motif].add(tractIndex);
                if (previousMotif != -1)
                {
                    tallies->motifPairs[{ previousMotif, motif }].add(tractIndex);
                }
            }
            previousMotif = motif;
        }
    }

    // --- Motif IDs: numbered from 1 by locus-wide occurrences, highest first (ties: motif sequence). Motifs never
    // counted get no ID.
    vector<int> countedMotifs;
    for (size_t index = 0; index != finalMotifs.size(); ++index)
    {
        if (locusTallies.motifs[index].occurrences > 0)
        {
            countedMotifs.push_back(index);
        }
    }
    std::sort(
        countedMotifs.begin(), countedMotifs.end(),
        [&](int left, int right)
        {
            const int leftCount = locusTallies.motifs[left].occurrences;
            const int rightCount = locusTallies.motifs[right].occurrences;
            return leftCount != rightCount ? leftCount > rightCount : finalMotifs[left] < finalMotifs[right];
        });

    // Which motifs the reference has depends on the known-motif list the splitter uses (it decides where to shift),
    // so split the reference again with the final list: a motif found there is a reference motif, not a new one.
    {
        splitTract(
            locus.referenceRepeatSequence, locus.referenceRepeatSequence, catalogMotif, finalMotifList, frameOffset,
            true, true, boost::none, sequenceSubstrings);
        for (const SequenceSubstring& sequenceSubstring : sequenceSubstrings)
        {
            if (sequenceSubstring.type != SequenceSubstringType::kGap)
            {
                referenceMotifs.insert(locus.referenceRepeatSequence.substr(
                    sequenceSubstring.offsetWithinRepeatTract, sequenceSubstring.length));
            }
        }
    }

    MotifComposition composition;
    vector<int> motifIds(finalMotifs.size(), 0);
    bool hasNonReferenceMotif = false;
    for (size_t rank = 0; rank != countedMotifs.size(); ++rank)
    {
        const string& motif = finalMotifs[countedMotifs[rank]];
        motifIds[countedMotifs[rank]] = rank + 1;
        composition.motifs.push_back(motif);
        hasNonReferenceMotif = hasNonReferenceMotif || referenceMotifs.count(motif) == 0;
    }
    if (onlyLociWithNonRefMotifs && !hasNonReferenceMotif)
    {
        return boost::none;
    }

    composition.locus = encodeCounts(locusTallies, motifIds);
    if (isHeterozygousWithDistinctAlleles && anyReadAssigned)
    {
        composition.hasAlleleBlocks = true;
        composition.allele1 = encodeCounts(allele1Tallies, motifIds);
        composition.allele2 = encodeCounts(allele2Tallies, motifIds);
    }
    return composition;
}

}
