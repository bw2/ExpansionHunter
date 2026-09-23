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

// Motif composition of a repeat locus (--output-motif-composition): counts of each motif and of each pair of
// adjacent motifs in the reads' repeat sequence, for the whole locus and, when the two alleles differ enough in
// length, for each allele. Counting uses the reads' original BAM/CRAM alignments rather than their graph
// alignments, so reads made of a motif other than the catalog motif are kept.

#pragma once

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include <boost/optional.hpp>

#include "genotyping/RepeatGenotype.hh"

namespace ehunter
{

struct FullReadPair;

// Counts for one group of reads (the whole locus, or the reads assigned to one allele). Each value is
// (occurrences, reads): how many times the motif or the pair was counted, and how many distinct reads hold at
// least one counted occurrence of it.
struct MotifCompositionCounts
{
    std::map<int, std::pair<int, int>> motifs; // motif ID -> (occurrences, reads)
    std::map<std::pair<int, int>, std::pair<int, int>> motifPairs; // (motif ID, next motif ID) -> (occurrences, reads)
};

struct MotifComposition
{
    // motifs[i] is the motif with ID i + 1. IDs are shared by every group at the locus and are numbered by
    // locus-wide occurrences, highest first, with ties broken by motif sequence.
    std::vector<std::string> motifs;
    MotifCompositionCounts locus;
    // Set only for heterozygous calls whose alleles differ by at least 2 repeat units, when at least one read
    // could be assigned to an allele. allele1 is the shorter allele, as in the Genotype field.
    bool hasAlleleBlocks = false;
    MotifCompositionCounts allele1;
    MotifCompositionCounts allele2;
};

// Everything the calculation needs to know about one repeat variant. All coordinates are 0-based half-open.
struct MotifCompositionLocus
{
    int32_t contigIndex = -1;
    int64_t referenceRepeatStart = 0; // S
    int64_t referenceRepeatEnd = 0; // E
    std::string catalogMotif; // catalog motif; may contain IUPAC codes (for example AARRG)
    std::string referenceRepeatSequence; // reference sequence of [S, E), upper case
    // Motifs the catalog lists as known at this locus ("KnownMotifs"), as returned by selectCatalogKnownMotifs. They
    // are matched by rotation: since the reads are cut in line with the catalog motif, they show the rotation that
    // actually occurs. For each listed motif, the candidate new motif among its rotations seen most often is trusted
    // without the error test new motifs need, once it has their minimum read-pair support; its other rotations are
    // frame shifts and are not counted. A trusted motif counts as a
    // non-reference motif unless the reference repeat has it. Listed motifs do not otherwise steer the splitter, and
    // other new motifs are still discovered from the reads.
    std::vector<std::string> catalogKnownMotifs;
    // Up to 30 reference bases just before S and just after E, upper case. Used to find where a read's repeat
    // sequence runs into the flank; empty if unknown.
    std::string leftFlankSequence;
    std::string rightFlankSequence;
    // Mean fragment length at the locus. A mate counts as anchored only if it starts within this distance of the
    // repeat. 0 means unknown, in which case computeMotifComposition's flankLength is used.
    int meanFragmentLength = 0;
};

// Motif lengths the calculation supports: at least 2, and at most a third of the read length.
bool isEligibleForMotifComposition(int motifLength, int typicalReadLength);

// The catalog's KnownMotifs for a locus, in upper case and in the catalog's order, without the entries the
// calculation cannot use: motifs whose length differs from the catalog motif's, motifs with bases other than A, C,
// G and T, and repeats of an earlier entry, including rotations of it (AAC after CAA). Each dropped entry is logged
// as a warning that names the locus.
std::vector<std::string> selectCatalogKnownMotifs(
    const std::vector<std::string>& knownMotifs, int motifLength, const std::string& locusId);

// Computes the motif composition of one repeat variant from the reads EH holds for its locus. The reads are
// only read, never modified. When onlyLociWithNonRefMotifs is true, returns boost::none unless some counted motif
// does not occur in the reference repeat (and is not the catalog motif itself).
//
// typicalReadLength is the genome-wide read length probed at startup. flankLength is how far from the repeat the
// locus's reads were collected (--region-extension-length): a read BWA placed confidently within this distance of
// the repeat, but not in it, is never an in-repeat read.
boost::optional<MotifComposition> computeMotifComposition(
    const MotifCompositionLocus& locus, int typicalReadLength, int flankLength,
    const std::vector<const FullReadPair*>& readPairs, const boost::optional<RepeatGenotype>& genotype,
    bool onlyLociWithNonRefMotifs);

// The building blocks below are exposed for unit testing.
namespace motifcomposition
{

// Fraction of positions whose base equals the base `lag` positions earlier (case-insensitive; N never matches).
double computePeriodScore(const std::string& sequence, int lag);

// True if the sequence repeats with period motifLength: its period score at lag motifLength is at least
// minScore, and beats the score at every proper divisor of motifLength (including 1) by at least 0.1, so that
// homopolymers and repeats of a shorter period are rejected.
bool passesPeriodTest(const std::string& sequence, int motifLength, double minScore);

// The offset in [0, motif length) at which tiling the motif over the sequence gives the fewest mismatches
// (ties: smallest offset). Every offset is scored over the same number of whole motif-sized windows.
int computeFrameOffset(const std::string& sequence, const std::string& motif);

// Number of soft-clipped bases to keep as repeat sequence, walking away from the aligned part. `sequence` is
// the read, the clip is [clipStart, clipEnd), and the aligned repeat bases next to it are [tractStart, clipStart)
// for a clip on the right (clipOnRight) or [clipEnd, tractEnd) for a clip on the left. Bases are kept while the
// fraction of bases equal to the base one motif length back (toward the aligned part) stays at or above
// minScore over a window of max(2 * motifLength, 12) comparisons.
int computeSoftClipBasesToKeep(
    const std::string& sequence, int tractStart, int clipStart, int clipEnd, int tractEnd, bool clipOnRight,
    int motifLength, double minScore);

enum class SequenceSubstringType
{
    kKnownMotif, // a motif from the list of known motifs
    kCatalogMotifMatch, // matches the catalog motif via IUPAC codes, but is not on the list of known motifs
    kNewMotifCandidate, // matches neither; a possible new motif
    kGap // bases that do not form a motif-sized substring (an indel, a partial unit, or a rejected candidate)
};

struct SequenceSubstring
{
    int offsetWithinRepeatTract;
    int length;
    SequenceSubstringType type;
};

// Splits a tract into motif-sized substrings and gaps. knownMotifs are concrete upper-case motifs
// ordered most common first; the shift search compares against the catalog motif and the first 8 of them.
// startOffset is where the first substring starts. Candidate new motifs that do not touch another substring or a
// trusted end of the tract on both sides are turned into gaps.
//
// Whether an end is trusted: startsAtRepeatEdge / endsAtRepeatEdge mean the read's alignment placed that end of the
// tract exactly at the repeat's edge in the reference. referenceEndPartial is the reference repeat's partial last
// unit (possibly empty); at an end placed at the repeat's edge the tract must end with a partial unit of the same
// length that resembles it (an empty one means the last substring must end exactly at the tract end). Pass
// boost::none when the tract is the reference repeat itself, whose end needs no check.
std::vector<SequenceSubstring> splitIntoMotifs(
    const std::string& tract, const std::string& catalogMotif, const std::vector<std::string>& knownMotifs,
    int startOffset, bool startsAtRepeatEdge, bool endsAtRepeatEdge,
    const boost::optional<std::string>& referenceEndPartial);

// P(X >= count) for X ~ Poisson(mean).
double computePoissonUpperTail(double mean, int count);

}

}
