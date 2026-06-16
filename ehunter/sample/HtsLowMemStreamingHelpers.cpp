#include "sample/HtsLowMemStreamingHelpers.hh"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "core/GenomicRegion.hh"
#include "core/HtsHelpers.hh"
#include "core/LocusStats.hh"
#include "core/Read.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/ParameterLoading.hh"
#include "locus/AlleleQualityMetrics.hh"
#include "reviewer/Metrics.hh"
#include "sample/MateExtractor.hh"
#include "spdlog/spdlog.h"


using std::string;

namespace ehunter
{


bool isRepeatGenotypeHomRef(const RepeatGenotype& genotype, int referenceSizeInUnits)
{
    return genotype.shortAlleleSizeInUnits() == referenceSizeInUnits
        && genotype.longAlleleSizeInUnits() == referenceSizeInUnits;
}


bool shouldFilterLocus(
    const LocusSpecification& locusSpec, const LocusFindings& locusFindings, bool skipHomRef, bool skipMissing)
{
    if (!skipHomRef && !skipMissing)
    {
        return false;
    }
    // A locus with no findings (decode failure, etc.) is a no-call we cannot classify; keep it.
    if (locusFindings.findingsForEachVariant.empty())
    {
        return false;
    }
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        const VariantFindings* findings = variantIdAndFindings.second.get();

        bool hasGenotype;
        bool isHomRef;
        const RepeatFindings* repeatFindings = dynamic_cast<const RepeatFindings*>(findings);
        if (repeatFindings != nullptr)
        {
            hasGenotype = static_cast<bool>(repeatFindings->optionalGenotype());
            if (hasGenotype)
            {
                const auto& variantSpec = locusSpec.getVariantSpecById(variantId);
                const auto repeatNodeId = variantSpec.nodes().front();
                const auto& repeatUnit = locusSpec.regionGraph().nodeSeq(repeatNodeId);
                if (repeatUnit.empty())
                {
                    throw std::runtime_error(
                        "Repeat unit sequence is empty for locus " + locusSpec.locusId() + " variant " + variantId);
                }
                const int referenceSizeInUnits = variantSpec.referenceLocus().length() / repeatUnit.length();
                isHomRef = isRepeatGenotypeHomRef(*repeatFindings->optionalGenotype(), referenceSizeInUnits);
            }
            else
            {
                isHomRef = false;
            }
        }
        else
        {
            const SmallVariantFindings* smallVariantFindings = dynamic_cast<const SmallVariantFindings*>(findings);
            if (smallVariantFindings == nullptr)
            {
                // Unclassifiable findings type: keep the locus.
                return false;
            }
            hasGenotype = static_cast<bool>(smallVariantFindings->optionalGenotype());
            isHomRef = hasGenotype && smallVariantFindings->optionalGenotype()->isHomRef();
        }

        const bool isMissing = !hasGenotype;
        const bool skippable = (skipHomRef && isHomRef) || (skipMissing && isMissing);
        if (!skippable)
        {
            return false;
        }
    }
    return true;
}


bool isAToTheLeftOfB(int32_t contigIdA, int64_t posA, int32_t contigIdB, int64_t posB, bool allowEquals) {
    if (contigIdA < contigIdB) {
        return true;
    } else if (contigIdA == contigIdB) {
         if (posA < posB) {
            return true;
         } else if (allowEquals && posA == posB) {
            return true;
         }
    }
    return false;
}

bool doesAcontainB(int32_t contigIdA, int64_t startA, int64_t endA, int32_t contigIdB, int64_t startB, int64_t endB) {
    if (contigIdA != contigIdB) {
        return false;
    }
    if (startA <= startB && endA >= endB) {
        return true;
    }
    return false;
}

bool doesAoverlapB(int32_t contigIdA, int64_t startA, int64_t endA, int32_t contigIdB, int64_t startB, int64_t endB) {
    if (contigIdA != contigIdB) {
        return false;
    }
    if (startA <= endB && endA >= startB) {
        return true;
    }
    return false;
}

bool intervalsOverlap(int start1, int end1, int start2, int end2) {
    return end1 >= start2 && start1 <= end2;
}

string minimalMotifUnderShift(const string& unit)
{
    string minimal_unit = unit;
    const string double_unit = unit + unit;
    for (size_t index = 0; index != unit.length(); ++index)
    {
        string current_unit = double_unit.substr(index, unit.length());
        if (current_unit < minimal_unit)
            minimal_unit = current_unit;
    }
    return minimal_unit;
}

string computeCanonicalMotif(const string& unit, int includeReverseComplement = false)
{
    const string minimal_unit = minimalMotifUnderShift(unit);

    if (includeReverseComplement) {
        const string unit_rc = graphtools::reverseComplement(unit);
        const string minimal_unit_rc = minimalMotifUnderShift(unit_rc);

        if (minimal_unit_rc < minimal_unit)
            return minimal_unit_rc;
    }

    return minimal_unit;
}

// Function to count the longest stretch of uninterrupted repeats in a sequence starting either from the
// left or from the right
int extendRepeatIntoSequence(std::string motif, std::string sequence, bool from_end = false) {
    // If searching from the end, reverse both the sequence and repeat unit
    if (from_end) {
        std::reverse(sequence.begin(), sequence.end());
        std::reverse(motif.begin(), motif.end());
    }

    const int repeat_length = motif.length();
    const int sequence_length = sequence.length();

    int num_pure_repeats = 0;
    for (int i = 0; i <= sequence_length - repeat_length; i += repeat_length) {
        // compare() matches a substring in place; substr() would heap-allocate a string per repeat unit.
        if (sequence.compare(i, repeat_length, motif) != 0) {
            break;
        }
        num_pure_repeats++;
    }

    return num_pure_repeats;
}

std::pair<int, bool> extendRepeatsIntoFlank(
    const std::string& locus_motif,
    const std::string& flanking_sequence,
    const std::string& soft_clip_sequence,
    bool from_end = true,
    bool is_spanning_read = false
) {
    const int locus_motif_size = locus_motif.size();
    const int flanking_sequence_size = flanking_sequence.size();
    const int soft_clip_sequence_size = soft_clip_sequence.size();
    int num_extra_repeats = 0;
    int num_flanking_repeats = 0;

    // The locus motif's canonical form is invariant for the whole call; compute it once instead of
    // recomputing it for both the flank and soft-clip comparisons below.
    const std::string canonical_locus_motif = computeCanonicalMotif(locus_motif, false);

    if (flanking_sequence_size >= locus_motif_size) {
        std::string flanking_motif = from_end
            ? flanking_sequence.substr(flanking_sequence_size - locus_motif_size)
            : flanking_sequence.substr(0, locus_motif_size);

        if (computeCanonicalMotif(flanking_motif, false) == canonical_locus_motif) {
            num_flanking_repeats = extendRepeatIntoSequence(flanking_motif, flanking_sequence, from_end);
        }
    }

    if (flanking_sequence_size < locus_motif_size ||
        flanking_sequence_size - num_flanking_repeats * locus_motif_size < locus_motif_size) {
        is_spanning_read = false;
    }

    if (soft_clip_sequence_size >= locus_motif_size &&
        flanking_sequence_size - num_flanking_repeats * locus_motif_size <= 2 * locus_motif_size + 1) {
        std::string soft_clip_motif = from_end
            ? soft_clip_sequence.substr(soft_clip_sequence_size - locus_motif_size)
            : soft_clip_sequence.substr(0, locus_motif_size);

        if (computeCanonicalMotif(soft_clip_motif, false) == canonical_locus_motif) {
            num_extra_repeats += extendRepeatIntoSequence(soft_clip_motif, soft_clip_sequence, from_end);
        }
    }

    return {num_extra_repeats, is_spanning_read};
}


// processRead consumes 0-based half-open locus coordinates (start inclusive, end exclusive),
// matching the convention used by GenomicRegion and elsewhere in the codebase.
FastReadAnalysisResult processRead(
    const FullRead& read,
    int64_t locus_start_0based,
    int64_t locus_end_0based,
    const std::string& locus_motif
) {
    FastReadAnalysisResult result;

    unsigned int read_start_0based = read.s.pos;

    unsigned int aligned_read_length = 0;
    for (const auto cigar_op : read.s.cigar) {
        int op = cigar_op & BAM_CIGAR_MASK;
        int length = cigar_op >> BAM_CIGAR_SHIFT;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP) {
            aligned_read_length += length;
        }
    }

    unsigned int read_end_0based = read_start_0based + aligned_read_length;

    if (read_start_0based >= locus_end_0based || read_end_0based <= locus_start_0based) {
        return result;  // Read does not overlap the repeat locus
    }
    result.overlaps_repeats = true;

    int locus_motif_size = locus_motif.size();
    bool spans_repeats_on_left = false;
    bool spans_repeats_on_right = false;

    int length_of_left_soft_clips = 0;  // None
    int length_of_right_soft_clips = 0;  // None

    int ref_offset = 0;
    int read_seq_position_0based = 0;

    int read_seq_position_0based_locus_start = -1; //None
    int read_seq_position_0based_locus_end = -1;  //None

    for (const auto& cigar_op : read.s.cigar) {
        int op = cigar_op & BAM_CIGAR_MASK;  // Extract operation type
        int op_length = cigar_op >> BAM_CIGAR_SHIFT;  // Extract length

        int op_ref_start_0based, op_ref_end_0based;  // half-open: [start, end)

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
            op_ref_start_0based = read_start_0based + ref_offset;
            op_ref_end_0based = op_ref_start_0based + op_length;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            op_ref_start_0based = op_ref_end_0based = read_start_0based + ref_offset;
        } else {
            if (op != BAM_CREF_SKIP && op != BAM_CHARD_CLIP) {
                spdlog::warn("Unexpected CIGAR operation '{}' in read {}", op, read.r.readId().toString());
            }
            continue;
        }

        const int min_overhang_for_spanning_read = 2;
        if (locus_start_0based - op_ref_start_0based >= min_overhang_for_spanning_read) {
            spans_repeats_on_left = true;
        }
        if (op_ref_end_0based - locus_end_0based >= min_overhang_for_spanning_read) {
            spans_repeats_on_right = true;
        }

        if (read_seq_position_0based_locus_start == -1 && op_ref_start_0based <= locus_start_0based && op_ref_end_0based >= locus_start_0based) {
            read_seq_position_0based_locus_start = read_seq_position_0based + locus_start_0based - op_ref_start_0based;
        }

        if (read_seq_position_0based_locus_end == -1 && op_ref_start_0based <= locus_end_0based && op_ref_end_0based >= locus_end_0based) {
            read_seq_position_0based_locus_end = read_seq_position_0based + locus_end_0based - op_ref_start_0based;
        }

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (op_ref_end_0based > locus_start_0based && op_ref_start_0based < locus_end_0based) {
                int op_first_base_within_locus_0based;
                if (op_ref_start_0based < locus_start_0based) {
                    op_first_base_within_locus_0based = read_seq_position_0based + (locus_start_0based - op_ref_start_0based);
                } else {
                    op_first_base_within_locus_0based = read_seq_position_0based;
                }

                int op_end_within_locus_0based;  // exclusive end (in read seq coordinates)
                if (op_ref_start_0based < locus_start_0based) {
                    op_end_within_locus_0based = std::min<int64_t>(
                        read_seq_position_0based + op_length,
                        op_first_base_within_locus_0based + (locus_end_0based - locus_start_0based)
                    );
                } else {
                    op_end_within_locus_0based = std::min<int64_t>(
                        read_seq_position_0based + op_length,
                        op_first_base_within_locus_0based + (locus_end_0based - op_ref_start_0based)
                    );
                }
                result.repeat_sequence_size_in_base_pairs += op_end_within_locus_0based - op_first_base_within_locus_0based;
            }
        } else if (op == BAM_CINS) {
            // intervalsOverlap uses closed-interval semantics; convert the half-open locus to closed bounds.
            // Padding is asymmetric: left bound by (motif_size + 1), right bound by motif_size — this matches
            // the half-open-to-closed conversion pinned by the InsertionJustBeyondPaddingIsNotCounted test.
            bool has_overlap = intervalsOverlap(op_ref_start_0based, op_ref_end_0based,
                locus_start_0based - locus_motif_size - 1, locus_end_0based + locus_motif_size);
            if (has_overlap && op_length % locus_motif_size == 0) {
                result.repeat_sequence_size_in_base_pairs += op_length;
            }
            // Inserted bases within the repeat (quality metric): count impurity insertions only —
            // those inside the locus whose length is NOT a whole number of motifs. Whole-motif
            // insertions are clean repeat-count changes (the allele-vs-reference size), already
            // captured by the genotype, so excluding them keeps this an interruption/impurity signal.
            if (op_ref_start_0based >= locus_start_0based && op_ref_start_0based <= locus_end_0based
                && op_length % locus_motif_size != 0) {
                result.inserted_bases_within_repeat += op_length;
            }
        } else if (op == BAM_CSOFT_CLIP) {
            if (ref_offset == 0) {
                length_of_left_soft_clips = op_length;
            } else {
                length_of_right_soft_clips = op_length;
            }
        } else if (op == BAM_CDEL) {
            // Deleted bases within the repeat (quality metric): impurity deletions only — the deletion's
            // overlap with the locus, counted only when the deletion length is NOT a whole number of
            // motifs. Whole-motif deletions are clean repeat-count changes already captured by the call.
            if (op_length % locus_motif_size != 0) {
                const int64_t overlap_start = std::max<int64_t>(op_ref_start_0based, locus_start_0based);
                const int64_t overlap_end = std::min<int64_t>(op_ref_end_0based, locus_end_0based);
                if (overlap_end > overlap_start) {
                    result.deleted_bases_within_repeat += static_cast<int>(overlap_end - overlap_start);
                }
            }
        }

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
            ref_offset += op_length;
        }

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            read_seq_position_0based += op_length;
        }
    }

    result.is_spanning_read = spans_repeats_on_left && spans_repeats_on_right;

    const std::string& readSequence = read.r.sequence();
    const int read_sequence_length = readSequence.length();

    std::string left_flank_bases = (read_seq_position_0based_locus_start > 0) ?
        readSequence.substr(length_of_left_soft_clips, read_seq_position_0based_locus_start - length_of_left_soft_clips)
        : "";

    std::string right_flank_bases = (read_seq_position_0based_locus_end >= 0 && read_seq_position_0based_locus_end <= read_sequence_length) ?
        (length_of_right_soft_clips > 0 ?
            readSequence.substr(read_seq_position_0based_locus_end, std::max(0, read_sequence_length - read_seq_position_0based_locus_end - length_of_right_soft_clips))
            : readSequence.substr(read_seq_position_0based_locus_end))
        : "";

    std::string left_soft_clip_sequence = length_of_left_soft_clips > 0 ? readSequence.substr(0, length_of_left_soft_clips) : "";
    std::string right_soft_clip_sequence = length_of_right_soft_clips > 0 ? readSequence.substr(read_sequence_length - length_of_right_soft_clips) : "";

    std::pair<int, bool> resultsFromExtendingToTheLeft = extendRepeatsIntoFlank(locus_motif, left_flank_bases,
        left_soft_clip_sequence, true, result.is_spanning_read);
    auto num_extra_repeats_left = resultsFromExtendingToTheLeft.first;
    auto is_spanning_read_left = resultsFromExtendingToTheLeft.second;

    std::pair<int, bool> resultsFromExtendingToTheRight = extendRepeatsIntoFlank(locus_motif, right_flank_bases,
        right_soft_clip_sequence, false, result.is_spanning_read);
    auto num_extra_repeats_right = resultsFromExtendingToTheRight.first;
    auto is_spanning_read_right = resultsFromExtendingToTheRight.second;

    result.is_spanning_read = is_spanning_read_left && is_spanning_read_right;
    if (num_extra_repeats_left > 0 || num_extra_repeats_right > 0) {
        result.soft_clipped_bases_contain_repetitive_sequence = true;
        result.repeat_sequence_size_in_base_pairs += (num_extra_repeats_left + num_extra_repeats_right) * locus_motif_size;
    }

    return result;
}


bool processLocusFast(
    const ProgramParameters& params, Reference& reference, LocusDescription& locusDescription,
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, IterativeJsonWriter& jsonWriter,
    IterativeVcfWriter& vcfWriter) {

    if (locusDescription.referenceRegions().size() != 1) {
        // fast processing is only implemented for loci with a single repeat region
        return false;
    }

    LocusSpecification locusSpec = decodeLocusSpecification(locusDescription, reference, params.heuristics(), false);
    if (locusSpec.variantSpecs().size() != 1) {
        return false;
    }
    if (locusSpec.variantSpecs().front().classification().type != VariantType::kRepeat) {
        // fast processing only handles repeat variants; SmallVariant loci must go through full genotyping
        return false;
    }
    const auto& locusReferenceRegion = locusDescription.referenceRegions().front();
    const auto repeatNodeId = locusSpec.variantSpecs().front().nodes().front();
    const string& locusMotif = locusSpec.regionGraph().nodeSeq(repeatNodeId);

    LocusStatsCalculatorFromReadAlignments locusStatsCalculator(locusDescription.chromType(), locusReferenceRegion);

	int totalReads = 0;
    int mappedReadCount = 0;
    float averageMapQAtLocus = 0.0;
    for (const auto& readPair : readPairs) {
        // numMatesSet() is 2 for a normal pair and 1 for a single-ended entry (mate unmapped/unavailable).
        totalReads += readPair->numMatesSet();
        if (readPair->firstMate && readPair->firstMate->s.isMapped) {
            averageMapQAtLocus += readPair->firstMate->s.mapq;
            ++mappedReadCount;
        }
        if (readPair->secondMate && readPair->secondMate->s.isMapped) {
            averageMapQAtLocus += readPair->secondMate->s.mapq;
            ++mappedReadCount;
        }
        locusStatsCalculator.inspect(*readPair);
    }
    if (mappedReadCount > 0) {
        averageMapQAtLocus /= mappedReadCount;
    }

    // iterate over each read pair
    std::vector<int> soft_clipped_read_repeat_sequence_sizes;
    std::map<int, int> allele_size_spanning_read_votes;
    std::map<int, int> allele_size_soft_clipped_read_votes;
    std::map<int, int> allele_size_flanking_read_votes;
    // Per-spanning-allele-size accumulators for the allele quality metrics (keyed by repeat size in bp).
    struct SpanningReadStats {
        int forwardReads = 0;
        int reverseReads = 0;
        long insertedBases = 0;
        long deletedBases = 0;
    };
    std::map<int, SpanningReadStats> allele_size_spanning_read_stats;
    for (const auto& readPair : readPairs) {

        // secondMate is empty for a single-ended entry (mate unmapped/unavailable); skip it below.
        const FullRead* reads[] = {
            readPair->firstMate ? &*readPair->firstMate : nullptr,
            readPair->secondMate ? &*readPair->secondMate : nullptr };

        // process the read and its alignment stats, followed by the mate and its alignment stats
        for (const auto* read : reads) {
            if (read == nullptr) {
                continue;
            }
            if (!read->s.isMapped || read->s.isSupplementaryAlignment || read->s.isSecondaryAlignment) {
                continue;
            }

            if (read->s.mapq <= 3 && averageMapQAtLocus >= 20) {
                // ignore reads with unusually low mapQ
                continue;
            }

            //print the read you're processing
            const FastReadAnalysisResult& readAnalysisResult = processRead(
               *read, locusReferenceRegion.start(), locusReferenceRegion.end(), locusMotif);

            if (!readAnalysisResult.overlaps_repeats) {
                continue;
            }


            if (readAnalysisResult.is_spanning_read && !readAnalysisResult.soft_clipped_bases_contain_repetitive_sequence) {
                const int spanningSizeBp = readAnalysisResult.repeat_sequence_size_in_base_pairs;
                allele_size_spanning_read_votes[spanningSizeBp]++;
                SpanningReadStats& stats = allele_size_spanning_read_stats[spanningSizeBp];
                if (read->r.isReversed()) {
                    stats.reverseReads++;
                } else {
                    stats.forwardReads++;
                }
                stats.insertedBases += readAnalysisResult.inserted_bases_within_repeat;
                stats.deletedBases += readAnalysisResult.deleted_bases_within_repeat;
            }
            if (readAnalysisResult.soft_clipped_bases_contain_repetitive_sequence) {
                soft_clipped_read_repeat_sequence_sizes.push_back(readAnalysisResult.repeat_sequence_size_in_base_pairs);
            }
            // Flanking read: overlaps the repeat but anchored on only one side (not spanning), and its
            // soft-clips carry no repeat sequence (those are handled separately above as larger-allele
            // support). repeat_sequence_size_in_base_pairs is a lower bound on this read's repeat content.
            if (!readAnalysisResult.is_spanning_read && !readAnalysisResult.soft_clipped_bases_contain_repetitive_sequence) {
                allele_size_flanking_read_votes[readAnalysisResult.repeat_sequence_size_in_base_pairs]++;
            }
        }
    }

    int num_soft_clipped_reads_supporting_larger_allele_than_any_spanning_read = 0;

    for (int soft_clipped_read_repeat_sequence_size : soft_clipped_read_repeat_sequence_sizes) {
        //
        std::vector<int> consistent_with_repeat_sequence_sizes;
        for (auto it = allele_size_spanning_read_votes.begin(); it != allele_size_spanning_read_votes.end(); ++it) {
            int s = it->first;
            if (soft_clipped_read_repeat_sequence_size <= s) {
                consistent_with_repeat_sequence_sizes.push_back(s);
            }
        }

        // if this soft clipped repeat sequence size is consistent with only one of the allele sizes from spanning reads,
        // record it as a vote for that allele size
        if (consistent_with_repeat_sequence_sizes.size() == 1) {
            allele_size_soft_clipped_read_votes[consistent_with_repeat_sequence_sizes.front()]++;
        } else if (consistent_with_repeat_sequence_sizes.empty()) {
            num_soft_clipped_reads_supporting_larger_allele_than_any_spanning_read++;
        }
    }

	if (num_soft_clipped_reads_supporting_larger_allele_than_any_spanning_read >= 2) {
        return false;
    }

    // combine spanning and soft-clipped read votes
    std::map<int, int> allele_size_votes;
    for (auto it = allele_size_spanning_read_votes.begin(); it != allele_size_spanning_read_votes.end(); ++it) {
        int s = it->first;
        int c = it->second;
        allele_size_votes[s] = c + allele_size_soft_clipped_read_votes[s];
    }

    // keep alleles with at least 2 votes
    std::vector<std::pair<int, int>> allele_size_votes_list;
    for (auto it = allele_size_votes.begin(); it != allele_size_votes.end(); ++it) {
        int s = it->first;
        int c = it->second;
        if (c >= 2) {
            allele_size_votes_list.emplace_back(s, c);
        }
    }

    if (allele_size_votes_list.empty()) {
        return false;
    }

    // Sort alleles first by vote count (descending), then by allele size (ascending)
    std::sort(allele_size_votes_list.begin(), allele_size_votes_list.end(),
        [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
            return (a.second > b.second) || (a.second == b.second && a.first < b.first);
        });


    const AlleleCount alleleCount = determineExpectedAlleleCount(locusDescription.chromType(), params.sample().sex());
    const int maxAlleles = alleleCount == AlleleCount::kOne ? 1 : 2;
    std::vector<std::pair<int, int>> top_genotypes;
    size_t top_genotypes_count = std::min<size_t>(maxAlleles, allele_size_votes_list.size());
    if (top_genotypes_count > 0) {
        top_genotypes.assign(allele_size_votes_list.begin(), allele_size_votes_list.begin() + top_genotypes_count);
    }

    if (top_genotypes.empty() && totalReads >= 5) {
        return false;
    }

    const int locusMotifSize = locusMotif.size();

    // STR stutter correction: PCR/sequencing slippage produces reads 1-2 repeat units off a homozygous
    // allele, which the naive top-2-by-votes can mis-call as a spurious heterozygous second allele. When
    // the second (minor-vote) allele is +/-1 or +/-2 units from the major and carries few votes, collapse
    // to a homozygous-major genotype. The vote-ratio thresholds below are motif- and gap-specific: they
    // were tuned against the HG002 truth set so that each motif size's collapses stay at least ~2.5:1
    // correct. Thresholds tighten with motif size because longer motifs host more genuine close
    // heterozygotes. Motifs >6bp are left untouched (stutter is negligible and real close hets dominate).
    // Runs after the fast path commits its allele selection, so the fast/full split and speed are unchanged.
    if (maxAlleles == 2 && top_genotypes.size() == 2 && locusMotifSize >= 1 && locusMotifSize <= 6) {
        const int cMajor = top_genotypes[0].second;
        const int cMinor = top_genotypes[1].second;
        const int gapUnits = (std::abs(top_genotypes[0].first - top_genotypes[1].first) + locusMotifSize / 2) / locusMotifSize;
        double stutterVoteRatioThreshold = 0.0;
        if (gapUnits == 1) {
            switch (locusMotifSize) {
                case 1: case 2: stutterVoteRatioThreshold = 0.30; break;
                case 3:         stutterVoteRatioThreshold = 0.28; break;
                case 4: case 5: stutterVoteRatioThreshold = 0.25; break;
                case 6:         stutterVoteRatioThreshold = 0.21; break;
            }
        } else if (gapUnits == 2) {
            switch (locusMotifSize) {
                case 1: case 2: stutterVoteRatioThreshold = 0.25; break;
                case 3:         stutterVoteRatioThreshold = 0.21; break;
                case 4:         stutterVoteRatioThreshold = 0.19; break;
                case 5:         stutterVoteRatioThreshold = 0.18; break;
                case 6:         stutterVoteRatioThreshold = 0.12; break;
            }
        }
        if (stutterVoteRatioThreshold > 0.0 && cMinor < stutterVoteRatioThreshold * cMajor) {
            top_genotypes.resize(1); // collapse stutter -> homozygous major allele
        }
    }

    std::vector<int> num_repeats;
    for (const auto& g : top_genotypes) {
        num_repeats.push_back(g.first / locusMotifSize);
    }
    if (num_repeats.size() > 1) {
   		//sort alleles by size (ascending)
		std::sort(num_repeats.begin(), num_repeats.end());
   	} else if (num_repeats.size() == 1 && alleleCount == AlleleCount::kTwo) {
		num_repeats.push_back(num_repeats.front()); // make it a diploid homozygous genotype
	}

    const LocusStats& locusStats = locusStatsCalculator.estimate(params.sample().sex());
    LocusFindings locusFindings(locusStats);

    const auto& variantId = locusDescription.variantIds().front(); // Access the first and only variant ID

	CountTable countsOfSpanningReads;
	CountTable countsOfFlankingReads;
	CountTable countsOfInrepeatReads;
	// Initialize spanning reads count table
	for (const auto& alleleVote : allele_size_spanning_read_votes) {
		countsOfSpanningReads.setCountOf(alleleVote.first / locusMotifSize, alleleVote.second);
	}
	// Initialize flanking reads count table (one-flank-anchored reads; see read loop above)
	for (const auto& alleleVote : allele_size_flanking_read_votes) {
		countsOfFlankingReads.setCountOf(alleleVote.first / locusMotifSize, alleleVote.second);
	}

	boost::optional<RepeatGenotype> repeatGenotype;
	if (!num_repeats.empty()) {
		repeatGenotype = RepeatGenotype(locusMotifSize, num_repeats);
	}

	// Derive a per-allele genotype confidence interval from the spread of high-quality reads. Both
	// allele_size_spanning_read_votes and soft_clipped_read_repeat_sequence_sizes are populated only
	// from reads that passed the fast-path quality filter in the read loop above (mapped, primary,
	// mapQ-screened), so the CI is built exclusively from high-quality reads. Each supported spanning
	// size (>=1 read) is assigned to its nearest called allele and widens that allele's CI to span
	// [min, max] of its cluster; a soft-clipped read longer than the long allele (a lower bound on a
	// possible expansion) extends the long allele's upper bound.
	if (repeatGenotype) {
		const int shortAlleleUnits = repeatGenotype->shortAlleleSizeInUnits();
		const int longAlleleUnits = repeatGenotype->longAlleleSizeInUnits();
		int shortLo = shortAlleleUnits, shortHi = shortAlleleUnits;
		int longLo = longAlleleUnits, longHi = longAlleleUnits;
		for (const auto& alleleVote : allele_size_spanning_read_votes) {
			if (alleleVote.second <= 0) {
				continue;
			}
			const int units = alleleVote.first / locusMotifSize;
			if (std::abs(units - shortAlleleUnits) <= std::abs(units - longAlleleUnits)) {
				shortLo = std::min(shortLo, units);
				shortHi = std::max(shortHi, units);
			} else {
				longLo = std::min(longLo, units);
				longHi = std::max(longHi, units);
			}
		}
		// Homozygous call: both alleles share the single supported cluster.
		if (shortAlleleUnits == longAlleleUnits) {
			longLo = shortLo;
			longHi = shortHi;
		}
		// Soft-clipped reads longer than the long allele extend its upper bound (possible under-call).
		for (const int softClipSizeBp : soft_clipped_read_repeat_sequence_sizes) {
			longHi = std::max(longHi, softClipSizeBp / locusMotifSize);
		}
		repeatGenotype->setShortAlleleSizeInUnitsCi(shortLo, shortHi);
		repeatGenotype->setLongAlleleSizeInUnitsCi(longLo, longHi);
	}

    // Apply --skip-hom-ref before building findings, so skipped loci don't pay the allocation. Fast
    // genotyping always produces a called genotype here (missing genotypes fall back to full genotyping),
    // so --skip-missing-genotypes never applies on this path.
    if (params.skipHomRef() && repeatGenotype) {
        const int referenceSizeInUnits = locusReferenceRegion.length() / locusMotifSize;
        if (isRepeatGenotypeHomRef(*repeatGenotype, referenceSizeInUnits)) {
            return true;  // genotyped hom-ref: emit no record
        }
    }

	// Per-allele quality metrics from the high-quality spanning reads (the same mapQ-filtered reads as
	// the CI above). Computed in fast mode: depth (DP), high-quality-unambiguous read count, strand-bias
	// phred, mean inserted/deleted bases within the repeat, and CI-width / allele-size. NOT computed here
	// (left at 0): qd (needs per-base qualities or a graph-alignment score) and the flank-normalized
	// depths (need a flank-window depth tally). These approximate the full genotyper, which derives them
	// from graph realignment, so they will not match it exactly.
	RepeatAlleleQualityMetrics qualityMetrics;
	CountTable countsOfHighQualityUnambiguousReads;
	if (repeatGenotype) {
		const int shortAlleleUnits = repeatGenotype->shortAlleleSizeInUnits();
		const int longAlleleUnits = repeatGenotype->longAlleleSizeInUnits();
		const bool isHet = repeatGenotype->numAlleles() == 2 && shortAlleleUnits != longAlleleUnits;
		// clusterTarget: -1 = all reads (hom/hemi single allele), 0 = reads nearer the short allele,
		// 1 = reads nearer the long allele. Each spanning read is assigned to its nearest called allele,
		// matching the CI clustering above.
		auto buildAlleleMetrics = [&](int alleleNumber, int alleleUnits, const NumericInterval& ci,
			int clusterTarget) {
			AlleleMetrics allele;
			allele.alleleNumber = alleleNumber;
			allele.alleleSize = alleleUnits;
			int depth = 0;
			int forwardReads = 0;
			int reverseReads = 0;
			long insertedBases = 0;
			long deletedBases = 0;
			for (const auto& sizeAndVotes : allele_size_spanning_read_votes) {
				if (clusterTarget != -1) {
					const int units = sizeAndVotes.first / locusMotifSize;
					const bool nearerShort =
						std::abs(units - shortAlleleUnits) <= std::abs(units - longAlleleUnits);
					if ((nearerShort ? 0 : 1) != clusterTarget) {
						continue;
					}
				}
				depth += sizeAndVotes.second;
				const auto statsIt = allele_size_spanning_read_stats.find(sizeAndVotes.first);
				if (statsIt != allele_size_spanning_read_stats.end()) {
					forwardReads += statsIt->second.forwardReads;
					reverseReads += statsIt->second.reverseReads;
					insertedBases += statsIt->second.insertedBases;
					deletedBases += statsIt->second.deletedBases;
				}
			}
			allele.depth = depth;
			// Spanning reads are mapQ-filtered (high quality) and each pins exactly one size
			// (unambiguous), so they are this allele's high-quality-unambiguous reads.
			allele.highQualityUnambiguousReads = depth;
			allele.strandBiasBinomialPhred = reviewer::computeStrandBiasBinomialPhred(forwardReads, reverseReads);
			if (depth > 0) {
				allele.meanInsertedBasesWithinRepeats = static_cast<double>(insertedBases) / depth;
				allele.meanDeletedBasesWithinRepeats = static_cast<double>(deletedBases) / depth;
				countsOfHighQualityUnambiguousReads.setCountOf(alleleUnits, depth);
			}
			if (alleleUnits > 0) {
				allele.confidenceIntervalDividedByAlleleSize =
					static_cast<double>(ci.end() - ci.start()) / alleleUnits;
			}
			qualityMetrics.alleles.push_back(allele);
		};
		if (isHet) {
			buildAlleleMetrics(1, shortAlleleUnits, repeatGenotype->shortAlleleSizeInUnitsCi(), 0);
			buildAlleleMetrics(2, longAlleleUnits, repeatGenotype->longAlleleSizeInUnitsCi(), 1);
		} else {
			buildAlleleMetrics(1, shortAlleleUnits, repeatGenotype->shortAlleleSizeInUnitsCi(), -1);
		}
		qualityMetrics.variantId = variantId;
		qualityMetrics.hasMetrics = !qualityMetrics.alleles.empty();
	}

	auto repeatFindingsPtr = std::make_unique<RepeatFindings>(
		countsOfSpanningReads, countsOfFlankingReads, countsOfInrepeatReads,
		locusStats.alleleCount(), repeatGenotype, GenotypeFilter());
	repeatFindingsPtr->setQuickGenotype(true);  // genotyped via the fast path
	if (qualityMetrics.hasMetrics) {
		repeatFindingsPtr->setAlleleQualityMetrics(qualityMetrics);
		repeatFindingsPtr->setCountsOfHighQualityUnambiguousReads(countsOfHighQualityUnambiguousReads);
	}
	locusFindings.findingsForEachVariant.emplace(variantId, std::move(repeatFindingsPtr));

    jsonWriter.addRecord(locusSpec, locusFindings);
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        vcfWriter.addRecord(variantId, locusSpec, locusFindings);
    }

    return true;
}

bool writeZeroCoverageRecord(
    const ProgramParameters& params, Reference& reference, const LocusDescription& locusDescription,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter)
{
    // A zero-coverage locus has all-missing genotypes, so --skip-missing-genotypes excludes it entirely
    // (no record, and no graph build).
    if (params.skipMissingGenotypes())
    {
        return false;
    }
    try
    {
        // extendFlanks=false builds the graph from the locus structure alone, skipping the
        // (~regionExtensionLength) reference flank reads. The repeat unit comes from the structure regex,
        // so RepeatUnit and all other catalog-derived fields still match seeking's output exactly. This is
        // only safe for single-region (single-repeat) loci — the same loci processLocusFast handles with
        // extendFlanks=false. The flankless graph build is unsupported for multi-region loci (interruptions
        // or multiple repeats) and corrupts the blueprint there, so those fall back to the flank-extending
        // decode (a small reference read), which matches seeking exactly.
        const bool extendFlanks = locusDescription.referenceRegions().size() != 1;
        LocusSpecification locusSpec
            = decodeLocusSpecification(locusDescription, reference, params.heuristics(), extendFlanks);

        const AlleleCount alleleCount
            = determineExpectedAlleleCount(locusDescription.chromType(), params.sample().sex());
        LocusFindings locusFindings(LocusStats(alleleCount, 0, 0, 0));

        // Empty per-variant findings, byte-for-byte equal to what the analyzers return on low depth:
        // empty read counts, no genotype (-> ./.), and the LowDepth genotype filter.
        for (const VariantSpecification& variantSpec : locusSpec.variantSpecs())
        {
            std::unique_ptr<VariantFindings> variantFindingsPtr;
            if (variantSpec.classification().type == VariantType::kRepeat)
            {
                variantFindingsPtr = std::make_unique<RepeatFindings>(
                    CountTable(), CountTable(), CountTable(), alleleCount, boost::none, GenotypeFilter::kLowDepth);
            }
            else
            {
                variantFindingsPtr = std::make_unique<SmallVariantFindings>(
                    0, 0, AlleleCheckSummary(AlleleStatus::kUncertain, 0),
                    AlleleCheckSummary(AlleleStatus::kUncertain, 0), alleleCount, boost::none,
                    GenotypeFilter::kLowDepth);
            }
            locusFindings.findingsForEachVariant.emplace(variantSpec.id(), std::move(variantFindingsPtr));
        }

        jsonWriter.addRecord(locusSpec, locusFindings);
        vcfWriter.addRecords(locusSpec, locusFindings);
    }
    catch (const MissingContigError& e)
    {
        // Locus sits on a contig absent from the reference FASTA (e.g. a _fix patch contig present in the
        // read file header but not in the FASTA). This is benign and expected, so warn rather than error.
        spdlog::warn("Skipping zero-coverage locus {}: {}", locusDescription.locusId(), e.what());
        jsonWriter.addSkippedRecord(locusDescription.locusId(), "error");
    }
    catch (const std::exception& e)
    {
        // A malformed locus would previously have been emitted as a bare skipped record (it has no
        // coverage, so it never reached graph construction). Preserve that fallback instead of aborting
        // the whole run.
        spdlog::error("Error emitting zero-coverage record for {}: {}", locusDescription.locusId(), e.what());
        jsonWriter.addSkippedRecord(locusDescription.locusId(), "error");
    }
    return true;
}

}
