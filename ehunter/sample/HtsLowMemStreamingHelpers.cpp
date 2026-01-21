#include "sample/HtsLowMemStreamingHelpers.hh"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "core/GenomicRegion.hh"
#include "core/HtsHelpers.hh"
#include "core/LocusStats.hh"
#include "core/Read.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/ParameterLoading.hh"
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
        if (sequence.substr(i, repeat_length) != motif) {
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

    if (flanking_sequence_size >= locus_motif_size) {
        std::string flanking_motif = from_end
            ? flanking_sequence.substr(flanking_sequence_size - locus_motif_size)
            : flanking_sequence.substr(0, locus_motif_size);

        if (computeCanonicalMotif(flanking_motif, false) == computeCanonicalMotif(locus_motif, false)) {
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

        if (computeCanonicalMotif(soft_clip_motif, false) == computeCanonicalMotif(locus_motif, false)) {
            num_extra_repeats += extendRepeatIntoSequence(soft_clip_motif, soft_clip_sequence, from_end);
        }
    }

    return {num_extra_repeats, is_spanning_read};
}


FastReadAnalysisResult processRead(
    const FullRead& read,
    int64_t locus_start_1based,
    int64_t locus_end_1based,
    const std::string& locus_motif
) {
    FastReadAnalysisResult result;

    unsigned int read_start_1based = read.s.pos + 1;

    unsigned int aligned_read_length = 0;
    for (const auto cigar_op : read.s.cigar) {
        int op = cigar_op & BAM_CIGAR_MASK;
        int length = cigar_op >> BAM_CIGAR_SHIFT;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CREF_SKIP) {
            aligned_read_length += length;
        }
    }

    unsigned int read_end_1based = read.s.pos + aligned_read_length;

    if (read_start_1based > locus_end_1based || read_end_1based < locus_start_1based) {
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
    int read_seq_position_1based_locus_end = -1;  //None

    for (const auto& cigar_op : read.s.cigar) {
        int op = cigar_op & BAM_CIGAR_MASK;  // Extract operation type
        int op_length = cigar_op >> BAM_CIGAR_SHIFT;  // Extract length

        int op_ref_start_1based, op_ref_end_1based;

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CDEL) {
            op_ref_start_1based = read_start_1based + ref_offset;
            op_ref_end_1based = op_ref_start_1based + op_length;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            op_ref_start_1based = op_ref_end_1based = read_start_1based + ref_offset;
        } else {
            if (op != BAM_CREF_SKIP && op != BAM_CHARD_CLIP) {
                spdlog::warn("Unexpected CIGAR operation '{}' in read {}", op, read.r.readId().toString());
            }
            continue;
        }

        const int min_overhang_for_spanning_read = 2;
        if (locus_start_1based - op_ref_start_1based >= min_overhang_for_spanning_read) {
            spans_repeats_on_left = true;
        }
        if (op_ref_end_1based - locus_end_1based >= min_overhang_for_spanning_read) {
            spans_repeats_on_right = true;
        }

        if (read_seq_position_0based_locus_start == -1 && op_ref_start_1based <= locus_start_1based && op_ref_end_1based >= locus_start_1based) {
            read_seq_position_0based_locus_start = read_seq_position_0based + locus_start_1based - op_ref_start_1based;
        }

        if (read_seq_position_1based_locus_end == -1 && op_ref_start_1based <= locus_end_1based && op_ref_end_1based >= locus_end_1based) {
            read_seq_position_1based_locus_end = read_seq_position_0based + locus_end_1based - op_ref_start_1based + 1;
        }

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (op_ref_end_1based >= locus_start_1based && op_ref_start_1based <= locus_end_1based) {
                int op_first_base_within_locus_0based;
                if (op_ref_start_1based < locus_start_1based) {
                    op_first_base_within_locus_0based = read_seq_position_0based + (locus_start_1based - op_ref_start_1based);
                } else {
                    op_first_base_within_locus_0based = read_seq_position_0based;
                }

                int op_last_base_within_locus_1based;
                if (op_ref_start_1based < locus_start_1based) {
                    op_last_base_within_locus_1based = std::min<int64_t>(
                        read_seq_position_0based + op_length,
                        op_first_base_within_locus_0based + (locus_end_1based - locus_start_1based + 1)
                    );
                } else {
                    op_last_base_within_locus_1based = std::min<int64_t>(
                        read_seq_position_0based + op_length,
                        op_first_base_within_locus_0based + (locus_end_1based - op_ref_start_1based + 1)
                    );
                }
                result.repeat_sequence_size_in_base_pairs += op_last_base_within_locus_1based - op_first_base_within_locus_0based;
            }
        } else if (op == BAM_CINS) {
            bool has_overlap = intervalsOverlap(op_ref_start_1based, op_ref_end_1based,
                locus_start_1based - locus_motif_size - 1, locus_end_1based + locus_motif_size + 1);
            if (has_overlap && op_length % locus_motif_size == 0) {
                result.repeat_sequence_size_in_base_pairs += op_length;
            }
        } else if (op == BAM_CSOFT_CLIP) {
            if (ref_offset == 0) {
                length_of_left_soft_clips = op_length;
            } else {
                length_of_right_soft_clips = op_length;
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

    std::string right_flank_bases = (read_seq_position_1based_locus_end >= 0 && read_seq_position_1based_locus_end <= read_sequence_length) ?
        (length_of_right_soft_clips > 0 ?
            readSequence.substr(read_seq_position_1based_locus_end, read_sequence_length - read_seq_position_1based_locus_end - length_of_right_soft_clips)
            : readSequence.substr(read_seq_position_1based_locus_end))
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
    const auto& locusReferenceRegion = locusDescription.referenceRegions().front();
    const auto repeatNodeId = locusSpec.variantSpecs().front().nodes().front();
    const string& locusMotif = locusSpec.regionGraph().nodeSeq(repeatNodeId);

    LocusStatsCalculatorFromReadAlignments locusStatsCalculator(locusDescription.chromType(), locusReferenceRegion);

	int totalReads = 0;
    float averageMapQAtLocus = 0.0;
    for (const auto& readPair : readPairs) {
        totalReads += 2;
        if (readPair->firstMate->s.isMapped) {
            averageMapQAtLocus += readPair->firstMate->s.mapq;
        }
        if (readPair->secondMate->s.isMapped) {
            averageMapQAtLocus += readPair->secondMate->s.mapq;
        }
        locusStatsCalculator.inspect(*readPair);
    }
    if (readPairs.size() > 0) {
        averageMapQAtLocus /= totalReads;
    }

    // iterate over each read pair
    std::vector<int> soft_clipped_read_repeat_sequence_sizes;
    std::map<int, int> allele_size_spanning_read_votes;
    std::map<int, int> allele_size_soft_clipped_read_votes;
    for (const auto& readPair : readPairs) {

        std::vector<FullRead> reads = { *readPair->firstMate, *readPair->secondMate };

        // process the read and its alignment stats, followed by the mate and its alignment stats
        for (const auto& read : reads) {
            if (!read.s.isMapped || read.s.isSupplementaryAlignment || read.s.isSecondaryAlignment) {
                continue;
            }

            if (read.s.mapq <= 3 && averageMapQAtLocus >= 20) {
                // ignore reads with unusually low mapQ
                continue;
            }

            //print the read you're processing
            const FastReadAnalysisResult& readAnalysisResult = processRead(
               read, locusReferenceRegion.start(), locusReferenceRegion.end(), locusMotif);

            if (!readAnalysisResult.overlaps_repeats) {
                continue;
            }


            if (readAnalysisResult.is_spanning_read && !readAnalysisResult.soft_clipped_bases_contain_repetitive_sequence) {
                allele_size_spanning_read_votes[readAnalysisResult.repeat_sequence_size_in_base_pairs]++;
            }
            if (readAnalysisResult.soft_clipped_bases_contain_repetitive_sequence) {
                soft_clipped_read_repeat_sequence_sizes.push_back(readAnalysisResult.repeat_sequence_size_in_base_pairs);
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

	boost::optional<RepeatGenotype> repeatGenotype;
	if (!num_repeats.empty()) {
		repeatGenotype = RepeatGenotype(locusMotifSize, num_repeats);
	}

    // Check if genotype is hom-ref and skip output if --skip-hom-ref is enabled
    if (params.skipHomRef() && repeatGenotype) {
        const int referenceSizeInUnits = locusReferenceRegion.length() / locusMotifSize;
        if (isRepeatGenotypeHomRef(*repeatGenotype, referenceSizeInUnits)) {
            return true;  // Skip output for hom-ref loci
        }
    }

	std::unique_ptr<VariantFindings> variantFindingsPtr = std::make_unique<RepeatFindings>(
		countsOfSpanningReads, countsOfFlankingReads, countsOfInrepeatReads,
		locusStats.alleleCount(), repeatGenotype, GenotypeFilter());
	locusFindings.findingsForEachVariant.emplace(variantId, std::move(variantFindingsPtr));


    jsonWriter.addRecord(locusSpec, locusFindings);
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        vcfWriter.addRecord(variantId, locusSpec, locusFindings);
    }

    return true;
}

}
