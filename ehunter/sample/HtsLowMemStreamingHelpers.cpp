#include "sample/HtsLowMemStreamingHelpers.hh"

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "core/GenomicRegion.hh"
#include "core/HtsHelpers.hh"
#include "core/LocusStats.hh"
#include "core/Read.hh"
#include "core/RepeatPurity.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/LocusSpecDecoding.hh"
#include "io/ParameterLoading.hh"
#include "locus/AlleleQualityMetrics.hh"
#include "reviewer/ConsensusSequence.hh"
#include "reviewer/Metrics.hh"
#include "reviewer/ReviewerWorkflow.hh"
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

    // Read-sequence offsets of the two locus edges, tracked independently of the anchors above for the
    // consensus tract (repeat_tract_*). Three kinds of CIGAR operation can place an edge, each handled
    // in its own branch of the loop below:
    //   - an aligned match (M/=/X), which puts a real read base at the reference position;
    //   - a whole-motif insertion abutting the edge, whose bases the genotype vote already counts;
    //   - a deletion containing the edge, where the read simply has no base at those positions and the
    //     current read position is already the first base past (or last base before) the deletion.
    // A soft clip never places an edge -- it is not aligned, so it carries no reference position. Both
    // offsets stay -1 when nothing places the corresponding edge, which excludes the read. For a spanning
    // read that cannot actually happen: M/=/X and D ops tile the aligned reference span with no gaps, so
    // each edge always falls in one of them. The residual exclusion is a zero-length tract, i.e. an
    // alignment that deletes the whole locus.
    // The start edge uses half-open containment (locus_start < op_ref_end, not <=), so an operation that
    // merely ENDS at the locus start does not claim it -- the next operation, the one actually covering
    // the first repeat base, does. The end edge keeps the inclusive bound for matches, since the tract is
    // half-open on the right and so ends exactly at an operation's end.
    // Claiming is first-wins via the == -1 guards, with one deliberate exception: the whole-motif
    // insertion branch EXTENDS an end edge the preceding match already placed (it is guarded by
    // tract_end == read_seq_position instead, so it can only ever extend the edge it directly abuts).
    int read_seq_position_0based_tract_start = -1;
    int read_seq_position_0based_tract_end = -1;

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

        // Consensus-tract edges. One branch per operation that can place an edge -- match, whole-motif
        // insertion abutting the edge, deletion containing the edge; see the declarations above for the
        // reasoning behind each and for why the start bound is strict.
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (read_seq_position_0based_tract_start == -1 && op_ref_start_0based <= locus_start_0based
                && locus_start_0based < op_ref_end_0based) {
                read_seq_position_0based_tract_start
                    = read_seq_position_0based + locus_start_0based - op_ref_start_0based;
            }

            if (read_seq_position_0based_tract_end == -1 && op_ref_start_0based <= locus_end_0based
                && locus_end_0based <= op_ref_end_0based) {
                read_seq_position_0based_tract_end
                    = read_seq_position_0based + locus_end_0based - op_ref_start_0based;
            }
        } else if (op == BAM_CINS && op_length % locus_motif_size == 0) {
            // A whole-motif insertion sitting exactly on a locus edge is how an aligner represents an
            // expansion whose extra copies it pushed just outside the repeat region. The genotype vote
            // already counts those bases (the BAM_CINS branch below), so the tract has to take them too:
            // otherwise this read reports the REFERENCE-length sequence for an expanded allele, and its
            // tract stops matching the allele it voted for. Insertions further out are separated from
            // the tract by flank bases and cannot be spliced in, so only the abutting ones qualify.
            if (read_seq_position_0based_tract_start == -1 && op_ref_start_0based == locus_start_0based) {
                read_seq_position_0based_tract_start = read_seq_position_0based;
            }

            if (op_ref_start_0based == locus_end_0based
                && read_seq_position_0based_tract_end == read_seq_position_0based) {
                read_seq_position_0based_tract_end += op_length;
            }
        } else if (op == BAM_CDEL) {
            // A locus edge landing inside a deletion still has a well-defined tract boundary: the read
            // carries no base at those reference positions, and a deletion consumes no read bases, so the
            // read position here already IS the first base after (or the last base before) the deleted
            // stretch. The genotype vote treats deletions the same way -- a D adds nothing to
            // repeat_sequence_size_in_base_pairs -- so the tract and the vote stay in agreement.
            // This case is common rather than exotic: aligners left-align an STR deletion onto the locus
            // boundary, so refusing these reads drops the shorter allele of many heterozygous calls.
            if (read_seq_position_0based_tract_start == -1 && op_ref_start_0based <= locus_start_0based
                && locus_start_0based < op_ref_end_0based) {
                read_seq_position_0based_tract_start = read_seq_position_0based;
            }

            if (read_seq_position_0based_tract_end == -1 && op_ref_start_0based <= locus_end_0based
                && locus_end_0based < op_ref_end_0based) {
                read_seq_position_0based_tract_end = read_seq_position_0based;
            }
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

    // ReadRepeatPurity (fast path): compare the read's in-repeat tract (read-sequence coordinates
    // [locus_start, locus_end)) against a perfect consecutive repeat sequence of the motif anchored
    // at the repeat-region start.
    // Only spanning reads feed the metric downstream, but the counts are filled whenever both bounds
    // are known so the function stays self-contained and testable.
    if (read_seq_position_0based_locus_start >= 0
        && read_seq_position_0based_locus_end >= read_seq_position_0based_locus_start
        && read_seq_position_0based_locus_end <= read_sequence_length) {
        const std::string repeat_tract = readSequence.substr(
            read_seq_position_0based_locus_start,
            read_seq_position_0based_locus_end - read_seq_position_0based_locus_start);
        const RepeatSequencePurity purity = computeRepeatSequencePurity(repeat_tract, locus_motif);
        result.matched_bases_within_repeat = static_cast<int>(purity.matchedBases);
        result.repeat_read_bases = static_cast<int>(purity.totalBases);
    }

    // The read's repeat tract, exported for per-allele consensus building: the read bases between the
    // two locus edges. Deliberately computed from its own anchors rather than from the purity tract
    // above, so ReadRepeatPurity is unchanged by this.
    if (read_seq_position_0based_tract_start >= 0
        && read_seq_position_0based_tract_end > read_seq_position_0based_tract_start
        && read_seq_position_0based_tract_end <= read_sequence_length) {
        result.repeat_tract_start = read_seq_position_0based_tract_start;
        result.repeat_tract_length
            = read_seq_position_0based_tract_end - read_seq_position_0based_tract_start;
    }

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
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, bool reservoirSampled,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter, CapturedLocusOutput* captured) {

    // Per-locus fast-path timing (thread-CPU clock), mirroring the full-genotyper path in
    // HtsLowMemStreamingSampleAnalysis.cpp; covers graph decode + read processing + heuristic genotyping.
    // Only sampled under --output-genotype-timing. Started before the structural fall-back checks below; the
    // measured value is only read on the success path that actually emits a fast-path genotype.
    const bool recordGenotypingTime = params.outputGenotypeTiming();
    timespec genotypingStart{};
    if (recordGenotypingTime) { clock_gettime(CLOCK_THREAD_CPUTIME_ID, &genotypingStart); }

    // --plot-all requests a REViewer image for every locus. The fast path performs no read-to-graph
    // realignment, so it can neither render an image nor compute the graph-based metrics REViewer needs;
    // decline up front so the full genotyper handles (and plots) this locus. (--plot-all and
    // --disable-all-plots are mutually exclusive, so the guard is belt-and-suspenders.)
    // Skip this when --quick-heuristic-genotyping-only is set: full genotyping is disabled in that mode,
    // so declining would only drop the locus to a skipped record (no genotype, no plot) -- keep the fast
    // genotype instead.
    if (params.plotAll() && !params.disableAllPlots() && !params.heuristicGenotypingOnly()) {
        return false;
    }

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
    const int locusMotifSize = locusMotif.size();

    // Coverage/ReadLength/FragmentLength are measured over the reference flanks, the same way the full
    // genotyper measures them from graph alignments, so the emitted fields mean the same thing in both
    // modes. The flanks are regionExtensionLength bases on either side of the repeat -- exactly the
    // window this locus's reads were collected over, so every read the estimate needs is already cached.
    LocusStatsCalculatorFromReadAlignments locusStatsCalculator(
        locusDescription.chromType(), locusReferenceRegion, params.heuristics().regionExtensionLength(),
        reference.contigInfo().getContigSize(locusReferenceRegion.contigIndex()));

    int mappedReadCount = 0;
    float averageMapQAtLocus = 0.0;
    for (const auto& readPair : readPairs) {
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
        long matchedRepeatReadBases = 0; // ReadRepeatPurity numerator (bases matching the motif)
        long totalRepeatReadBases = 0;   // ReadRepeatPurity denominator (repeat-region read bases)
    };
    std::map<int, SpanningReadStats> allele_size_spanning_read_stats;

    // Per-allele consensus sequences (on unless --dont-output-consensus-sequences). Spanning-read repeat
    // tracts are collected in two nested levels, because the two quantities involved are NOT
    // interchangeable and each level needs a different one:
    //   - the OUTER key is the allele size in whole motif units the read voted for
    //     (repeat_sequence_size_in_base_pairs / motif size, exactly how num_repeats below derives the
    //     called allele). Keying on this guarantees every read ends up under the allele it actually
    //     supported, so a called allele can never come up empty while reads voted for it;
    //   - the INNER key is the TRACT LENGTH -- the number of read bases between the two locus edges.
    //     Grouping by it makes every tract in a group the same length, which is what lets the consensus
    //     be a plain column-wise majority with no realignment.
    // The two diverge when the locus reference length is not a whole multiple of the motif, and per-read
    // when the aligner places an indel at or near a locus boundary (a whole-motif insertion just outside
    // the locus counts toward the vote via the padding in the BAM_CINS branch of processRead, but
    // contributes no tract bases). Keying the outer level on the vote is what keeps those reads with
    // their own allele; only the tract length they contribute differs, which the inner level absorbs.
    const bool buildConsensus = params.enableConsensusSequences();
    // `bases` points into each read's own sequence buffer rather than copying: the locus cache owning the
    // reads outlives this call, so the bases stay valid for the whole function. Kept as one flat vector
    // so the read loop pays a single push_back per read; the grouping by allele and tract length happens
    // once afterwards, off the hot path (and, for --output-genotype-timing, outside the genotyping clock).
    struct RepeatTract { const char* bases; int length; int alleleUnits; };
    std::vector<RepeatTract> spanning_read_tracts;

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
                stats.matchedRepeatReadBases += readAnalysisResult.matched_bases_within_repeat;
                stats.totalRepeatReadBases += readAnalysisResult.repeat_read_bases;
                if (buildConsensus && readAnalysisResult.repeat_tract_length > 0) {
                    spanning_read_tracts.push_back(
                        RepeatTract{ read->r.sequence().data() + readAnalysisResult.repeat_tract_start,
                                     readAnalysisResult.repeat_tract_length,
                                     spanningSizeBp / locusMotifSize });
                }
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
	if (repeatGenotype && params.enableAlleleQualityMetrics()) {
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
			long matchedRepeatBases = 0;
			long purityMatchedBases = 0;  // ReadRepeatPurity numerator (repeat-region read bases matching the motif)
			long purityTotalBases = 0;    // ReadRepeatPurity denominator (repeat-region read bases)
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
				// Each spanning read traversed a repeat tract of sizeAndVotes.first bp, so that is its matched-
				// base count within the repeat (the fast-path analogue of the full genotyper's STR-node matches).
				matchedRepeatBases += static_cast<long>(sizeAndVotes.first) * sizeAndVotes.second;
				const auto statsIt = allele_size_spanning_read_stats.find(sizeAndVotes.first);
				if (statsIt != allele_size_spanning_read_stats.end()) {
					forwardReads += statsIt->second.forwardReads;
					reverseReads += statsIt->second.reverseReads;
					insertedBases += statsIt->second.insertedBases;
					deletedBases += statsIt->second.deletedBases;
					purityMatchedBases += statsIt->second.matchedRepeatReadBases;
					purityTotalBases += statsIt->second.totalRepeatReadBases;
				}
			}
			// DP is a base-level coverage estimate (matched repeat bases / allele length in bp), matching
			// the full genotyper's allele depth so the field carries the same units in both modes. For reads
			// that exactly span the called allele this equals the spanning-read count; it diverges when reads
			// observe a different size than the called allele.
			const long alleleLengthBp = static_cast<long>(alleleUnits) * locusMotifSize;
			allele.depth = alleleLengthBp > 0 ? static_cast<double>(matchedRepeatBases) / alleleLengthBp : 0.0;
			// Spanning reads are mapQ-filtered (high quality) and each pins exactly one size (unambiguous),
			// so the COUNT of them is this allele's high-quality-unambiguous read count.
			allele.highQualityUnambiguousReads = depth;
			allele.strandBiasBinomialPhred = reviewer::computeStrandBiasBinomialPhred(forwardReads, reverseReads);
			if (depth > 0) {
				allele.meanInsertedBasesWithinRepeats = static_cast<double>(insertedBases) / depth;
				allele.meanDeletedBasesWithinRepeats = static_cast<double>(deletedBases) / depth;
				countsOfHighQualityUnambiguousReads.setCountOf(alleleUnits, depth);
			}
			// Denominator is AlleleSize+1 (Laplace pseudo-count): always defined even at
			// alleleUnits==0, and matches the calibrator's ci_over_eh = ci_width/(eh+1).
			allele.confidenceIntervalDividedByAlleleSize =
				static_cast<double>(ci.end() - ci.start()) / (alleleUnits + 1);
			// ReadRepeatPurity from this allele's spanning reads; -1 when no repeat-region read bases.
			allele.readRepeatPurity = (purityTotalBases > 0)
				? static_cast<double>(purityMatchedBases) / purityTotalBases
				: -1.0;
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
	repeatFindingsPtr->setReservoirSampling(reservoirSampled);  // read set capped by --max-depth

	// Stop the genotyping clock at the same point it stopped before consensus output existed, and keep
	// the consensus build below outside it: that build is output annotation, not genotyping, and is on
	// by default, so timing it would silently redefine the field.
	// Two caveats on GenotypingTimeMillis, both small but worth stating rather than implying otherwise:
	//   - processRead's tract-edge branches and its tract export are NOT gated on consensus output, so
	//     they run for every read in every configuration and sit inside this interval. The field is
	//     therefore slightly higher than in releases predating this feature, with or without
	//     --dont-output-consensus-sequences. Gating them would mean threading the flag through
	//     processRead to save a handful of integer comparisons per CIGAR operation.
	//   - the one flag-gated part inside the interval is the per-read push_back into
	//     spanning_read_tracts, so the number does also shift a little with the flag.
	if (recordGenotypingTime) {
		timespec genotypingEnd{};
		clock_gettime(CLOCK_THREAD_CPUTIME_ID, &genotypingEnd);
		repeatFindingsPtr->setGenotypingTimeMillis(
			(genotypingEnd.tv_sec - genotypingStart.tv_sec) * 1e3
			+ (genotypingEnd.tv_nsec - genotypingStart.tv_nsec) / 1e6);
	}
	if (qualityMetrics.hasMetrics) {
		repeatFindingsPtr->setAlleleQualityMetrics(qualityMetrics);
		repeatFindingsPtr->setCountsOfHighQualityUnambiguousReads(countsOfHighQualityUnambiguousReads);
	}

	// Per-allele consensus sequences: a column-wise majority vote over the spanning reads that voted for
	// this allele. Reuses the full genotyper's PositionVotes / AlleleConsensus so both paths emit
	// identically formatted strings ('N' for uncovered positions, the same read-support digit string).
	// The CONTENT will not match the full path, which votes with graph-aligned anchor fragments (flanking
	// reads included) rather than exactly-spanning reads only -- the same approximation already noted for
	// the allele quality metrics above.
	std::vector<std::string> consensusSequences;
	std::vector<std::string> consensusReadSupport;
	if (buildConsensus && repeatGenotype && repeatGenotype->longAlleleSizeInUnits() > 0) {
		// The full genotyper hardcodes this same weight for every base, so with a constant weight both
		// paths reduce to an unweighted majority vote.
		constexpr double kConsensusQualityWeight = 0.99;
		auto buildAlleleConsensus = [&](int alleleIndex, int alleleUnits) {
			// This allele's reads are exactly those that voted for its unit count. Among them only
			// equal-length tracts can be stacked column-wise, so take the best-supported tract length;
			// std::map iterates in ascending key order and the comparison is strict, so ties resolve to
			// the shortest tract and the choice is deterministic.
			std::map<int, int> readsPerTractLength;
			for (const auto& tract : spanning_read_tracts) {
				if (tract.alleleUnits == alleleUnits) {
					readsPerTractLength[tract.length]++;
				}
			}
			int bestTractLength = 0;
			int bestTractReads = 0;
			for (const auto& lengthAndCount : readsPerTractLength) {
				if (lengthAndCount.second > bestTractReads) {
					bestTractReads = lengthAndCount.second;
					bestTractLength = lengthAndCount.first;
				}
			}

			reviewer::AlleleConsensus consensus;
			consensus.alleleIndex = alleleIndex;
			// The number of bases these reads actually carry between the locus edges. Equals
			// alleleUnits * locusMotifSize on a clean locus, but differs when the locus length is not a
			// whole number of motifs or the aligner put an indel at an edge. With no supporting tract,
			// fall back to the called allele length so the all-'N' string still has a meaningful size.
			consensus.repeatLength = (bestTractReads > 0) ? bestTractLength : alleleUnits * locusMotifSize;
			consensus.positions.resize(consensus.repeatLength);
			consensus.anchorReadCount = bestTractReads;
			for (const auto& tract : spanning_read_tracts) {
				// bestTractLength stays 0 when this allele has no tract at all, and a stored tract is
				// never zero-length, so this correctly votes nothing and leaves an all-'N' consensus.
				if (tract.alleleUnits != alleleUnits || tract.length != bestTractLength) {
					continue;
				}
				for (int position = 0; position < tract.length; ++position) {
					if (!consensus.positions[position]) {
						consensus.positions[position] = reviewer::PositionVotes();
					}
					consensus.positions[position]->addVote(tract.bases[position], kConsensusQualityWeight);
				}
			}
			consensusSequences.push_back(consensus.toString());
			consensusReadSupport.push_back(consensus.toReadSupportString());
		};
		const int shortAlleleUnits = repeatGenotype->shortAlleleSizeInUnits();
		const int longAlleleUnits = repeatGenotype->longAlleleSizeInUnits();
		buildAlleleConsensus(0, shortAlleleUnits);
		// Homozygous calls get a single consensus, matching the full genotyper.
		if (repeatGenotype->numAlleles() == 2 && shortAlleleUnits != longAlleleUnits) {
			buildAlleleConsensus(1, longAlleleUnits);
		}
	}

	if (!consensusSequences.empty()) {
		repeatFindingsPtr->setConsensusSequences(consensusSequences);
		repeatFindingsPtr->setConsensusReadSupport(consensusReadSupport);
	}
	locusFindings.findingsForEachVariant.emplace(variantId, std::move(repeatFindingsPtr));

	// Conditional plotting: the fast path has now produced the genotype + metrics, so evaluate the
	// catalog plot conditions against these findings. If any fires, this locus needs a REViewer image,
	// which the fast path can't render (no graph realignment) -- decline so the full genotyper re-runs
	// and plots it. (Respects --disable-all-plots; the --plot-all case was handled at function entry.)
	// As above, skip when --quick-heuristic-genotyping-only is set: declining would only drop the locus
	// to a skipped record instead of running full genotyping, so keep the fast genotype.
	if (!params.disableAllPlots() && !params.heuristicGenotypingOnly()
		&& reviewer::shouldPlotReadVisualization(locusSpec, locusFindings)) {
		return false;
	}

    jsonWriter.addRecord(locusSpec, locusFindings, captured ? &captured->json : nullptr);
    for (const auto& variantIdAndFindings : locusFindings.findingsForEachVariant)
    {
        const string& variantId = variantIdAndFindings.first;
        vcfWriter.addRecord(variantId, locusSpec, locusFindings, captured ? &captured->vcf : nullptr);
    }

    return true;
}

bool writeZeroCoverageRecord(
    const ProgramParameters& params, Reference& reference, const LocusDescription& locusDescription,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter, CapturedLocusOutput* captured)
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

        jsonWriter.addRecord(locusSpec, locusFindings, captured ? &captured->json : nullptr);
        vcfWriter.addRecords(locusSpec, locusFindings, captured ? &captured->vcf : nullptr);
    }
    catch (const MissingContigError& e)
    {
        // Locus sits on a contig absent from the reference FASTA (e.g. a _fix patch contig present in the
        // read file header but not in the FASTA). This is benign and expected, so warn rather than error.
        spdlog::warn("Skipping zero-coverage locus {}: {}", locusDescription.locusId(), e.what());
        jsonWriter.addSkippedRecord(locusDescription.locusId(), "error", captured ? &captured->json : nullptr);
    }
    catch (const std::exception& e)
    {
        // A malformed locus would previously have been emitted as a bare skipped record (it has no
        // coverage, so it never reached graph construction). Preserve that fallback instead of aborting
        // the whole run.
        spdlog::error("Error emitting zero-coverage record for {}: {}", locusDescription.locusId(), e.what());
        jsonWriter.addSkippedRecord(locusDescription.locusId(), "error", captured ? &captured->json : nullptr);
    }
    return true;
}

}
