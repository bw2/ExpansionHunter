#pragma once

#include "core/GenomicRegion.hh"
#include "core/Parameters.hh"
#include "core/Read.hh"
#include "genotyping/RepeatGenotype.hh"
#include "io/IterativeJsonWriter.hh"
#include "io/IterativeVcfWriter.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"


namespace ehunter
{

struct FastReadAnalysisResult {
    int repeat_sequence_size_in_base_pairs = 0;
    bool overlaps_repeats = false;
    bool is_spanning_read = false;
    bool soft_clipped_bases_contain_repetitive_sequence = false;
};


// Check if region A is to the left of region B
bool isAToTheLeftOfB(int32_t contigIdA, int64_t posA, int32_t contigIdB, int64_t posB, bool allowEquals = false);

// Check if region A fully contains region B
bool doesAcontainB(int32_t contigIdA, int64_t startA, int64_t endA, int32_t contigIdB, int64_t startB, int64_t endB);

// Check if regions A and B overlap
bool doesAoverlapB(int32_t contigIdA, int64_t startA, int64_t endA, int32_t contigIdB, int64_t startB, int64_t endB);


// Check if a RepeatGenotype is homozygous reference
bool isRepeatGenotypeHomRef(const RepeatGenotype& genotype, int referenceSizeInUnits);

// Decide whether a locus should be excluded from the VCF/JSON output given the --skip-hom-ref and
// --skip-missing-genotypes flags. The locus is filtered out only if it has at least one variant and EVERY
// variant is "skippable": homozygous reference (when skipHomRef) or missing/no-genotype (when skipMissing).
// A locus with any non-reference call, or any unclassifiable variant, is always kept. With skipMissing=false
// this reproduces the original --skip-hom-ref behavior (skip iff all variants are called homozygous ref);
// with both flags set, loci whose variants are any mixture of hom-ref and missing are also filtered out.
bool shouldFilterLocus(
    const LocusSpecification& locusSpec, const LocusFindings& locusFindings, bool skipHomRef, bool skipMissing);

// Analyze a single read against a repeat locus given in 0-based half-open coordinates
// (start inclusive, end exclusive). Exposed here so it can be unit tested directly.
FastReadAnalysisResult processRead(
    const FullRead& read, int64_t locus_start_0based, int64_t locus_end_0based, const std::string& locus_motif);

bool processLocusFast(
	const ProgramParameters& params, Reference& reference, LocusDescription& locusDescription,
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, IterativeJsonWriter& jsonWriter,
    IterativeVcfWriter& vcfWriter);

// Emit a no-call record for a locus with zero coverage, identical to what seeking/streaming mode produce,
// without running genotyping. For single-region loci it builds a flankless stub LocusSpecification
// (extendFlanks=false) so the reference flanks are not read and the only reference access is the single
// left-flank base the VCF REF column requires. Multi-region loci instead use the flank-extending decode
// (extendFlanks=true), which reads ~regionExtensionLength reference flank bases per side; see the
// implementation comment in HtsLowMemStreamingHelpers.cpp for why. It then writes empty findings whose
// per-variant values mirror what RepeatAnalyzer/SmallVariantAnalyzer return on low depth (empty read
// counts, no genotype, LowDepth filter). The stub spec is local and freed as soon as the record is
// written, preserving the streaming mode's low-memory profile.
// Returns true if a no-call record was written, false if the locus was suppressed (--skip-missing-genotypes).
bool writeZeroCoverageRecord(
    const ProgramParameters& params, Reference& reference, const LocusDescription& locusDescription,
    IterativeJsonWriter& jsonWriter, IterativeVcfWriter& vcfWriter);
}
