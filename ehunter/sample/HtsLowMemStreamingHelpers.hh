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

bool processLocusFast(
	const ProgramParameters& params, Reference& reference, LocusDescription& locusDescription,
    const std::vector<std::shared_ptr<FullReadPair>>& readPairs, IterativeJsonWriter& jsonWriter,
    IterativeVcfWriter& vcfWriter);
}
