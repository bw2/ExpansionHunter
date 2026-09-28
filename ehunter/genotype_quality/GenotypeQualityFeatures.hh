//
// ExpansionHunter
//
// Per-allele feature assembly for the genotype-quality model. This is the C++
// port of genotype_quality/eh_json_features.extract_variant_rows + features.py:
// it turns the typed EH findings for one called allele into the float feature
// vector (in features.py order) that the model evaluates, and routes the allele
// to its genotyping genotypingRegime. Kept as pure functions over lightweight inputs
// (CountTable + AlleleMetrics + primitives) so it is unit-testable without the
// full analysis machinery. See GENOTYPE_QUALITY_CPP_INTEGRATION_PLAN.md.
//

#pragma once

#include <limits>
#include <string>
#include <vector>

#include "core/CountTable.hh"
#include "locus/AlleleQualityMetrics.hh"

namespace ehunter
{
namespace gq
{

// One model-set per genotyping genotypingRegime (mirrors features.REGIMES).
enum class GenotypingRegime
{
    Quick,
    FullSpanning,
    FullNonspanning
};

// Locus-level inputs shared by every allele of a variant.
struct LocusFeatureContext
{
    int motifSize = 0;
    int refSizeBp = 0;
    const CountTable& spanningReads;
    const CountTable& flankingReads;
    const CountTable& hqUnambiguousReads;
    // The fraction of bases in the reference repeat sequence that matches a perfect repeat sequence. 
    // -1.0 = not computed, mapped to NaN in the feature vector (the model was trained on NaN for missing).
    double referenceRepeatPurity = -1.0;
    // Read depth over the locus. Must be the SAME value JsonWriter emits as LocusResults.Coverage
    // (rounded to 2 decimals), since that rounded number is what the model was trained on. NaN = not
    // supplied; production always supplies it, so NaN only arises in tests (it makes a tree take its
    // learned missing-value branch rather than treating some sentinel as a real depth).
    double coverage = std::numeric_limits<double>::quiet_NaN();
    // Shape of the genotype this allele came from: total alleles called (1 = hemizygous, 2 = diploid)
    // and how many of them are distinct sizes (1 = hom, 2 = het). 0 = not set, mapped to NaN the same
    // way; also test-only, since these are derived from a genotype the caller already has.
    int numAlleles = 0;
    int numDistinctAlleles = 0;
};

// Number of DISTINCT called sizes in a genotype -- the model's `n_distinct_alleles` feature
// (1 = hom or hemizygous, 2 = het). Mirrors the training side's `len(set(genotype))` in
// eh_json.extract_variant_rows. `isHomozygous` is a front-vs-back comparison, so a 1-allele
// genotype correctly yields 1.
int numDistinctAllelesOf(bool isHomozygous, int numAlleles);

// Returns the genotypingRegime for an allele (mirrors features.genotyping_regime_of): the quick
// path is always Quick; otherwise >=1 spanning read at the called size is
// FullSpanning, else FullNonspanning.
GenotypingRegime genotypingRegimeOf(bool quickGenotype, int spanningAtCalled);

// The canonical feature names this assembler produces, in the same order as the
// values returned by assembleFeatures for the given genotypingRegime (features.py QUICK_FEATURES
// for Quick, FULL_FEATURES for the full genotyping_regimes). The annotator zips these with the
// assembled values to resolve a model's declared feature_names by name.
const std::vector<std::string>& featureNamesForGenotypingRegime(GenotypingRegime genotypingRegime);

// Builds the model feature vector for one allele, in features.py order
// (QUICK_FEATURES = 27 entries for Quick, FULL_FEATURES = 29 for the full genotyping_regimes).
// `eh` is the called allele size in repeat units; `ciStart`/`ciEnd` its confidence
// interval; `aqm` the matching per-allele quality metrics.
std::vector<double> assembleFeatures(
    const LocusFeatureContext& ctx, int alleleRank, int eh, int ciStart, int ciEnd,
    const AlleleMetrics& aqm, GenotypingRegime genotypingRegime);

} // namespace gq
} // namespace ehunter
