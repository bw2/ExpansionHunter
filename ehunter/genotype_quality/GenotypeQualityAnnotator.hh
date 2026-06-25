//
// ExpansionHunter
//
// Glue that composes the genotype-quality feature assembler and model evaluator
// into a single per-allele prediction. The JSON writer calls predictAllele for
// each called allele and emits PredictedLengthCorrectionFactor / pTooShort / pTooLong.
// Kept separate so this composition is unit-testable without the writer.
//

#pragma once

#include "genotype_quality/GenotypeQualityFeatures.hh"
#include "genotype_quality/GenotypeQualityModel.hh"

namespace ehunter
{
namespace gq
{

// Routes one allele to its genotypingRegime, assembles its feature vector, and runs both
// model heads. `quickGenotype` is the fast-path flag; `eh` the called allele size
// in repeat units; `ciStart`/`ciEnd` its confidence interval; `aqm` its quality
// metrics.
AllelePrediction predictAllele(
    const GenotypeQualityModel& model, bool quickGenotype, const LocusFeatureContext& ctx,
    int alleleRank, int eh, int ciStart, int ciEnd, const AlleleMetrics& aqm);

} // namespace gq
} // namespace ehunter
