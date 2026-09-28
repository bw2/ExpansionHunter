//
// ExpansionHunter
//
// Unit tests for the genotype-quality per-allele feature assembler.
//

#include "genotype_quality/GenotypeQualityFeatures.hh"

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "gmock/gmock.h"

using namespace ehunter;
using namespace ehunter::gq;

namespace
{

AlleleMetrics makeAqm()
{
    AlleleMetrics aqm;
    aqm.alleleNumber = 1;
    aqm.alleleSize = 20;
    aqm.depth = 30.5;
    aqm.qd = 4.0;
    aqm.meanInsertedBasesWithinRepeats = 0.3;
    aqm.meanDeletedBasesWithinRepeats = 0.1;
    aqm.readRepeatPurity = 0.97;
    aqm.strandBiasBinomialPhred = 1.2;
    aqm.leftFlankNormalizedDepth = 0.9;
    aqm.rightFlankNormalizedDepth = 1.1;
    aqm.highQualityUnambiguousReads = 12;
    aqm.confidenceIntervalDividedByAlleleSize = 0.285;
    return aqm;
}

} // namespace

TEST(GenotypeQualityFeatures, GenotypingRegimeRouting)
{
    EXPECT_EQ(genotypingRegimeOf(true, 0), GenotypingRegime::Quick);  // quick path -> Quick regardless of spanning
    EXPECT_EQ(genotypingRegimeOf(true, 5), GenotypingRegime::Quick);
    EXPECT_EQ(genotypingRegimeOf(false, 1), GenotypingRegime::FullSpanning);
    EXPECT_EQ(genotypingRegimeOf(false, 8), GenotypingRegime::FullSpanning);
    EXPECT_EQ(genotypingRegimeOf(false, 0), GenotypingRegime::FullNonspanning);
}

TEST(GenotypeQualityFeatures, FullGenotypingRegimeVectorMatchesFeaturesPyOrder)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}, {22, 1}});
    const CountTable flanking(std::map<int32_t, int32_t>{{18, 1}, {25, 2}});
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq, /*referenceRepeatPurity=*/0.92,
                                  /*coverage=*/41.5, /*numAlleles=*/2, /*numDistinctAlleles=*/2};

    const std::vector<double> f = assembleFeatures(ctx, /*rank=*/0, /*eh=*/20, /*ciStart=*/18,
                                                   /*ciEnd=*/24, makeAqm(), GenotypingRegime::FullSpanning);

    ASSERT_EQ(f.size(), 29u);
    EXPECT_DOUBLE_EQ(f[0], (3.0));             // motif_size
    EXPECT_DOUBLE_EQ(f[1], (10.0));            // num_repeats_in_reference = 30/3
    EXPECT_DOUBLE_EQ(f[2], (30.0));            // ref_size_bp
    EXPECT_DOUBLE_EQ(f[3], (20.0));            // eh
    EXPECT_DOUBLE_EQ(f[4], (10.0));            // eh_minus_ref = 20 - 10
    EXPECT_DOUBLE_EQ(f[5], (0.0));             // allele_rank
    EXPECT_DOUBLE_EQ(f[6], (2.0));             // n_alleles
    EXPECT_DOUBLE_EQ(f[7], (2.0));             // n_distinct_alleles (het)
    EXPECT_DOUBLE_EQ(f[8], (6.0));             // ci_width = 24 - 18
    EXPECT_DOUBLE_EQ(f[9], ((4.0 - 2.0) / 7.0)); // ci_asymmetry = ((24-20)-(20-18))/(6+1)
    EXPECT_DOUBLE_EQ(f[10], (6.0 / 21.0));     // ci_over_eh = 6/(20+1)
    EXPECT_DOUBLE_EQ(f[11], (9.0));            // spanning_total = 8 + 1
    EXPECT_DOUBLE_EQ(f[12], (5.0));            // hq_unamb_total
    EXPECT_DOUBLE_EQ(f[13], (3.0));            // flanking_total = 1 + 2
    EXPECT_DOUBLE_EQ(f[14], (8.0));            // spanning_at_called = count at 20
    EXPECT_DOUBLE_EQ(f[15], (1.0));            // spanning_above_called = count of 22
    EXPECT_DOUBLE_EQ(f[16], (2.0));            // flanking_above_called = count of 25
    EXPECT_DOUBLE_EQ(f[17], (8.0 / 9.0));      // support_frac
    EXPECT_DOUBLE_EQ(f[18], (3.0 / 12.0));     // flanking_frac = 3 / (3 flanking + 9 spanning)
    EXPECT_DOUBLE_EQ(f[19], (41.5));           // coverage
    EXPECT_DOUBLE_EQ(f[20], (30.5));           // depth
    EXPECT_DOUBLE_EQ(f[21], (12.0));           // hq_unambiguous_reads
    EXPECT_DOUBLE_EQ(f[22], (1.2));            // strand_bias_phred
    EXPECT_DOUBLE_EQ(f[23], (0.3));            // mean_inserted_bases
    EXPECT_DOUBLE_EQ(f[24], (0.1));            // mean_deleted_bases
    EXPECT_DOUBLE_EQ(f[25], (0.92));           // reference_repeat_purity
    EXPECT_DOUBLE_EQ(f[26], (0.97));           // read_repeat_purity
    EXPECT_DOUBLE_EQ(f[27], (0.9));            // left_flank_norm_depth
    EXPECT_DOUBLE_EQ(f[28], (1.1));            // right_flank_norm_depth
}

TEST(GenotypeQualityFeatures, NumDistinctAllelesOfMatchesLenSetGenotype)
{
    // The (isHomozygous, numAlleles) pairs RepeatGenotype actually produces, against the training
    // side's len(set(genotype)) in eh_json.extract_variant_rows.
    EXPECT_EQ(numDistinctAllelesOf(/*isHomozygous=*/true, /*numAlleles=*/1), 1);   // "20"     -> {20}
    EXPECT_EQ(numDistinctAllelesOf(/*isHomozygous=*/true, /*numAlleles=*/2), 1);   // "20/20"  -> {20}
    EXPECT_EQ(numDistinctAllelesOf(/*isHomozygous=*/false, /*numAlleles=*/2), 2);  // "20/97"  -> {20,97}
}

TEST(GenotypeQualityFeatures, HomozygousAndHemizygousDistinctAlleleCounts)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    // Diploid hom call (20/20): two alleles, one distinct size.
    const LocusFeatureContext hom{3, 30, spanning, flanking, hq, 0.92, 41.5, 2, 1};
    const std::vector<double> fHom = assembleFeatures(hom, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_DOUBLE_EQ(fHom[6], 2.0);
    EXPECT_DOUBLE_EQ(fHom[7], 1.0);
    // Hemizygous call (20): one allele, one distinct size.
    const LocusFeatureContext hemi{3, 30, spanning, flanking, hq, 0.92, 41.5, 1, 1};
    const std::vector<double> fHemi = assembleFeatures(hemi, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_DOUBLE_EQ(fHemi[6], 1.0);
    EXPECT_DOUBLE_EQ(fHemi[7], 1.0);
}

TEST(GenotypeQualityFeatures, UnsetLocusContextFieldsMapToNaN)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    // Defaults: coverage NaN, numAlleles/numDistinctAlleles 0 (not set).
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_TRUE(std::isnan(f[6]));  // n_alleles
    EXPECT_TRUE(std::isnan(f[7]));  // n_distinct_alleles
    EXPECT_TRUE(std::isnan(f[19])); // coverage
}

TEST(GenotypeQualityFeatures, QuickGenotypingRegimeOmitsFlankDepths)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_EQ(f.size(), 27u); // QUICK_FEATURES has no flank-normalized depths
}

TEST(GenotypeQualityFeatures, RepeatPuritySentinelMapsToNaN)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    // Default ctx.referenceRepeatPurity (-1.0, not computed) and an aqm with the default
    // (unset) readRepeatPurity (-1.0) should both map to NaN, not the raw sentinel.
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};
    AlleleMetrics aqm = makeAqm();
    aqm.readRepeatPurity = -1.0;

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, aqm, GenotypingRegime::Quick);
    EXPECT_TRUE(std::isnan(f[25])); // reference_repeat_purity
    EXPECT_TRUE(std::isnan(f[26])); // read_repeat_purity
}

TEST(GenotypeQualityFeatures, SupportFracIsNaNWhenNoSpanningReads)
{
    const CountTable spanning; // empty
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_TRUE(std::isnan(f[17])); // support_frac undefined with 0 spanning reads
    EXPECT_DOUBLE_EQ(f[11], 0.0);   // spanning_total
    EXPECT_DOUBLE_EQ(f[14], 0.0);   // spanning_at_called
    EXPECT_DOUBLE_EQ(f[13], 0.0);   // flanking_total
    EXPECT_DOUBLE_EQ(f[18], 0.0);   // flanking_frac is 0, not NaN, with no reads of either kind
}

TEST(GenotypeQualityFeatures, ValuesKeepFullDoublePrecisionForTheModel)
{
    // The training side keeps float64 end to end, so inference must hand the model the full double.
    // An earlier pipeline downcast the training parquets to float32 and this assembler reproduced that
    // downcast; both are gone. 47.47 and 21/53 are not exactly representable in float32, so a
    // reinstated round-trip on this side alone shows up here.
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 32}});
    const CountTable flanking(std::map<int32_t, int32_t>{{25, 21}});
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq, 0.92, /*coverage=*/47.47, 2, 2};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_DOUBLE_EQ(f[19], 47.47);          // coverage, not its float32 image
    EXPECT_DOUBLE_EQ(f[18], 21.0 / 53.0);    // flanking_frac, not its float32 image
    EXPECT_NE(f[19], static_cast<double>(static_cast<float>(47.47)));
    EXPECT_NE(f[18], static_cast<double>(static_cast<float>(21.0 / 53.0)));
    EXPECT_DOUBLE_EQ(f[13], 21.0);           // flanking_total: integral, identical either way
    EXPECT_TRUE(std::isnan(assembleFeatures(
        LocusFeatureContext{3, 30, spanning, flanking, hq}, 0, 20, 18, 24, makeAqm(),
        GenotypingRegime::Quick)[19]));      // an unsupplied coverage is still NaN
}

TEST(GenotypeQualityFeatures, AlleleMetricsAreRoundedTo3DecimalsLikeTheEmittedJson)
{
    // The training parquets are parsed from JsonWriter's output, where every AlleleQualityMetrics
    // field goes through round3. Inference must round the same way, or it feeds the model a value up
    // to 5e-4 away from the one training saw. makeAqm()'s values are already 3-decimal, so this test
    // supplies longer ones.
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq, /*referenceRepeatPurity=*/0.9166666, 41.5, 2, 2};
    AlleleMetrics aqm = makeAqm();
    aqm.depth = 22.3333333;
    aqm.strandBiasBinomialPhred = 5.4321987;
    aqm.readRepeatPurity = 0.9594999;
    aqm.leftFlankNormalizedDepth = 0.7776666;

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, aqm, GenotypingRegime::FullSpanning);
    EXPECT_DOUBLE_EQ(f[20], (22.333));  // depth
    EXPECT_DOUBLE_EQ(f[22], (5.432));   // strand_bias_phred
    EXPECT_DOUBLE_EQ(f[25], (0.917));   // reference_repeat_purity
    EXPECT_DOUBLE_EQ(f[26], (0.959));   // read_repeat_purity
    EXPECT_DOUBLE_EQ(f[27], (0.778));   // left_flank_norm_depth (full-only)
    EXPECT_NE(f[20], (22.3333333));     // the unrounded value must NOT have reached the model
}

TEST(GenotypeQualityFeatures, CanonicalNamesMatchAssemblerOrder)
{
    const std::vector<std::string>& quick = featureNamesForGenotypingRegime(GenotypingRegime::Quick);
    const std::vector<std::string>& full = featureNamesForGenotypingRegime(GenotypingRegime::FullSpanning);
    ASSERT_EQ(quick.size(), 27u);
    ASSERT_EQ(full.size(), 29u);
    EXPECT_EQ(quick.front(), "motif_size");           // value index 0 in the asserts above
    EXPECT_EQ(quick[6], "n_alleles");                 // value index 6
    EXPECT_EQ(quick[7], "n_distinct_alleles");        // value index 7
    EXPECT_EQ(quick[13], "flanking_total");           // value index 13
    EXPECT_EQ(quick[18], "flanking_frac");            // value index 18
    EXPECT_EQ(quick[19], "coverage");                 // value index 19
    EXPECT_EQ(quick[20], "depth");                    // value index 20
    EXPECT_EQ(quick[24], "mean_deleted_bases");       // value index 24
    EXPECT_EQ(quick[25], "reference_repeat_purity");  // value index 25
    EXPECT_EQ(quick.back(), "read_repeat_purity");    // value index 26
    EXPECT_EQ(full[27], "left_flank_norm_depth");    // full-only, value index 27
    EXPECT_EQ(full[28], "right_flank_norm_depth");   // full-only, value index 28
    EXPECT_TRUE(std::equal(quick.begin(), quick.end(), full.begin())); // full is quick + 2
    EXPECT_EQ(featureNamesForGenotypingRegime(GenotypingRegime::FullNonspanning), full); // both full genotyping_regimes share the list
}
