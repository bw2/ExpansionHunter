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
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, /*rank=*/0, /*eh=*/20, /*ciStart=*/18,
                                                   /*ciEnd=*/24, makeAqm(), GenotypingRegime::FullSpanning);

    ASSERT_EQ(f.size(), 22u);
    EXPECT_DOUBLE_EQ(f[0], 3.0);             // motif_size
    EXPECT_DOUBLE_EQ(f[1], 10.0);            // num_repeats_in_reference = 30/3
    EXPECT_DOUBLE_EQ(f[2], 30.0);            // ref_size_bp
    EXPECT_DOUBLE_EQ(f[3], 20.0);            // eh
    EXPECT_DOUBLE_EQ(f[4], 10.0);            // eh_minus_ref = 20 - 10
    EXPECT_DOUBLE_EQ(f[5], 0.0);             // allele_rank
    EXPECT_DOUBLE_EQ(f[6], 6.0);             // ci_width = 24 - 18
    EXPECT_DOUBLE_EQ(f[7], (4.0 - 2.0) / 7.0); // ci_asymmetry = ((24-20)-(20-18))/(6+1)
    EXPECT_DOUBLE_EQ(f[8], 6.0 / 21.0);      // ci_over_eh = 6/(20+1)
    EXPECT_DOUBLE_EQ(f[9], 9.0);             // spanning_total = 8 + 1
    EXPECT_DOUBLE_EQ(f[10], 5.0);            // hq_unamb_total
    EXPECT_DOUBLE_EQ(f[11], 8.0);            // spanning_at_called = count at 20
    EXPECT_DOUBLE_EQ(f[12], 1.0);            // spanning_above_called = count of 22
    EXPECT_DOUBLE_EQ(f[13], 2.0);            // flanking_above_called = count of 25
    EXPECT_DOUBLE_EQ(f[14], 8.0 / 9.0);      // support_frac
    EXPECT_DOUBLE_EQ(f[15], 30.5);           // depth
    EXPECT_DOUBLE_EQ(f[16], 12.0);           // hq_unambiguous_reads
    EXPECT_DOUBLE_EQ(f[17], 1.2);            // strand_bias_phred
    EXPECT_DOUBLE_EQ(f[18], 0.3);            // mean_inserted_bases
    EXPECT_DOUBLE_EQ(f[19], 0.1);            // mean_deleted_bases
    EXPECT_DOUBLE_EQ(f[20], 0.9);            // left_flank_norm_depth
    EXPECT_DOUBLE_EQ(f[21], 1.1);            // right_flank_norm_depth
}

TEST(GenotypeQualityFeatures, QuickGenotypingRegimeOmitsFlankDepths)
{
    const CountTable spanning(std::map<int32_t, int32_t>{{20, 8}});
    const CountTable flanking;
    const CountTable hq(std::map<int32_t, int32_t>{{20, 5}});
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_EQ(f.size(), 20u); // QUICK_FEATURES has no flank-normalized depths
}

TEST(GenotypeQualityFeatures, SupportFracIsNaNWhenNoSpanningReads)
{
    const CountTable spanning; // empty
    const CountTable flanking;
    const CountTable hq;
    const LocusFeatureContext ctx{3, 30, spanning, flanking, hq};

    const std::vector<double> f = assembleFeatures(ctx, 0, 20, 18, 24, makeAqm(), GenotypingRegime::Quick);
    EXPECT_TRUE(std::isnan(f[14])); // support_frac undefined with 0 spanning reads
    EXPECT_DOUBLE_EQ(f[9], 0.0);    // spanning_total
    EXPECT_DOUBLE_EQ(f[11], 0.0);   // spanning_at_called
}

TEST(GenotypeQualityFeatures, CanonicalNamesMatchAssemblerOrder)
{
    const std::vector<std::string>& quick = featureNamesForGenotypingRegime(GenotypingRegime::Quick);
    const std::vector<std::string>& full = featureNamesForGenotypingRegime(GenotypingRegime::FullSpanning);
    ASSERT_EQ(quick.size(), 20u);
    ASSERT_EQ(full.size(), 22u);
    EXPECT_EQ(quick.front(), "motif_size");          // value index 0 in the asserts above
    EXPECT_EQ(quick[15], "depth");                   // value index 15
    EXPECT_EQ(quick.back(), "mean_deleted_bases");   // value index 19
    EXPECT_EQ(full[20], "left_flank_norm_depth");   // full-only, value index 20
    EXPECT_EQ(full[21], "right_flank_norm_depth");  // full-only, value index 21
    EXPECT_TRUE(std::equal(quick.begin(), quick.end(), full.begin())); // full is quick + 2
    EXPECT_EQ(featureNamesForGenotypingRegime(GenotypingRegime::FullNonspanning), full); // both full genotyping_regimes share the list
}
