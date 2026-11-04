# Genotype quality model data

`genotype_quality_model_from_HG002_and_CHM1_CHM13_plus50.20261104.json.gz` specifies a pretrained gradient-boosted decision tree model 
that ExpansionHunter uses to generate the `PredictedLengthCorrectionFactor`,
`pOk`, `pTooShort`, and `pTooLong` fields within the output JSON `AlleleQualityMetrics` section.

To override the default model, specify `--genotype-quality-model PATH`.

The pipeline used to train this model is available in:

https://github.com/bw2/ExpansionHunterGenotypeQualityModel

It was trained on HG002 (10x, 20x and 31x), CHM1_CHM13 (46x) and 50 1kGP/HPRC samples.

The truth data was generated from T2T assemblies as described in [[Weisburd 2023](https://pubmed.ncbi.nlm.nih.gov/37214979/)] using the pipelines in 

https://github.com/broadinstitute/str-truth-set-v2

Also, simulated expanded alleles are added at the known disease loci, because healthy genomes
don't carry alleles in the pathogenic range. The simulated set covers carriers at every disease locus, including
simulations of two expanded alleles at the autosomal recessive loci.



