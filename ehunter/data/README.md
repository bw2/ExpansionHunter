# Genotype quality model data

`genotype_quality_model_from_HG002_and_CHM1_CHM13.20260701.json.gz` specifies a pretrained gradient-boosted decision tree model 
that ExpansionHunter uses to generate the `PredictedLengthCorrectionFactor`,
`pOk`, `pTooShort`, and `pTooLong` fields within the output JSON `AlleleQualityMetrics`  section.

To override the default model, specify  `--genotype-quality-model PATH` .

The pipeline used to train this model is available in:

https://github.com/bw2/ExpansionHunterGenotypeQualityModel

The truth data was generated from T2T assemblies as described in [[Weisburd 2023](https://pubmed.ncbi.nlm.nih.gov/37214979/)] using the pipelines in 

https://github.com/broadinstitute/str-truth-set-v2



