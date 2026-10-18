# Genotype quality model data

`genotype_quality_model_from_HG002_and_CHM1_CHM13_plus50.20261017.json.gz` specifies a pretrained gradient-boosted decision tree model 
that ExpansionHunter uses to generate the `PredictedLengthCorrectionFactor`,
`pOk`, `pTooShort`, and `pTooLong` fields within the output JSON `AlleleQualityMetrics` section.

To override the default model, specify `--genotype-quality-model PATH`.

The pipeline used to train this model is available in:

https://github.com/bw2/ExpansionHunterGenotypeQualityModel

It was trained on HG002 (10x, 20x and 31x), CHM1_CHM13 (46x) and 50 1kGP/HPRC samples.

The truth data was generated from T2T assemblies as described in [[Weisburd 2023](https://pubmed.ncbi.nlm.nih.gov/37214979/)] using the pipelines in 

https://github.com/broadinstitute/str-truth-set-v2

The truth genotypes count repeat insertions that DipCall's left-alignment placed up to one motif length before a
locus (str-analysis bd68a2a). The earlier `..._plus50.20261006` model was trained on truth that reported the
reference length there.

The full genotyper's models don't use the allele's `StrandBiasBinomialPhred`. Until 2026-10-17, the full genotyper in
streaming and optimized-streaming modes counted the strand of whichever read of a pair it decoded last instead of read
1, and the training data was made that way, so the previous model (`..._plus50.20261008`) learned values that the fixed
binary no longer produces. The fast-path model still uses strand bias, which the fast path computes with its own code.



