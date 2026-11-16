# Multi-sample utilities

Scripts for turning many samples' ExpansionHunter output JSONs into one cohort table, after leaving out the samples
and loci that fail quality control (QC). The intended order is:

1. **QC stats**: `generate_multisample_qc_stats.py` combines the JSONs of every sample into one QC stats file.
2. **Choose what to exclude**: `start_multisample_qc_dashboard.py` opens that file in an interactive dashboard, where
   you set thresholds, review outliers, and export the samples and loci to exclude.
3. **Combine**: `generate_multisample_tsv.py` combines the same JSONs into one table with a row per repeat, leaving out
   the exported samples and loci.

Each script's `--help` describes its options and output columns in full.

## 1. Compute QC stats: `generate_multisample_qc_stats.py`

```
python3 generate_multisample_qc_stats.py --threads 8 --show-progress-bar -o cohort --input-paths-file json_paths.txt
```

- **Input**: ExpansionHunter output JSONs (`.json`, `.json.gz`, `.json.bgz`) from any of its analysis modes, as local
  paths or `gs://` paths, which are streamed rather than downloaded. `--input-paths-file` takes a file with one path
  per line.
- **Output**: `<prefix>.qc_stats.tsv.gz`, a gzipped text file with:
  - one row per repeat: totals across the samples, such as how many samples were genotyped and the allele sizes seen
  - one row per sample: metrics computed from that sample's JSON, such as call rate, coverage, fragment length, chrX and
    chrY coverage, heterozygosity, and read-based estimates of contamination and allele balance
- **Memory** does not grow with the number of samples. The inputs are read side by side and merged one repeat at a
  time. With `--threads N`, N processes each merge a share of the inputs, and their results are then merged.
- **Large cohorts**: a QC stats file can itself be an input, so you can run one job per batch of JSONs and then combine
  the batch files. The result is the same as one run over all the JSONs:

  ```
  python3 generate_multisample_qc_stats.py --threads 16 -o batch1 --input-paths-file batch1_json_paths.txt
  python3 generate_multisample_qc_stats.py --threads 16 -o batch2 --input-paths-file batch2_json_paths.txt
  python3 generate_multisample_qc_stats.py -o cohort batch1.qc_stats.tsv.gz batch2.qc_stats.tsv.gz
  ```

- A sample ID that appears in more than one input is an error.

## 2. Choose samples and loci to exclude: `start_multisample_qc_dashboard.py`

```
python3 start_multisample_qc_dashboard.py cohort.qc_stats.tsv.gz
```

The script reads the QC stats file and starts a local web server. It prints a URL to open in your browser.

- **What it shows**: one plot per metric, for samples and for loci. Each plot has a point per sample or locus and
  faint suggested thresholds, and points beyond a suggested threshold are highlighted in orange.
- **Thresholds**: drag a plot's red threshold lines to exclude the points beyond them, or click "Use suggested
  thresholds". A sample or locus is excluded if any plot's threshold excludes it. "Reset thresholds" turns a plot's
  thresholds off, and "Reset all sample thresholds" and "Reset all locus thresholds", under the excluded samples and
  loci tables, turn off all of them. "Reset manually excluded samples" (or loci) undoes every choice made by hand.
  - A plot only has thresholds on the sides where its values point to a problem. For example, the sample and locus
    call rates have only a lower threshold.
  - The median fragment length and the observed vs. expected (HWE) heterozygosity plots have no thresholds. Each
    plot's help says why.
  - Read lengths are shown as a table of the number of samples with each read length. Its two lines snap to the
    boundaries between rows and exclude the rows above the upper line or below the lower line.
- **Overriding the thresholds**: click, shift-click, or box and lasso select points, then right-click and choose
  "Toggle whether these samples are on the Excluded list" (or this sample, this locus, these loci). This includes or
  excludes those points regardless of the thresholds.
- **Reasons**: each exclusion gets a reason, such as "Call rate 0.62 < 0.9". You can edit it in the exclusion tables
  (pencil icon) before exporting.
- **Exports**: the "Export excluded samples" and "Export excluded loci" buttons write:
  - `<prefix>.samples_to_exclude.tsv`, with columns `SampleId`, `ReasonsForExclusion`
  - `<prefix>.loci_to_exclude.tsv`, with columns `LocusId`, `VariantId`, `ReasonsForExclusion`
- **Saved state**: thresholds, toggles, edited reasons and selections are saved to `<prefix>.qc_dashboard_state.json`
  after every change. Starting the dashboard again with the same prefix restores them.
- **Default prefix**: the QC stats file's path without `.qc_stats.tsv.gz`. Use `-o` to change it.

Things to know:

- **Loci with more than one repeat are left out of the dashboard.** These are compound locus definitions, such as
  HTT's CAG and CCG repeats. `generate_multisample_tsv.py --loci-to-exclude` leaves out whole loci, so excluding one
  of their repeats would also drop the others. To exclude such a locus, add its LocusId to the loci TSV by hand.
- **The locus metrics include every sample in the QC stats file**, including samples you exclude in the dashboard. To
  recompute them without the excluded samples:
  1. Rerun step 1 with `--samples-to-exclude cohort.samples_to_exclude.tsv` and a new `-o`.
  2. Open that file in the dashboard and choose the loci there.
  3. Use the samples TSV from the first dashboard and the loci TSV from the second in step 3.
- **The server listens on 127.0.0.1 only** and answers only requests that carry the token in the printed URL. To use
  the dashboard on a remote machine, forward the port over SSH (`ssh -L <port>:127.0.0.1:<port> <host>`).
- **Requirements**: the `numpy` and `plotly` packages.

## 3. Combine the samples: `generate_multisample_tsv.py`

```
python3 generate_multisample_tsv.py --threads 8 --show-progress-bar -o cohort --input-paths-file json_paths.txt \
    --samples-to-exclude cohort.samples_to_exclude.tsv --loci-to-exclude cohort.loci_to_exclude.tsv
```

- **Output**: `<prefix>.tsv.gz`, with one row per repeat. It has allele size histograms across the samples, and
  optional columns for outlier sample IDs, motif composition and consensus sequences.
- **`--samples-to-exclude`** takes a TSV with a `SampleId` column (like the dashboard's export) or one ID per line.
  - It can only leave samples out of JSON inputs, since an earlier run's TSV has its samples merged already.
  - When you combine a large cohort in batches, pass it to each batch run.
- **`--loci-to-exclude`** takes a TSV with a `LocusId` column or one LocusId per line. It matches loci by LocusId only
  and leaves out every repeat of each listed locus.
- **Warnings**: listed samples or loci that aren't in the inputs are reported.
- **Large cohorts**: as in step 1, the outputs can be combined again, so a cohort can be run in batches:

  ```
  python3 generate_multisample_tsv.py -o batch1 --samples-to-exclude ... --loci-to-exclude ... --input-paths-file batch1_json_paths.txt
  python3 generate_multisample_tsv.py -o cohort batch1.tsv.gz batch2.tsv.gz
  ```

## Other files

- `expansion_hunter_json_streaming_utils.py`: shared code used by `generate_multisample_qc_stats.py` and
  `generate_multisample_tsv.py`. It reads local and `gs://` inputs, puts each JSON's repeats in genomic order, splits
  the inputs among processes, and reads the exclusion lists.
- `multisample_qc_dashboard.html`: the dashboard page that `start_multisample_qc_dashboard.py` serves.
- `*_tests.py`: unit tests. Run them from this folder with
  `python3 -m unittest generate_multisample_qc_stats_tests start_multisample_qc_dashboard_tests generate_multisample_tsv_tests`.

The scripts need the `ijson` and `tqdm` packages. Reading `gs://` paths also needs the `google-cloud-storage` package
and application default credentials (`gcloud auth application-default login`).
