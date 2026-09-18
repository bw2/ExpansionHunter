# `--resume`: design notes

Status: implemented. Written 2026-09-17.

User-facing documentation lives in [docs/03_Usage.md](docs/03_Usage.md#resuming-an-interrupted-run).
This file records why the feature is built the way it is.

Goal: let an ExpansionHunter run that was interrupted (SIGKILL, OOM kill, Ctrl-C, VM preemption,
crash) be restarted so it genotypes only the loci it had not finished, instead of starting over.

## 1. Why the obvious approach does not work

Neither place that could hold "work done so far" is usable as-is:

- **`seeking` and `streaming` modes** build `SampleFindings` in memory and call `writeToFile` once, at
  the very end (`ehunter/app/ExpansionHunter.cpp`). An interrupted run has written nothing, so these
  modes are out of scope: `--resume` there logs a warning and changes nothing.
- **`low-mem-streaming` / `optimized-streaming` at `--threads > 1`** write per-contig temp files
  `<prefix>.contig<N>.json` / `.vcf` and only produce the final `<prefix>.json[.gz]` / `.vcf[.gz]`
  after every worker has joined, via `mergeRegionJsonFiles` / `mergeRegionVcfFiles`. An interrupted
  run has no final output at all, and `TempFileRemover` deletes every temp on any exit path,
  including exceptions.
- Only **`--threads 1`** in a streaming mode writes records straight to the final files as it goes.

So resume needs its own crash-survivable record of completed loci.

## 2. Design

### 2.1 Checkpoint files

With `--resume`, a background thread maintains three files. They are written during genotyping, read
at startup by a later run, and deleted once the final outputs are complete.

| File | Contents | Compressed |
| --- | --- | --- |
| `<outputPaths.json()>.unfinished` | Every finished locus's JSON record, in completion order (unsorted) | Follows `-z` |
| `<outputPaths.vcf()>.unfinished` | Every finished locus's VCF lines, in completion order (unsorted) | Follows `-z` |
| `<prefix>.processed_loci.unfinished` | `<LocusId>\t<jsonRecords>\t<vcfLines>` per line, appended after that locus is durable in both files above | No |

Deriving the first two names from `outputPaths.json()` / `outputPaths.vcf()` makes them come out as
`<prefix>.json.gz.unfinished` with `-z` and `<prefix>.json.unfinished` without, matching whatever the
run's real output is.

The two record files use the same format the writers already produce (JSON: header plus
`"LocusResults"` records separated by `", "`; VCF: header plus body lines), so one scanner reads a
checkpoint, a per-contig temp, and a final output file alike, and the files can be inspected directly.

When compressed, each flush closes the current gzip member, making the file a multi-member gzip. That
is valid gzip: `gzcat`, `zcat`, python's `gzip`, and boost's `gzip_decompressor` all handle it
(`boost/iostreams/filter/gzip.hpp:520-530` resets to `s_start` when bytes remain after a footer). The
effect is that the file on disk is complete and readable after every flush rather than only at the end.

### 2.2 Why the third file exists

Three categories of loci produce **no record in either output**: loci dropped by `--skip-hom-ref`,
loci dropped by `--skip-missing-genotypes`, and loci that threw during genotyping
(`LocusOutput::Kind::kError` and `kFilteredOut` in `writeOutput`). Without a separate list of finished
loci, resume would re-genotype every one of them. On a large catalog run with `--skip-hom-ref` that is
most of the catalog, which would make the feature nearly worthless.

It also makes the JSON/VCF consistency rule exact instead of heuristic. The checkpoint thread appends
a locus to `processed_loci` only **after** both record files have been flushed past it, so the list is
by construction the set of loci durably present in both. Resume then trims both record files to
exactly that set, which is the "trim both back to their common prefix" behavior, handled by
construction rather than by diffing the two files.

That matters because the JSON and VCF writers have independent `ofstream` and gzip buffers: after a
kill they are essentially never at the same locus, so a rule of "if the two files disagree, start
over" would fire on nearly every real crash and throw away the whole run.

The per-locus counts record how many records the locus wrote to each file. They let a resume tell "this
locus produced no record" apart from "this locus's records are missing from the file", so a checkpoint
that lost buffered record writes (a machine crash rather than a killed process) is detected and those
loci are genotyped again instead of vanishing from the output. The VCF side is a line count rather than
a flag because a multi-variant locus writes one line per variant, and losing one of several lines has to
be caught too.

### 2.3 Feeding the checkpoint

The tee is driven from the call sites that write a locus, not from inside the writers, because a locus
may write to one file, both, or neither:

| Call site | JSON | VCF |
| --- | --- | --- |
| `writeOutput` / `kGenotyped` | record | one line per variant |
| `writeOutput` / `kNoCoverage` -> `writeZeroCoverageRecord` | record | lines |
| `writeOutput` / `kHeuristicOnlySkip` | skipped record | nothing |
| `writeOutput` / `kError`, `kFilteredOut` | nothing | nothing |
| `processLocusFast` (fast path) | record | lines |

`IterativeJsonWriter::addRecord` / `addSkippedRecord` and `IterativeVcfWriter::addRecord` /
`addRecords` take an optional `capturedText` out-parameter, so the exact bytes they write can be handed
to the checkpoint without re-formatting. Each call site then makes one
`checkpoint->recordLocus(locusId, json, vcf)` call, with either or both possibly empty.

The checkpoint thread drains a bounded queue, appends to both record files, flushes them (closing the
gzip member), then appends the batch's locus ids to `processed_loci` and flushes that. No `fsync`: a
process-level kill leaves OS buffers intact, and the flag-based verification above covers the cases
that do not.

### 2.4 Resume at startup

When `--resume` finds existing checkpoint files (`htsLowMemStreamingSampleAnalysis`):

1. **Validate the run signature** (section 2.5). On mismatch, fail with a clear error rather than
   silently discarding the checkpoint. If the checkpoint cannot be read at all, warn, discard it, and
   start from the beginning.
2. Read `processed_loci`.
3. **Verify** it against the record files, truncating at the first locus whose declared records are
   missing (section 2.2).
4. **Restrict to what this run's slices can append to** (`firstUnusableProcessedLocus`). Each slice's
   rebuilt temp is appended to by a genotyping writer, which can only add records at the end, so what a
   slice needs is that the finished loci are a prefix of *the order that slice emits records in*.

   That order is **not** catalog order. A locus is released for genotyping once the reads pass the end of
   its flank window, so a locus nested inside a wider one finishes first even though it comes second in
   the position-sorted catalog. The checkpoint records the order the interrupted run emitted in, so within
   one slice any leading run of it is already a valid prefix and nothing is checked: at `--threads > 1`,
   where a slice is one contig, every entry is reusable.

   What does need checking is the slice layout changing between runs. At `--threads 1` one slice covers
   every contig in a single coordinate sweep, so its rebuilt temp must be a prefix of that sweep: a contig
   may only be left behind once every one of its loci is finished. A `--threads > 1` checkpoint generally
   fails that (it can finish chr3 while chr1 is still running), and such a resume is **refused before any
   file is rewritten**. Trimming instead would cut the checkpoint down to the reusable part, destroying the
   rest of the interrupted run's work for good, and the advice to re-run with `--threads > 1` would be
   useless by the time it was read.
5. Rewrite all three checkpoint files to the surviving set, via a sibling temp plus rename, so they
   agree exactly and later appends start from a clean stream.
6. Rebuild the genotyping temps from the kept records: `<prefix>.contig<N>.{json,vcf}` at
   `--threads > 1`, or `<prefix>.resume_part.{json,vcf}` at `--threads 1`. They are written directly
   rather than through the writer classes so they stay header-plus-records: the genotyping writer
   appends to them and closes them, which is what adds the JSON footer. A contig whose loci are all
   finished has no slice, so nothing would close it; the orchestrator closes those explicitly before
   the merge.
7. Filter the catalog by the finished set and log `Resuming: N of M loci already genotyped, K remaining`.

### 2.5 Run signature

A `ResumeInfo` object in the JSON checkpoint's header records everything that would make two runs'
records incompatible: sample id, sex, analysis mode, `kCommitSha`, catalog path and size,
genotype-quality model version, and the flags that change what is emitted. A mismatch is a hard error
naming the field, with the remedy. Erroring is deliberate: silently restarting would throw away hours of
work over a typo in a path.

The catalog is keyed by size rather than by modification time on purpose: a pipeline that localizes the
catalog fresh for every attempt, as a cloud batch job does, would otherwise fail every resume over a
byte-identical file.

Two things are rejected up front rather than through the signature. A catalog with duplicate `LocusId`s
cannot be resumed at all, since resume tells loci apart by their id and finishing one occurrence would
mark both done (such a catalog is already malformed: the JSON output is keyed by `LocusId`, so the
entries collide into one record even without `--resume`). And a checkpoint written with the opposite
`-z` setting is reported explicitly, because the two record files' names carry the output's `.gz`
suffix, so flipping `-z` would otherwise look like "no checkpoint here" and quietly start over.

`--threads` is deliberately **not** part of the signature. Temp files are keyed by contig index rather
than by worker, so a `--threads 1` checkpoint is reusable at any thread count, and a `--threads > 1`
checkpoint at any count above 1: a `--threads 8` run can be resumed with `--threads 2`. The one
combination step 4 refuses is a `--threads > 1` checkpoint resumed at `--threads 1`, and only when that
run actually finished contigs out of order.

### 2.6 Finalization

- `--threads > 1`: after the workers join, `mergeRegionJsonFiles` / `mergeRegionVcfFiles` run over the
  per-contig temps (now including the rebuilt ones). **The final output bytes are unchanged**, verified
  against a build of the base revision across every analysis mode and thread count. The merge itself was
  reworked to locate each region file's header and footer by reading only those two windows and then
  stream the body through in chunks: it previously read each region file into memory whole, which under
  `--resume` at `--threads 1` would have meant the entire genome-wide output at once. Then the temps and
  the checkpoint files are deleted.
- `--threads 1` under `--resume`: the same merge runs over the single `resume_part` temp. Without
  `--resume` this path is untouched and still writes the final files directly.

### 2.7 Resumed output is byte-identical

The only order-dependent randomness in the analysis path is the `--max-depth` reservoir sampling, and
its seed is an FNV-1a hash of the `LocusId` (`ehunter/sample/HtsLowMemStreamingSampleAnalysis.cpp:94`),
not a global RNG. Nothing else in `core/`, `locus/`, `sample/`, `io/`, `genotyping/` or `alignment/`
calls `rand()`. So a locus genotyped in a resumed process produces the same bytes as in an
uninterrupted run, and the tests assert it.

## 3. Known limitations

- Resuming skips the genotyping work but not the read scanning: the BAM/CRAM is still read from the
  beginning, and the far-away-mate prepass still runs. On large catalogs genotyping dominates.
- `--enable-bamlet-output` is rejected with `--resume`: the bamlet is opened with `hts_open(path,
  "wb")` on every run and carries no record of which loci it covers, so a resumed run would replace it
  with one holding only that run's loci.
- Every record is written three times over the run (per-slice temp, checkpoint, final merge). At
  `--threads > 1` the per-slice temps are the per-contig files the merge already used; at `--threads 1`,
  `--resume` adds one temp covering the whole output, where without it the writers went straight to the
  final files. The temps are never compressed, so under `-z` the peak extra disk is roughly one
  uncompressed copy of the output plus one compressed copy.
- A run that fails after writing checkpoint files and is then re-run *without* `--resume` leaves the
  `.unfinished` and `.resume_part.*` files behind; they are only cleaned up by a completing `--resume`
  run.
- The catalog is fingerprinted by size, not content, so an edit that happens to preserve the file size
  is not caught by the signature check (a locus that is no longer in the catalog still is).

## 4. Where the code lives

| File | Role |
| --- | --- |
| `ehunter/io/ResumeCheckpoint.{hh,cpp}` | checkpoint writer and thread; tolerant scanner for plain, single-member gzip, multi-member gzip and truncated input; verification, trimming and slice rebuilding; run signature |
| `ehunter/io/IterativeJsonWriter.{hh,cpp}`, `IterativeVcfWriter.{hh,cpp}` | `capturedText` out-parameter, append-mode construction, extracted document headers |
| `ehunter/sample/HtsLowMemStreamingSampleAnalysis.cpp` | resume at startup, checkpoint lifetime, per-slice temp wiring, finalization |
| `ehunter/sample/HtsLowMemStreamingHelpers.{hh,cpp}` | capture at the fast-path and zero-coverage call sites |
| `ehunter/io/StreamingOutputMerge.cpp` | merge reworked to stream region files rather than read them whole; write failures now propagate |
| `ehunter/io/ParameterLoading.cpp`, `ehunter/core/Parameters.hh` | `--resume`, `--internal-abort-after-loci`, validation |
| `ehunter/tests/ResumeTest.cpp` | unit tests for the scanner, verification, ordering rules and signature |
| `resume_tests.py` | end-to-end tests: interrupt, resume, compare against an uninterrupted run |

`--internal-abort-after-loci N` is a test-only hook in the "Internal options" section: it calls
`std::_Exit` once N loci have been checkpointed, leaving the process exactly as an external kill would.
