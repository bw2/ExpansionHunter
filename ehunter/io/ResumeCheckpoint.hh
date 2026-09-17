//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
//

// --resume support.
//
// An interrupted ExpansionHunter run leaves nothing usable behind: seeking/streaming modes write their
// output only at the end, and low-mem/optimized streaming with --threads > 1 produces the final files
// only once every worker has joined. So --resume maintains its own crash-survivable record of which loci
// are finished, in three files alongside the output:
//
//   <output>.json[.gz].unfinished   every finished locus's JSON record, in completion order (unsorted)
//   <output>.vcf[.gz].unfinished    every finished locus's VCF lines, in completion order (unsorted)
//   <prefix>.processed_loci.unfinished   one "<LocusId>\t<jsonRecords>\t<vcfLines>" line each
//
// The first two are written in the same format as a normal (truncated) ExpansionHunter output file, so
// they can be inspected directly. When compressed, each flush closes the current gzip member, which makes
// the file complete and readable on disk after every flush; multi-member gzip is understood by gzip/zcat,
// python, and boost.
//
// The third file is what makes resuming exact. A LocusId is appended to it only after both record files
// have been flushed past that locus, so it is by construction the set of loci durably present in both,
// which is how the "trim both files back to their common prefix" rule is enforced. It also records loci
// that produce no output record at all (dropped by --skip-hom-ref / --skip-missing-genotypes, or failed
// during genotyping), which would otherwise be re-genotyped on every resume. Its per-locus record counts
// let a resume tell a locus that produced no record apart from one whose records are missing because the
// file lost buffered writes.
//
// On resume the two record files are trimmed to exactly the processed set and demultiplexed into the
// per-contig temp files the normal end-of-run merge consumes, so the final output is produced by the
// existing merge code and its bytes are unchanged.

#pragma once

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <vector>

#include "core/Parameters.hh"
#include "locus/LocusSpecification.hh"

namespace ehunter
{

// The output text a locus produced, captured while it is written to the genotyping temp files so the
// checkpoint can replay it verbatim. Both members stay empty for a locus that produced no record at all
// (dropped by --skip-hom-ref / --skip-missing-genotypes, or failed during genotyping); such a locus is
// still checkpointed, so that it is not genotyped again.
struct CapturedLocusOutput
{
    std::string json;
    std::string vcf;
};

struct ResumeCheckpointPaths
{
    std::string json;
    std::string vcf;
    std::string processedLoci;
};

// A checkpoint that cannot be parsed at all, as opposed to one that parses but belongs to a different
// run. The caller warns and starts from the beginning for the former, and stops with an error for the
// latter rather than silently discarding an interrupted run's work.
class ResumeCheckpointCorruptError : public std::runtime_error
{
public:
    explicit ResumeCheckpointCorruptError(const std::string& message)
        : std::runtime_error(message)
    {
    }
};

ResumeCheckpointPaths resumeCheckpointPaths(const OutputPaths& outputPaths);

// The same output paths with the opposite --compress-output-files setting, used to notice a checkpoint
// left by a run that had -z the other way round (its record files carry, or lack, the .gz suffix).
OutputPaths otherCompressionOutputPaths(const OutputPaths& outputPaths);

// Paths of the plain-text per-slice temp files the genotyping writers target. The per-contig form is what
// --threads > 1 has always used; the single-slice form is used at --threads 1 under --resume, where one
// writer covers every contig in genomic order.
std::string contigTempJsonPath(const std::string& outputPrefix, std::int32_t contigIndex);
std::string contigTempVcfPath(const std::string& outputPrefix, std::int32_t contigIndex);
std::string singleSliceTempJsonPath(const std::string& outputPrefix);
std::string singleSliceTempVcfPath(const std::string& outputPrefix);

// True when all three checkpoint files are present.
bool resumeCheckpointExists(const ResumeCheckpointPaths& paths);

void removeCheckpointFiles(const ResumeCheckpointPaths& paths);

// Remove every resume temp file for this output prefix: the single-slice pair and every per-contig pair,
// not just the layout the current run used. A run that resumes an interrupted one with a different
// --threads setting rebuilds the other layout, so cleaning up only its own would leave those behind.
void removeResumeTempFiles(const std::string& outputPrefix);

// A JSON object holding every parameter that has to match for a checkpoint to be resumable: the sample,
// the analysis mode, the build, the catalog, and each flag that changes which loci are emitted or what is
// written for them. --threads is deliberately excluded, since temp files are keyed by contig index rather
// than by worker, so a run may be resumed with a different thread count.
std::string buildRunSignature(const ProgramParameters& params);

// Throws when the catalog holds duplicate LocusIds, which resume cannot tell apart.
void assertCatalogIsResumable(const LocusDescriptionCatalog& catalog);

// Returns a human-readable description of the first field that differs between two signatures, or an
// empty string when they match.
std::string describeSignatureMismatch(const std::string& expected, const std::string& found);

// The ResumeInfo signature stored in an existing checkpoint's JSON header. Throws
// ResumeCheckpointCorruptError when the checkpoint holds no readable signature.
std::string readCheckpointRunSignature(const ResumeCheckpointPaths& paths);

// Key used in ResumeLoadResult::rebuiltSlices for the single-slice (--threads 1) temp files, which cover
// every contig rather than one.
constexpr std::int32_t kSingleSliceKey = -1;

struct RebuiltSlice
{
    // Whether this slice's rebuilt JSON temp holds at least one record, which decides whether the next
    // record appended to it needs a ", " separator.
    bool hasJsonRecords = false;
};

struct ResumeLoadResult
{
    // Loci that must not be genotyped again.
    std::unordered_set<std::string> doneLocusIds;
    // Slices whose temp files were rebuilt from the checkpoint, keyed by contig index (or kSingleSliceKey).
    // The caller opens these in append mode and must include them in the end-of-run merge, even when the
    // filtered catalog no longer has any loci left on that contig.
    std::map<std::int32_t, RebuiltSlice> rebuiltSlices;
    // Whether the trimmed JSON checkpoint holds at least one record (decides whether the next record
    // appended to it needs a ", " separator).
    bool hasExistingJsonRecords = false;
};

// Reads the checkpoint, trims both record files to the loci listed in processed_loci, and rebuilds the
// genotyping temp files from the kept records. With `demultiplexPerContig` the records are split into
// <outputPrefix>.contig<N>.{json,vcf}; otherwise they all go to the single-slice temp, which is correct
// only when the checkpoint is already in genomic order (--threads 1).
//
// Throws ResumeCheckpointCorruptError when the checkpoint cannot be parsed, and plain std::runtime_error
// when it parses but does not belong to this run (a locus that is not in the catalog), which the caller
// must report rather than silently discard.
ResumeLoadResult loadResumeCheckpoint(
    const ResumeCheckpointPaths& paths, const SampleParameters& sampleParams,
    const LocusDescriptionCatalog& catalog, const std::string& outputPrefix, bool demultiplexPerContig);

// Drops the loci that are already done, preserving the order of the rest.
void removeCompletedLoci(LocusDescriptionCatalog& catalog, const std::unordered_set<std::string>& doneLocusIds);

// Receives each finished locus's output text and appends it to the checkpoint files from a background
// thread, so genotyping threads are never blocked on checkpoint I/O.
class ResumeCheckpointWriter
{
public:
    // `append` continues an existing (already trimmed) checkpoint instead of creating a new one.
    ResumeCheckpointWriter(
        ResumeCheckpointPaths paths, const SampleParameters& sampleParams, const std::string& runSignature,
        bool append, bool hasExistingJsonRecords, std::size_t abortAfterLoci);
    ~ResumeCheckpointWriter();

    ResumeCheckpointWriter(const ResumeCheckpointWriter&) = delete;
    ResumeCheckpointWriter& operator=(const ResumeCheckpointWriter&) = delete;

    // Called from the genotyping threads once a locus is fully written to its temp files. Either text may
    // be empty; a locus that produced no output at all is still recorded, so it is not genotyped again.
    void recordLocus(std::string locusId, std::string jsonText, std::string vcfText);

    // Drain the queue, flush, and join the background thread. Idempotent. Rethrows a write failure from
    // the background thread.
    void close();

private:
    struct QueuedLocus
    {
        std::string locusId;
        std::string jsonText;
        std::string vcfText;
    };

    void run();
    void writeBatch(const std::vector<QueuedLocus>& batch);

    ResumeCheckpointPaths paths_;
    bool compressJson_;
    bool compressVcf_;
    bool firstJsonRecord_;
    std::size_t abortAfterLoci_;
    std::size_t checkpointedLocusCount_ = 0;

    std::mutex mutex_;
    std::condition_variable queueNotEmpty_;
    std::condition_variable queueNotFull_;
    std::deque<QueuedLocus> queue_;
    bool stopRequested_ = false;
    std::exception_ptr writeError_;
    std::thread thread_;
    bool closed_ = false;
};

}
