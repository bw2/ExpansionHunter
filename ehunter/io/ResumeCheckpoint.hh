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
// An interrupted ExpansionHunter run leaves nothing usable behind by default: seeking/streaming modes write
// their output only at the end, and low-mem/optimized streaming genotypes into per-contig temp files
// (<prefix>.contig<N>.json / .vcf) that are deleted on any exit. With --resume, low-mem/optimized streaming
// always genotypes through those per-contig temp files (at --threads 1 as well, with a single worker),
// flushes them after every locus, keeps them when the run fails, and maintains one more file beside them:
//
//   <prefix>.processed_loci.txt
//
// Its first line is "#" followed by the run signature (see buildRunSignature). Every other line is
// "<LocusId>\t<jsonFileSize>\t<vcfFileSize>": a finished locus, followed by the sizes in bytes of its
// contig's JSON and VCF temp files right after that locus's output was flushed to them. A line is appended
// only after that flush, so a killed process always leaves the listed bytes on disk. Loci that wrote nothing
// at all (dropped by --skip-hom-ref / --skip-missing-genotypes, or failed during genotyping) are listed too,
// so a resume does not genotype them again.
//
// On resume each contig's temp files are truncated to the sizes listed for that contig's last finished
// locus, which also drops a record the interrupted run was part way through writing, and the genotyping
// writers append from there. The final output comes from the same end-of-run merge as a run without
// --resume, so its bytes are unchanged.

#pragma once

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <mutex>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "core/Parameters.hh"
#include "core/ReferenceContigInfo.hh"
#include "locus/LocusSpecification.hh"

namespace ehunter
{

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

// <outputPrefix>.processed_loci.txt
std::string processedLociPath(const std::string& outputPrefix);

// Paths of the plain-text per-contig temp files the genotyping writers target.
std::string contigTempJsonPath(const std::string& outputPrefix, std::int32_t contigIndex);
std::string contigTempVcfPath(const std::string& outputPrefix, std::int32_t contigIndex);

// True when an interrupted run left a processed-loci list for this output prefix.
bool resumeCheckpointExists(const std::string& outputPrefix);

// A single-line JSON object holding every parameter that has to match for a checkpoint to be resumable:
// the sample, the analysis mode, the build, the catalog, and each flag that changes which loci are emitted
// or what is written for them. --threads and -z are deliberately excluded: the temp files are keyed by
// contig rather than by worker and are always uncompressed, so either may change between the two runs.
std::string buildRunSignature(const ProgramParameters& params);

// Throws when the catalog holds duplicate LocusIds, which resume cannot tell apart.
void assertCatalogIsResumable(const LocusDescriptionCatalog& catalog);

// Returns a human-readable description of the first field that differs between two signatures, or an
// empty string when they match.
std::string describeSignatureMismatch(const std::string& expected, const std::string& found);

// The run signature on the first line of an existing processed-loci list. Throws
// ResumeCheckpointCorruptError when the list has no readable signature line.
std::string readCheckpointRunSignature(const std::string& outputPrefix);

struct ResumedContig
{
    // Whether the contig's JSON temp holds at least one record, which decides whether the next record
    // appended to it needs a ", " separator.
    bool hasJsonRecords = false;
};

struct ResumeLoadResult
{
    // Loci that must not be genotyped again.
    std::unordered_set<std::string> doneLocusIds;
    // Contigs whose temp files hold finished loci, keyed by contig index. The caller opens these in append
    // mode and must include them in the end-of-run merge, even when the filtered catalog no longer has any
    // loci left on that contig.
    std::map<std::int32_t, ResumedContig> resumedContigs;
};

// Reads the processed-loci list, truncates each contig's temp files to the sizes listed for that contig's
// last finished locus, and rewrites the list to the loci it kept.
//
// A listed locus whose sizes are not all present in its contig's temp files lost writes to a crashing
// machine (nothing here forces data onto the disk, so a machine crash, unlike a killed process, can lose
// bytes the list already counts). It and every later locus on that contig are genotyped again.
//
// Throws ResumeCheckpointCorruptError when the list cannot be parsed, and plain std::runtime_error when it
// parses but does not belong to this run (a locus that is not in the catalog), which the caller must report
// rather than silently discard.
//
// `catalog` must be the whole catalog, before the finished loci are removed from it, and `contigInfo` the
// BAM/CRAM header's contigs: together they give the VCF header the temp files start with.
ResumeLoadResult loadResumeCheckpoint(
    const std::string& outputPrefix, const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo,
    const LocusDescriptionCatalog& catalog);

// Drops the loci that are already done, preserving the order of the rest.
void removeCompletedLoci(LocusDescriptionCatalog& catalog, const std::unordered_set<std::string>& doneLocusIds);

// Appends finished loci to the processed-loci list. Safe to call from several genotyping threads at once.
class ResumeCheckpointWriter
{
public:
    // `append` continues an existing (already trimmed) list instead of starting a new one headed by
    // `runSignature`, which must be a single line (as buildRunSignature's output is).
    ResumeCheckpointWriter(
        const std::string& outputPrefix, const std::string& runSignature, bool append, std::size_t abortAfterLoci);

    ResumeCheckpointWriter(const ResumeCheckpointWriter&) = delete;
    ResumeCheckpointWriter& operator=(const ResumeCheckpointWriter&) = delete;

    // Called once a locus's output has been flushed to its contig's temp files, with those files' sizes
    // (IterativeJsonWriter / IterativeVcfWriter::flushAndGetFileSize). The line is flushed before this
    // returns, and a failed write throws.
    void recordLocus(const std::string& locusId, std::uintmax_t jsonFileSize, std::uintmax_t vcfFileSize);

private:
    std::string path_;
    std::ofstream file_;
    std::mutex mutex_;
    std::size_t abortAfterLoci_;
    std::size_t recordedLocusCount_ = 0;
};

}
