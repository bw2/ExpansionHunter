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

#include "io/ResumeCheckpoint.hh"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <unordered_map>
#include <vector>

#include <boost/filesystem.hpp>

#include "spdlog/spdlog.h"

#include "app/Version.hh"
#include "core/Common.hh"
#include "genotype_quality/GenotypeQualityModel.hh"
#include "io/IterativeJsonWriter.hh"
#include "io/IterativeVcfWriter.hh"
#include "io/StringUtils.hh"

#include "thirdparty/json/json.hpp"

namespace fs = boost::filesystem;

namespace ehunter
{

using Json = nlohmann::json;

namespace
{

// One processed-loci entry: a finished locus and the sizes of its contig's JSON and VCF temp files right
// after its output was flushed to them.
struct ProcessedLocus
{
    std::string locusId;
    std::uintmax_t jsonFileSize = 0;
    std::uintmax_t vcfFileSize = 0;
};

struct ProcessedLociList
{
    std::string runSignature;
    std::vector<ProcessedLocus> loci;
};

std::string formatProcessedLocusLine(const std::string& locusId, std::uintmax_t jsonFileSize, std::uintmax_t vcfFileSize)
{
    return locusId + "\t" + std::to_string(jsonFileSize) + "\t" + std::to_string(vcfFileSize) + "\n";
}

// Reads the signature line, or throws ResumeCheckpointCorruptError. Leaves `file` positioned at the first
// locus line.
std::string readSignatureLine(std::ifstream& file, const std::string& path)
{
    std::string line;
    if (!std::getline(file, line) || file.eof() || line.empty() || line[0] != '#')
    {
        throw ResumeCheckpointCorruptError("no run signature on the first line of " + path);
    }
    return line.substr(1);
}

ProcessedLociList readProcessedLoci(const std::string& path)
{
    std::ifstream file(path, std::ios::in | std::ios::binary);
    if (!file)
    {
        throw ResumeCheckpointCorruptError("cannot open " + path);
    }

    ProcessedLociList list;
    list.runSignature = readSignatureLine(file, path);

    std::string line;
    while (std::getline(file, line))
    {
        if (file.eof())
        {
            // No newline after it: the interrupted run was killed part way through writing this line. Its
            // locus is simply not finished.
            break;
        }

        const std::size_t jsonSizeStart = line.find('\t');
        const std::size_t vcfSizeStart
            = jsonSizeStart == std::string::npos ? std::string::npos : line.find('\t', jsonSizeStart + 1);
        if (jsonSizeStart == 0 || vcfSizeStart == std::string::npos)
        {
            throw ResumeCheckpointCorruptError("malformed processed-loci entry in " + path + ": " + line);
        }

        ProcessedLocus locus;
        locus.locusId = line.substr(0, jsonSizeStart);
        try
        {
            std::size_t parsedLength = 0;
            const std::string jsonSize = line.substr(jsonSizeStart + 1, vcfSizeStart - jsonSizeStart - 1);
            locus.jsonFileSize = std::stoull(jsonSize, &parsedLength);
            if (parsedLength != jsonSize.size())
            {
                throw std::invalid_argument(jsonSize);
            }
            const std::string vcfSize = line.substr(vcfSizeStart + 1);
            locus.vcfFileSize = std::stoull(vcfSize, &parsedLength);
            if (parsedLength != vcfSize.size())
            {
                throw std::invalid_argument(vcfSize);
            }
        }
        catch (const std::exception&)
        {
            throw ResumeCheckpointCorruptError("malformed processed-loci entry in " + path + ": " + line);
        }
        list.loci.push_back(std::move(locus));
    }
    return list;
}

std::uintmax_t fileSizeOrZero(const std::string& path)
{
    boost::system::error_code error;
    const std::uintmax_t size = fs::file_size(path, error);
    return error ? 0 : size;
}

void renameOver(const std::string& from, const std::string& to)
{
    if (std::rename(from.c_str(), to.c_str()) != 0)
    {
        throw std::runtime_error(
            "Failed to replace " + to + " with " + from + " (" + std::strerror(errno) + ")");
    }
}

// The catalog's size, used to notice that it changed between an interrupted run and its resume.
// Deliberately not its modification time: a pipeline that localizes the catalog fresh for every attempt
// (as a cloud batch job does) would otherwise fail every resume over a byte-identical file.
std::string catalogStamp(const std::string& path)
{
    boost::system::error_code sizeError;
    const std::uintmax_t size = fs::file_size(path, sizeError);
    return sizeError ? "unknown" : std::to_string(size);
}

std::string sortCatalogByToString(SortCatalogBy sortCatalogBy)
{
    switch (sortCatalogBy)
    {
    case SortCatalogBy::kPosition: return "position";
    case SortCatalogBy::kLocusId: return "id";
    case SortCatalogBy::kNone: return "none";
    }
    return "unknown";
}

}  // namespace

std::string processedLociPath(const std::string& outputPrefix)
{
    return outputPrefix + ".processed_loci.txt";
}

std::string contigTempJsonPath(const std::string& outputPrefix, std::int32_t contigIndex)
{
    return outputPrefix + ".contig" + std::to_string(contigIndex) + ".json";
}

std::string contigTempVcfPath(const std::string& outputPrefix, std::int32_t contigIndex)
{
    return outputPrefix + ".contig" + std::to_string(contigIndex) + ".vcf";
}

bool resumeCheckpointExists(const std::string& outputPrefix) { return fs::exists(processedLociPath(outputPrefix)); }

void assertCatalogIsResumable(const LocusDescriptionCatalog& catalog)
{
    // Resume identifies a locus solely by its LocusId, so a catalog holding the same id twice cannot be
    // resumed: finishing one occurrence would mark both done and drop the other from the output. Such a
    // catalog is already malformed (the JSON output is keyed by LocusId, so the two entries collide into a
    // single record even without --resume), so refuse it up front rather than silently losing a locus.
    std::unordered_set<std::string> seenLocusIds;
    seenLocusIds.reserve(catalog.size());
    for (const LocusDescription& locusDescription : catalog)
    {
        if (!seenLocusIds.insert(locusDescription.locusId()).second)
        {
            throw std::runtime_error(
                "--resume cannot be used with this catalog: it contains more than one locus with the id '"
                + locusDescription.locusId() + "', and resume tells loci apart by their LocusId. Give each "
                "locus a unique LocusId, or run without --resume.");
        }
    }
}

std::string buildRunSignature(const ProgramParameters& params)
{
    Json signature;
    signature["Version"] = kCommitSha;
    signature["SampleId"] = params.sample().id();
    signature["Sex"] = streamToString(params.sample().sex());
    signature["AnalysisMode"] = analysisModeToString(params.analysisMode());
    signature["Reads"] = params.inputPaths().htsFile();
    signature["Reference"] = params.inputPaths().reference();
    signature["Catalog"] = params.inputPaths().catalog();
    signature["CatalogSize"] = catalogStamp(params.inputPaths().catalog());
    signature["SortCatalogBy"] = sortCatalogByToString(params.sortCatalogBy());
    signature["Locus"] = params.locus();
    signature["Region"] = params.region();
    signature["StartWith"] = params.startWith();
    signature["NLoci"] = params.nLoci();
    signature["SkipHomRef"] = params.skipHomRef();
    signature["SkipMissingGenotypes"] = params.skipMissingGenotypes();
    signature["HeuristicGenotypingOnly"] = params.heuristicGenotypingOnly();
    signature["MaxDepth"] = params.maxDepth();
    signature["QualityMetrics"] = params.enableAlleleQualityMetrics();
    signature["ConsensusSequences"] = params.enableConsensusSequences();
    signature["CopyCatalogFields"] = params.copyCatalogFields();
    signature["PlotAll"] = params.plotAll();
    signature["DisableAllPlots"] = params.disableAllPlots();
    signature["RegionExtensionLength"] = params.heuristics().regionExtensionLength();
    signature["MinLocusCoverage"] = params.heuristics().minLocusCoverage();
    signature["QualityCutoffForGoodBaseCall"] = params.heuristics().qualityCutoffForGoodBaseCall();
    signature["AlignerType"] = static_cast<int>(params.heuristics().alignerType());
    signature["OutputGenotypeTiming"] = params.outputGenotypeTiming();
    signature["MotifComposition"] = motifCompositionModeToString(params.motifCompositionMode());
    signature["GenotypeQualityModelVersion"]
        = params.genotypeQualityModel() ? params.genotypeQualityModel()->version : std::string();
    // No indentation, so the signature fits on the processed-loci list's first line.
    return signature.dump();
}

std::string describeSignatureMismatch(const std::string& expected, const std::string& found)
{
    Json expectedJson;
    Json foundJson;
    try
    {
        expectedJson = Json::parse(expected);
        foundJson = Json::parse(found);
    }
    catch (const std::exception&)
    {
        return "the checkpoint's run signature could not be parsed";
    }

    for (auto it = expectedJson.begin(); it != expectedJson.end(); ++it)
    {
        if (!foundJson.contains(it.key()))
        {
            return it.key() + " is missing from the checkpoint";
        }
        if (foundJson[it.key()] != it.value())
        {
            return it.key() + ": checkpoint has " + foundJson[it.key()].dump() + ", this run has "
                + it.value().dump();
        }
    }
    for (auto it = foundJson.begin(); it != foundJson.end(); ++it)
    {
        if (!expectedJson.contains(it.key()))
        {
            return it.key() + " is in the checkpoint but not in this run";
        }
    }
    return "";
}

std::string readCheckpointRunSignature(const std::string& outputPrefix)
{
    const std::string path = processedLociPath(outputPrefix);
    std::ifstream file(path, std::ios::in | std::ios::binary);
    if (!file)
    {
        throw ResumeCheckpointCorruptError("cannot open " + path);
    }
    return readSignatureLine(file, path);
}

ResumeLoadResult loadResumeCheckpoint(
    const std::string& outputPrefix, const SampleParameters& sampleParams, const ReferenceContigInfo& contigInfo,
    const LocusDescriptionCatalog& catalog)
{
    const std::string listPath = processedLociPath(outputPrefix);
    const ProcessedLociList list = readProcessedLoci(listPath);

    std::unordered_map<std::string, std::int32_t> contigOfLocus;
    contigOfLocus.reserve(catalog.size());
    for (const LocusDescription& locusDescription : catalog)
    {
        contigOfLocus[locusDescription.locusId()] = locusDescription.locusContigIndex();
    }

    // Every listed locus has to be in the catalog: without its contig there is no temp file to find its
    // output in, and it would silently vanish from the final output.
    for (const ProcessedLocus& locus : list.loci)
    {
        if (contigOfLocus.find(locus.locusId) == contigOfLocus.end())
        {
            throw std::runtime_error(
                "checkpoint locus " + locus.locusId + " is not in the catalog, so the checkpoint does not "
                "belong to this run");
        }
    }

    // A contig's first record is appended right after the header its writer starts the file with, so the
    // listed sizes can never be smaller than the headers.
    const std::uintmax_t jsonHeaderSize = jsonDocumentHeader(sampleParams).size();
    const std::uintmax_t vcfHeaderSize
        = vcfDocumentHeader(sampleParams.id(), contigInfo, contigsWithLoci(catalog)).size();

    struct ContigProgress
    {
        std::uintmax_t jsonSizeOnDisk = 0;
        std::uintmax_t vcfSizeOnDisk = 0;
        // Sizes listed for the last locus kept so far (the headers until one is kept).
        std::uintmax_t jsonFileSize = 0;
        std::uintmax_t vcfFileSize = 0;
        bool anyLocusKept = false;
        bool stillKeeping = true;
    };
    std::map<std::int32_t, ContigProgress> contigs;

    std::vector<const ProcessedLocus*> keptLoci;
    keptLoci.reserve(list.loci.size());
    for (const ProcessedLocus& locus : list.loci)
    {
        const std::int32_t contigIndex = contigOfLocus.at(locus.locusId);
        auto contigIt = contigs.find(contigIndex);
        if (contigIt == contigs.end())
        {
            ContigProgress progress;
            progress.jsonSizeOnDisk = fileSizeOrZero(contigTempJsonPath(outputPrefix, contigIndex));
            progress.vcfSizeOnDisk = fileSizeOrZero(contigTempVcfPath(outputPrefix, contigIndex));
            progress.jsonFileSize = jsonHeaderSize;
            progress.vcfFileSize = vcfHeaderSize;
            contigIt = contigs.emplace(contigIndex, progress).first;
        }
        ContigProgress& progress = contigIt->second;

        // A contig's files only ever grow, so its listed sizes never shrink from one locus to the next.
        // Once one locus's sizes are not all on disk, that locus and every later one on its contig are
        // redone: its writer can only append after the last locus that is.
        progress.stillKeeping = progress.stillKeeping && locus.jsonFileSize >= progress.jsonFileSize
            && locus.vcfFileSize >= progress.vcfFileSize && locus.jsonFileSize <= progress.jsonSizeOnDisk
            && locus.vcfFileSize <= progress.vcfSizeOnDisk;
        if (!progress.stillKeeping)
        {
            continue;
        }
        progress.jsonFileSize = locus.jsonFileSize;
        progress.vcfFileSize = locus.vcfFileSize;
        progress.anyLocusKept = true;
        keptLoci.push_back(&locus);
    }

    if (keptLoci.size() < list.loci.size())
    {
        spdlog::warn(
            "Resume: {} of the {} loci listed in {} are missing output from their contig's temp files and "
            "will be genotyped again",
            add_commas_at_thousands(list.loci.size() - keptLoci.size()), add_commas_at_thousands(list.loci.size()),
            listPath);
    }

    // Rewrite the list to exactly the kept loci, via a sibling file and a rename: it must not go on naming a
    // locus whose output is about to be cut off, and a partly written last line has to be gone before the
    // resumed run appends after it.
    const std::string trimmedListPath = listPath + ".trimmed";
    {
        std::ofstream trimmed(trimmedListPath, std::ios::out | std::ios::binary | std::ios::trunc);
        if (!trimmed)
        {
            throw std::runtime_error("Failed to open file for writing: " + trimmedListPath);
        }
        trimmed << "#" << list.runSignature << "\n";
        for (const ProcessedLocus* locus : keptLoci)
        {
            trimmed << formatProcessedLocusLine(locus->locusId, locus->jsonFileSize, locus->vcfFileSize);
        }
        trimmed.flush();
        if (!trimmed)
        {
            throw std::runtime_error("Failed to write " + trimmedListPath + " (" + std::strerror(errno) + ")");
        }
    }
    renameOver(trimmedListPath, listPath);

    ResumeLoadResult result;
    for (const ProcessedLocus* locus : keptLoci)
    {
        result.doneLocusIds.insert(locus->locusId);
    }
    for (const auto& contigAndProgress : contigs)
    {
        const std::int32_t contigIndex = contigAndProgress.first;
        const ContigProgress& progress = contigAndProgress.second;
        if (!progress.anyLocusKept)
        {
            continue;  // every locus on this contig is genotyped again, into a freshly truncated file
        }
        // Cuts off whatever followed the last kept locus: a record the interrupted run was part way through
        // writing, the records of loci dropped above, or the closing footer a writer adds when it unwinds.
        fs::resize_file(contigTempJsonPath(outputPrefix, contigIndex), progress.jsonFileSize);
        fs::resize_file(contigTempVcfPath(outputPrefix, contigIndex), progress.vcfFileSize);
        result.resumedContigs[contigIndex] = ResumedContig{ progress.jsonFileSize > jsonHeaderSize };
    }
    return result;
}

void removeCompletedLoci(LocusDescriptionCatalog& catalog, const std::unordered_set<std::string>& doneLocusIds)
{
    if (doneLocusIds.empty())
    {
        return;
    }
    catalog.erase(
        std::remove_if(
            catalog.begin(), catalog.end(),
            [&doneLocusIds](const LocusDescription& locusDescription) {
                return doneLocusIds.find(locusDescription.locusId()) != doneLocusIds.end();
            }),
        catalog.end());
}

ResumeCheckpointWriter::ResumeCheckpointWriter(
    const std::string& outputPrefix, const std::string& runSignature, bool append, std::size_t abortAfterLoci)
    : path_(processedLociPath(outputPrefix))
    , abortAfterLoci_(abortAfterLoci)
{
    if (runSignature.find('\n') != std::string::npos)
    {
        throw std::logic_error("the run signature must be a single line");
    }

    // A fresh list discards anything left by an earlier run that this one is not resuming.
    file_.open(path_, std::ios::out | std::ios::binary | (append ? std::ios::app : std::ios::trunc));
    if (!file_)
    {
        throw std::runtime_error("Failed to open checkpoint file for writing: " + path_);
    }
    if (!append)
    {
        file_ << "#" << runSignature << "\n";
        file_.flush();
        if (!file_)
        {
            throw std::runtime_error("Failed to write " + path_ + " (" + std::strerror(errno) + ")");
        }
    }
}

void ResumeCheckpointWriter::recordLocus(
    const std::string& locusId, std::uintmax_t jsonFileSize, std::uintmax_t vcfFileSize)
{
    std::lock_guard<std::mutex> lock(mutex_);
    file_ << formatProcessedLocusLine(locusId, jsonFileSize, vcfFileSize);
    file_.flush();
    if (!file_)
    {
        throw std::runtime_error("Failed to write " + path_ + " (" + std::strerror(errno) + ")");
    }

    ++recordedLocusCount_;
    if (abortAfterLoci_ != 0 && recordedLocusCount_ >= abortAfterLoci_)
    {
        // Test hook (--internal-abort-after-loci): leave the process exactly as an external kill would,
        // with the checkpoint on disk and nothing unwound.
        spdlog::warn("Aborting after {} checkpointed loci (--internal-abort-after-loci)", recordedLocusCount_);
        std::_Exit(137);
    }
}

}
