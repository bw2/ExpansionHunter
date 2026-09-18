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
#include <cctype>
#include <limits>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <memory>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include <zlib.h>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

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

const char* const kCheckpointSuffix = ".unfinished";
const char* const kLocusResultsMarker = "\"LocusResults\": {";
const char* const kResumeInfoMarker = "\"ResumeInfo\": ";

bool endsWithGz(const std::string& path)
{
    return path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

// True for a checkpoint path whose payload is gzip-compressed, i.e. one derived from a compressed output
// path (<prefix>.json.gz.unfinished).
bool checkpointIsCompressed(const std::string& path)
{
    const std::string suffix(kCheckpointSuffix);
    if (path.size() > suffix.size() && path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0)
    {
        return endsWithGz(path.substr(0, path.size() - suffix.size()));
    }
    return endsWithGz(path);
}

// Reads a possibly-truncated (and possibly gzip-compressed) file, handing out every byte that could be
// decoded before the damage.
//
// Decompression goes through zlib directly rather than boost's gzip filter because of what happens at a
// truncated tail. boost's filter throws from the stream buffer, and std::istream::read then reports
// nothing at all for that call -- not even the bytes it had already decoded -- so a checkpoint whose last
// gzip member is cut short would lose everything decoded alongside it, which for a small checkpoint is the
// entire run's progress. zlib's inflate() instead reports exactly how much output it produced before it
// ran out of input, so nothing recoverable is thrown away.
//
// Multi-member input is expected: the checkpoint writer closes a gzip member on every flush, so the file
// is a chain of them. Each Z_STREAM_END resets the stream for the next member.
class TolerantInput
{
public:
    TolerantInput(const std::string& path, bool compressed)
        : path_(path)
        , compressed_(compressed)
        , inBuffer_(kInBufferSize)
    {
        file_.open(path, std::ios::in | std::ios::binary);
        if (!file_)
        {
            throw std::runtime_error("Failed to open file: " + path);
        }
        if (compressed_ && !startMember())
        {
            finished_ = true;
        }
    }

    ~TolerantInput()
    {
        if (streamActive_)
        {
            inflateEnd(&stream_);
        }
    }

    TolerantInput(const TolerantInput&) = delete;
    TolerantInput& operator=(const TolerantInput&) = delete;

    // Returns the number of bytes decoded into `buffer`; 0 means end of stream, or a tail that could not
    // be decoded at all.
    std::size_t read(char* buffer, std::size_t count)
    {
        if (!compressed_)
        {
            file_.read(buffer, static_cast<std::streamsize>(count));
            return static_cast<std::size_t>(file_.gcount());
        }

        if (finished_ || count == 0)
        {
            return 0;
        }

        stream_.next_out = reinterpret_cast<Bytef*>(buffer);
        stream_.avail_out = static_cast<uInt>(count);

        while (stream_.avail_out > 0)
        {
            if (stream_.avail_in == 0 && !refill())
            {
                // Out of input. Anything already inflated into the caller's buffer still counts; the
                // records it holds are kept and whatever was cut short is simply never reported.
                finished_ = true;
                break;
            }

            const int status = inflate(&stream_, Z_NO_FLUSH);
            if (status == Z_STREAM_END)
            {
                // End of one gzip member, its CRC verified by zlib. Another may follow; if what comes next
                // is not a member header, the inflate() after this reset reports a data error.
                if (inflateReset(&stream_) != Z_OK)
                {
                    finished_ = true;
                    break;
                }
                continue;
            }
            if (status != Z_OK)
            {
                // A zlib data error means the bytes themselves are damaged, not merely cut short: whatever
                // inflate() produced before noticing is not trustworthy, and unlike a truncated tail it
                // cannot be salvaged by keeping the complete members (the damage may sit anywhere in this
                // one). Treat the checkpoint as unreadable so the caller starts over rather than replaying
                // corrupted records as if they were real results.
                throw ResumeCheckpointCorruptError(
                    "corrupted gzip data in " + path_ + " (" + (stream_.msg ? stream_.msg : "inflate error")
                    + ")");
            }
        }

        return count - static_cast<std::size_t>(stream_.avail_out);
    }

private:
    static const std::size_t kInBufferSize = 64 * 1024;

    bool startMember()
    {
        stream_ = z_stream{};
        // 15 window bits plus 16 selects gzip framing rather than raw zlib.
        if (inflateInit2(&stream_, 15 + 16) != Z_OK)
        {
            return false;
        }
        streamActive_ = true;
        return true;
    }

    bool refill()
    {
        file_.read(inBuffer_.data(), static_cast<std::streamsize>(inBuffer_.size()));
        const std::streamsize got = file_.gcount();
        if (got <= 0)
        {
            return false;
        }
        stream_.next_in = reinterpret_cast<Bytef*>(inBuffer_.data());
        stream_.avail_in = static_cast<uInt>(got);
        return true;
    }

    std::string path_;
    std::ifstream file_;
    bool compressed_;
    std::vector<char> inBuffer_;
    z_stream stream_{};
    bool streamActive_ = false;
    bool finished_ = false;
};

// Writes a file, compressing it iff `compressed`. Used to rewrite trimmed checkpoints in one pass.
class OutputFile
{
public:
    OutputFile(const std::string& path, bool compressed)
        : path_(path)
    {
        file_.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
        if (!file_)
        {
            throw std::runtime_error("Failed to open file for writing: " + path);
        }
        if (compressed)
        {
            stream_.push(boost::iostreams::gzip_compressor());
        }
        stream_.push(file_);
    }

    ~OutputFile()
    {
        try
        {
            close();
        }
        catch (...)
        {
            // Never throw from a destructor.
        }
    }

    void write(const std::string& text) { stream_ << text; }

    void close()
    {
        if (closed_)
        {
            return;
        }
        closed_ = true;
        stream_.flush();
        stream_.reset();
        file_.flush();
        const bool failed = !file_;
        if (file_.is_open())
        {
            file_.close();
        }
        // A silently truncated rewrite would look like a checkpoint whose tail was lost, quietly dropping
        // the loci it could not write; fail the run instead.
        if (failed || !file_)
        {
            throw std::runtime_error("Failed to write " + path_ + " (" + std::strerror(errno) + ")");
        }
    }

private:
    std::string path_;
    std::ofstream file_;
    boost::iostreams::filtering_ostream stream_;
    bool closed_ = false;
};

// Appends `text` to `path` as one self-contained gzip member (or as raw bytes when not compressing), then
// flushes it to the OS. Concatenated gzip members form a valid gzip file, so the checkpoint is complete and
// readable on disk after every call.
void appendToFile(const std::string& path, const std::string& text, bool compressed)
{
    if (text.empty())
    {
        return;
    }

    std::ofstream file(path, std::ios::out | std::ios::binary | std::ios::app);
    if (!file)
    {
        throw std::runtime_error("Failed to open checkpoint file for appending: " + path);
    }

    {
        boost::iostreams::filtering_ostream stream;
        if (compressed)
        {
            stream.push(boost::iostreams::gzip_compressor());
        }
        stream.push(file);
        stream << text;
        stream.flush();
        stream.reset();
    }

    file.flush();
    if (!file)
    {
        throw std::runtime_error("Failed to append to checkpoint file: " + path);
    }
    file.close();
}

// Index just past the JSON string starting at `start` (which must be the opening quote), or npos when the
// buffer ends before the closing quote. `value` receives the unescaped content.
std::size_t parseJsonString(const std::string& buffer, std::size_t start, std::string& value)
{
    if (start >= buffer.size() || buffer[start] != '"')
    {
        return std::string::npos;
    }

    value.clear();
    for (std::size_t i = start + 1; i < buffer.size(); ++i)
    {
        const char c = buffer[i];
        if (c == '\\')
        {
            if (i + 1 >= buffer.size())
            {
                return std::string::npos;
            }
            const char escaped = buffer[i + 1];
            switch (escaped)
            {
            case 'n': value += '\n'; break;
            case 't': value += '\t'; break;
            case 'r': value += '\r'; break;
            case 'b': value += '\b'; break;
            case 'f': value += '\f'; break;
            case 'u':
                // A \uXXXX escape cannot appear in a LocusId that ExpansionHunter itself wrote, so keeping
                // the escape verbatim is enough: such an id would simply not match any catalog id.
                if (i + 5 >= buffer.size())
                {
                    return std::string::npos;
                }
                value.append(buffer, i, 6);
                i += 4;
                break;
            default: value += escaped; break;
            }
            ++i;
        }
        else if (c == '"')
        {
            return i + 1;
        }
        else
        {
            value += c;
        }
    }
    return std::string::npos;
}

// Index just past the JSON object starting at `start` (which must be '{'), or npos when the buffer ends
// before the object closes.
std::size_t findObjectEnd(const std::string& buffer, std::size_t start)
{
    if (start >= buffer.size() || buffer[start] != '{')
    {
        return std::string::npos;
    }

    int depth = 0;
    bool inString = false;
    bool escaped = false;
    for (std::size_t i = start; i < buffer.size(); ++i)
    {
        const char c = buffer[i];
        if (escaped)
        {
            escaped = false;
            continue;
        }
        if (inString)
        {
            if (c == '\\')
            {
                escaped = true;
            }
            else if (c == '"')
            {
                inString = false;
            }
            continue;
        }
        if (c == '"')
        {
            inString = true;
        }
        else if (c == '{')
        {
            ++depth;
        }
        else if (c == '}')
        {
            if (--depth == 0)
            {
                return i + 1;
            }
        }
    }
    return std::string::npos;
}

// Streams the LocusResults entries of an ExpansionHunter JSON document, tolerating a truncated tail: a
// record that is cut off part-way through is simply not reported. `onRecord` receives the locus id and the
// exact record bytes (without the ", " separator), so they can be replayed verbatim.
// Returns false when the document header was never found, which means the file is not a usable checkpoint.
bool scanJsonRecords(
    const std::string& path, bool compressed,
    const std::function<void(const std::string&, const std::string&)>& onRecord)
{
    const std::size_t kChunkSize = 64 * 1024;
    TolerantInput input(path, compressed);

    std::string buffer;
    std::vector<char> chunk(kChunkSize);
    bool foundHeader = false;
    std::size_t pos = 0;
    bool atEnd = false;

    while (true)
    {
        const std::size_t bytesRead = input.read(chunk.data(), kChunkSize);
        if (bytesRead == 0)
        {
            atEnd = true;
        }
        else
        {
            buffer.append(chunk.data(), bytesRead);
        }

        if (!foundHeader)
        {
            const std::size_t markerPos = buffer.find(kLocusResultsMarker);
            if (markerPos == std::string::npos)
            {
                if (atEnd)
                {
                    return false;
                }
                continue;
            }
            foundHeader = true;
            pos = markerPos + std::strlen(kLocusResultsMarker);
        }

        // Parse as many complete records as the buffer holds, then drop what has been consumed so the
        // buffer stays the size of a record rather than of the file.
        bool documentEnded = false;
        while (true)
        {
            std::size_t cursor = pos;
            while (cursor < buffer.size() && (buffer[cursor] == ',' || buffer[cursor] == ' '))
            {
                ++cursor;
            }
            const std::size_t recordStart = cursor;

            while (cursor < buffer.size() && std::isspace(static_cast<unsigned char>(buffer[cursor])))
            {
                ++cursor;
            }
            if (cursor >= buffer.size())
            {
                break;  // need more data
            }
            if (buffer[cursor] != '"')
            {
                // Either the closing '}' of LocusResults, or trailing garbage from a truncated write.
                documentEnded = true;
                break;
            }

            std::string locusId;
            const std::size_t afterKey = parseJsonString(buffer, cursor, locusId);
            if (afterKey == std::string::npos)
            {
                break;  // need more data
            }

            std::size_t valueStart = afterKey;
            while (valueStart < buffer.size()
                && (std::isspace(static_cast<unsigned char>(buffer[valueStart])) || buffer[valueStart] == ':'))
            {
                ++valueStart;
            }
            if (valueStart >= buffer.size())
            {
                break;  // need more data
            }
            if (buffer[valueStart] != '{')
            {
                documentEnded = true;  // malformed tail
                break;
            }

            const std::size_t recordEnd = findObjectEnd(buffer, valueStart);
            if (recordEnd == std::string::npos)
            {
                break;  // need more data
            }

            onRecord(locusId, buffer.substr(recordStart, recordEnd - recordStart));
            pos = recordEnd;
        }

        if (documentEnded)
        {
            return true;
        }

        if (pos > 0)
        {
            buffer.erase(0, pos);
            pos = 0;
        }

        if (atEnd)
        {
            return true;
        }
    }
}

// Streams the complete lines of a possibly-truncated text file. A final line without a newline is a
// partial write and is skipped. Each line is passed with its trailing newline so it can be replayed
// verbatim.
void scanLines(const std::string& path, bool compressed, const std::function<void(const std::string&)>& onLine)
{
    const std::size_t kChunkSize = 64 * 1024;
    TolerantInput input(path, compressed);

    std::string buffer;
    std::vector<char> chunk(kChunkSize);
    while (true)
    {
        const std::size_t bytesRead = input.read(chunk.data(), kChunkSize);
        if (bytesRead > 0)
        {
            buffer.append(chunk.data(), bytesRead);
        }

        std::size_t lineStart = 0;
        std::size_t newlinePos = buffer.find('\n', lineStart);
        while (newlinePos != std::string::npos)
        {
            onLine(buffer.substr(lineStart, newlinePos - lineStart + 1));
            lineStart = newlinePos + 1;
            newlinePos = buffer.find('\n', lineStart);
        }
        buffer.erase(0, lineStart);

        if (bytesRead == 0)
        {
            return;
        }
    }
}

// The VARID value from a VCF record line's INFO column, or an empty string when the line has no INFO
// column or no VARID key.
std::string variantIdOfVcfLine(const std::string& line)
{
    const int kInfoColumnIndex = 7;
    std::size_t fieldStart = 0;
    for (int column = 0; column < kInfoColumnIndex; ++column)
    {
        const std::size_t tabPos = line.find('\t', fieldStart);
        if (tabPos == std::string::npos)
        {
            return "";
        }
        fieldStart = tabPos + 1;
    }

    std::size_t fieldEnd = line.find('\t', fieldStart);
    if (fieldEnd == std::string::npos)
    {
        fieldEnd = line.size();
    }

    const std::string info = line.substr(fieldStart, fieldEnd - fieldStart);
    const std::string key = "VARID=";
    std::size_t keyPos = info.find(key);
    while (keyPos != std::string::npos && keyPos != 0 && info[keyPos - 1] != ';')
    {
        keyPos = info.find(key, keyPos + 1);  // not a key boundary, e.g. "REPID=..." would not match anyway
    }
    if (keyPos == std::string::npos)
    {
        return "";
    }

    const std::size_t valueStart = keyPos + key.size();
    const std::size_t valueEnd = info.find(';', valueStart);
    return info.substr(valueStart, valueEnd == std::string::npos ? std::string::npos : valueEnd - valueStart);
}

std::string checkpointJsonHeader(const SampleParameters& sampleParams, const std::string& runSignature)
{
    const std::string indentedSignature = std::regex_replace(runSignature, std::regex("\n"), "\n  ");
    const std::string base = jsonDocumentHeader(sampleParams);

    // jsonDocumentHeader ends with ",\n  \"LocusResults\": {"; the ResumeInfo object goes just before it,
    // reusing the comma that already separates it from SampleParameters.
    const std::string marker = std::string(",\n  ") + kLocusResultsMarker;
    const std::size_t markerPos = base.rfind(marker);
    if (markerPos == std::string::npos)
    {
        throw std::logic_error("jsonDocumentHeader no longer ends with the LocusResults marker");
    }
    return base.substr(0, markerPos) + ",\n  \"ResumeInfo\": " + indentedSignature + marker;
}

// The ResumeInfo object stored in a checkpoint's header, as JSON text.
std::string readRunSignatureFromFile(const std::string& jsonCheckpointPath, bool compressed)
{
    // Read the header in growing chunks rather than one fixed window: the signature quotes options
    // verbatim, and a long --locus list makes it far bigger than any window picked up front. Getting this
    // wrong would not be subtle, since failing to read the signature discards the whole checkpoint.
    const std::size_t kReadChunk = 64 * 1024;
    const std::size_t kMaxHeaderSize = 64u * 1024 * 1024;

    TolerantInput input(jsonCheckpointPath, compressed);
    std::string header;
    std::vector<char> chunk(kReadChunk);
    std::size_t objectStart = std::string::npos;
    std::size_t objectEnd = std::string::npos;

    while (true)
    {
        const std::size_t bytesRead = input.read(chunk.data(), kReadChunk);
        if (bytesRead > 0)
        {
            header.append(chunk.data(), bytesRead);
        }

        const std::size_t markerPos = header.find(kResumeInfoMarker);
        if (markerPos != std::string::npos)
        {
            objectStart = markerPos + std::strlen(kResumeInfoMarker);
            objectEnd = findObjectEnd(header, objectStart);
            if (objectEnd != std::string::npos)
            {
                break;
            }
        }

        if (bytesRead == 0)
        {
            throw ResumeCheckpointCorruptError(
                (markerPos == std::string::npos ? "no ResumeInfo record in " : "truncated ResumeInfo record in ")
                + jsonCheckpointPath);
        }
        if (header.size() >= kMaxHeaderSize)
        {
            throw ResumeCheckpointCorruptError("ResumeInfo record is implausibly large in " + jsonCheckpointPath);
        }
    }

    // Re-serialize rather than returning the stored text: checkpointJsonHeader indents whatever it is
    // given, and this same text is handed back to it when the checkpoint is rewritten, so returning it
    // as-is would add two more spaces to every line on each successive resume.
    try
    {
        return Json::parse(header.substr(objectStart, objectEnd - objectStart)).dump(2);
    }
    catch (const std::exception&)
    {
        throw ResumeCheckpointCorruptError("unparseable ResumeInfo record in " + jsonCheckpointPath);
    }
}

// One processed_loci entry: a finished locus and how many records it wrote to each file, as
// "<LocusId>\t<jsonRecords>\t<vcfLines>". The counts let a resume tell "this locus produced no record"
// apart from "this locus's records are missing from the file", so a checkpoint whose record files lost
// buffered writes is detected rather than silently dropping those loci from the output. The VCF count is
// per line rather than a flag because a multi-variant locus writes one line per variant, and losing just
// one of them has to be caught too.
struct ProcessedLocus
{
    std::string locusId;
    std::size_t jsonRecordCount = 0;
    std::size_t vcfLineCount = 0;
};

std::string formatProcessedLocusLine(
    const std::string& locusId, std::size_t jsonRecordCount, std::size_t vcfLineCount)
{
    return locusId + "\t" + std::to_string(jsonRecordCount) + "\t" + std::to_string(vcfLineCount) + "\n";
}

std::size_t countLines(const std::string& text)
{
    return static_cast<std::size_t>(std::count(text.begin(), text.end(), '\n'));
}

std::vector<ProcessedLocus> readProcessedLoci(const std::string& path)
{
    std::vector<ProcessedLocus> processedLoci;
    scanLines(path, false, [&processedLoci](const std::string& line) {
        std::string trimmed = line;
        while (!trimmed.empty() && (trimmed.back() == '\n' || trimmed.back() == '\r'))
        {
            trimmed.pop_back();
        }
        if (trimmed.empty())
        {
            return;
        }

        const std::size_t jsonCountStart = trimmed.find('\t');
        const std::size_t vcfCountStart
            = jsonCountStart == std::string::npos ? std::string::npos : trimmed.find('\t', jsonCountStart + 1);
        if (vcfCountStart == std::string::npos)
        {
            throw ResumeCheckpointCorruptError("malformed processed-loci entry: " + trimmed);
        }

        ProcessedLocus processedLocus;
        processedLocus.locusId = trimmed.substr(0, jsonCountStart);
        try
        {
            processedLocus.jsonRecordCount = std::stoul(
                trimmed.substr(jsonCountStart + 1, vcfCountStart - jsonCountStart - 1));
            processedLocus.vcfLineCount = std::stoul(trimmed.substr(vcfCountStart + 1));
        }
        catch (const std::exception&)
        {
            throw ResumeCheckpointCorruptError("malformed processed-loci entry: " + trimmed);
        }
        processedLoci.push_back(std::move(processedLocus));
    });
    return processedLoci;
}

// Index of the first processed-loci entry that cannot be reused, given how this run slices its output.
//
// Each slice's rebuilt temp is appended to by a genotyping writer, which can only add records at the end,
// so what a slice needs is that the finished loci are a *prefix of the order that slice emits records in*.
// That is not catalog order: a locus is released for genotyping once the reads pass the end of its flank
// window, so a locus nested inside a wider one finishes first even though it comes second in the
// position-sorted catalog. The checkpoint records the order the interrupted run emitted in, so within one
// slice any leading run of it is by construction a valid prefix and nothing needs checking.
//
// What does need checking is the slice layout changing between the two runs. At --threads > 1 a slice is
// one contig, and a contig's entries keep their relative order no matter which run wrote them, so every
// entry is reusable and this returns the full count. At --threads 1 a single slice covers every contig in
// one coordinate sweep, so its rebuilt temp has to be a prefix of that sweep: contigs may only move
// forwards, and a contig may only be left behind once every one of its loci is finished. A --threads > 1
// checkpoint generally fails that (it can finish chr3 while chr1 is still running), and then this returns
// a count smaller than the list; loadResumeCheckpoint treats that as a refusal and stops before rewriting
// any file, rather than trimming away work the operator can still use by re-running with --threads > 1.
std::size_t firstUnusableProcessedLocus(
    const std::vector<ProcessedLocus>& processedLoci, const LocusDescriptionCatalog& catalog,
    const std::unordered_map<std::string, std::int32_t>& contigOfLocus, bool perContig)
{
    if (perContig)
    {
        return processedLoci.size();
    }

    std::map<std::int32_t, std::size_t> lociPerContig;
    for (const LocusDescription& locusDescription : catalog)
    {
        ++lociPerContig[locusDescription.locusContigIndex()];
    }

    std::vector<std::int32_t> contigsWithLoci;
    contigsWithLoci.reserve(lociPerContig.size());
    for (const auto& contigAndCount : lociPerContig)
    {
        contigsWithLoci.push_back(contigAndCount.first);
    }

    std::map<std::int32_t, std::size_t> finishedPerContig;
    // Index into contigsWithLoci of the first contig that is not yet fully finished; everything before it
    // is complete, which is what lets the sweep move past those contigs.
    std::size_t firstUnfinishedContig = 0;
    auto advanceFirstUnfinishedContig = [&]() {
        while (firstUnfinishedContig < contigsWithLoci.size()
            && finishedPerContig[contigsWithLoci[firstUnfinishedContig]]
                == lociPerContig[contigsWithLoci[firstUnfinishedContig]])
        {
            ++firstUnfinishedContig;
        }
    };
    advanceFirstUnfinishedContig();

    std::size_t usableCount = 0;
    for (const ProcessedLocus& processedLocus : processedLoci)
    {
        const std::int32_t contigIndex = contigOfLocus.at(processedLocus.locusId);
        const std::size_t contigPosition = static_cast<std::size_t>(std::distance(contigsWithLoci.begin(),
            std::lower_bound(contigsWithLoci.begin(), contigsWithLoci.end(), contigIndex)));
        // The sweep can only be working on the earliest contig that still has unfinished loci.
        if (contigPosition != firstUnfinishedContig)
        {
            break;
        }
        ++finishedPerContig[contigIndex];
        advanceFirstUnfinishedContig();
        ++usableCount;
    }
    return usableCount;
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

// Per-slice temp files rebuilt from a checkpoint. Writing them directly (rather than through
// IterativeJsonWriter / IterativeVcfWriter) keeps them header-only-plus-records: the genotyping writers
// append to them afterwards and close them, which is what adds the JSON footer.
class SliceRebuilder
{
public:
    SliceRebuilder(std::string outputPrefix, const SampleParameters& sampleParams, bool perContig)
        : outputPrefix_(std::move(outputPrefix))
        , sampleParams_(sampleParams)
        , perContig_(perContig)
    {
    }

    void addJsonRecord(std::int32_t contigIndex, const std::string& recordText)
    {
        Slice& slice = sliceFor(contigIndex);
        if (slice.hasJsonRecords)
        {
            *slice.jsonFile << ", ";
        }
        *slice.jsonFile << recordText;
        slice.hasJsonRecords = true;
    }

    void addVcfLine(std::int32_t contigIndex, const std::string& line) { *sliceFor(contigIndex).vcfFile << line; }

    std::map<std::int32_t, RebuiltSlice> finish()
    {
        std::map<std::int32_t, RebuiltSlice> result;
        for (auto& keyAndSlice : slices_)
        {
            for (std::ofstream* file : { keyAndSlice.second.jsonFile.get(), keyAndSlice.second.vcfFile.get() })
            {
                file->flush();
                const bool failed = !*file;
                file->close();
                // A temp file silently truncated here would drop its loci from the merged output while the
                // checkpoint still counts them as finished, so fail the run instead.
                if (failed || !*file)
                {
                    throw std::runtime_error(
                        "Failed to write a resume temp file for slice " + std::to_string(keyAndSlice.first) + " ("
                        + std::strerror(errno) + ")");
                }
            }
            result[keyAndSlice.first] = RebuiltSlice{ keyAndSlice.second.hasJsonRecords };
        }
        slices_.clear();
        return result;
    }

private:
    struct Slice
    {
        std::unique_ptr<std::ofstream> jsonFile;
        std::unique_ptr<std::ofstream> vcfFile;
        bool hasJsonRecords = false;
    };

    Slice& sliceFor(std::int32_t contigIndex)
    {
        const std::int32_t key = perContig_ ? contigIndex : kSingleSliceKey;
        auto it = slices_.find(key);
        if (it != slices_.end())
        {
            return it->second;
        }

        const std::string jsonPath
            = perContig_ ? contigTempJsonPath(outputPrefix_, key) : singleSliceTempJsonPath(outputPrefix_);
        const std::string vcfPath
            = perContig_ ? contigTempVcfPath(outputPrefix_, key) : singleSliceTempVcfPath(outputPrefix_);

        Slice slice;
        slice.jsonFile = std::make_unique<std::ofstream>(jsonPath, std::ios::out | std::ios::binary | std::ios::trunc);
        slice.vcfFile = std::make_unique<std::ofstream>(vcfPath, std::ios::out | std::ios::binary | std::ios::trunc);
        if (!*slice.jsonFile || !*slice.vcfFile)
        {
            throw std::runtime_error("Failed to open temp file for writing: " + jsonPath + " / " + vcfPath);
        }
        *slice.jsonFile << jsonDocumentHeader(sampleParams_);
        *slice.vcfFile << vcfDocumentHeader(sampleParams_.id());

        return slices_.emplace(key, std::move(slice)).first->second;
    }

    std::string outputPrefix_;
    const SampleParameters& sampleParams_;
    bool perContig_;
    std::map<std::int32_t, Slice> slices_;
};

}  // namespace

ResumeCheckpointPaths resumeCheckpointPaths(const OutputPaths& outputPaths)
{
    return ResumeCheckpointPaths{ outputPaths.json() + kCheckpointSuffix, outputPaths.vcf() + kCheckpointSuffix,
        outputPaths.outputPrefix() + ".processed_loci" + kCheckpointSuffix };
}

OutputPaths otherCompressionOutputPaths(const OutputPaths& outputPaths)
{
    auto flipGz = [](const std::string& path) {
        return endsWithGz(path) ? path.substr(0, path.size() - 3) : path + ".gz";
    };
    return OutputPaths(flipGz(outputPaths.vcf()), flipGz(outputPaths.json()), outputPaths.bamlet(),
        outputPaths.timing(), outputPaths.outputPrefix());
}

std::string contigTempJsonPath(const std::string& outputPrefix, std::int32_t contigIndex)
{
    return outputPrefix + ".contig" + std::to_string(contigIndex) + ".json";
}

std::string contigTempVcfPath(const std::string& outputPrefix, std::int32_t contigIndex)
{
    return outputPrefix + ".contig" + std::to_string(contigIndex) + ".vcf";
}

std::string singleSliceTempJsonPath(const std::string& outputPrefix) { return outputPrefix + ".resume_part.json"; }

std::string singleSliceTempVcfPath(const std::string& outputPrefix) { return outputPrefix + ".resume_part.vcf"; }

bool resumeCheckpointExists(const ResumeCheckpointPaths& paths)
{
    return fs::exists(paths.json) && fs::exists(paths.vcf) && fs::exists(paths.processedLoci);
}

void removeCheckpointFiles(const ResumeCheckpointPaths& paths)
{
    for (const std::string* path : { &paths.json, &paths.vcf, &paths.processedLoci })
    {
        if (std::remove(path->c_str()) != 0 && errno != ENOENT)
        {
            spdlog::warn("Could not remove checkpoint file {}", *path);
        }
    }
}

void removeResumeTempFiles(const std::string& outputPrefix)
{
    std::vector<std::string> paths{ singleSliceTempJsonPath(outputPrefix), singleSliceTempVcfPath(outputPrefix) };

    // Per-contig temps are named <prefix>.contig<N>.{json,vcf}; which contigs exist depends on the run that
    // wrote them, so scan the output directory rather than assuming this run's contig set.
    const fs::path prefixPath(outputPrefix);
    const fs::path directory = prefixPath.parent_path().empty() ? fs::path(".") : prefixPath.parent_path();
    const std::string stem = prefixPath.filename().string() + ".contig";
    boost::system::error_code directoryError;
    for (fs::directory_iterator entry(directory, directoryError), end; entry != end; ++entry)
    {
        const std::string name = entry->path().filename().string();
        if (name.compare(0, stem.size(), stem) != 0)
        {
            continue;
        }
        const std::string remainder = name.substr(stem.size());
        const std::size_t dotPos = remainder.find('.');
        if (dotPos == 0 || dotPos == std::string::npos)
        {
            continue;
        }
        const std::string extension = remainder.substr(dotPos);
        if ((extension != ".json" && extension != ".vcf")
            || remainder.find_first_not_of("0123456789") != dotPos)
        {
            continue;
        }
        paths.push_back(entry->path().string());
    }

    for (const std::string& path : paths)
    {
        if (std::remove(path.c_str()) != 0 && errno != ENOENT)
        {
            spdlog::warn("Could not remove temporary file {}", path);
        }
    }
}

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

    // Checkpointed VCF lines are matched back to their locus through the catalog's variant ids (a catalog
    // may set VariantId explicitly), so those have to be unique across loci too: a collision would file one
    // locus's lines under another, and both the verification counts and the rebuilt temp files would follow
    // the wrong locus. Such a catalog already emits two VCF records carrying the same VARID without
    // --resume.
    std::unordered_set<std::string> seenVariantIds;
    for (const LocusDescription& locusDescription : catalog)
    {
        for (const std::string& variantId : locusDescription.variantIds())
        {
            if (!seenVariantIds.insert(variantId).second)
            {
                throw std::runtime_error(
                    "--resume cannot be used with this catalog: the VariantId '" + variantId
                    + "' appears on more than one locus, and resume matches VCF records back to their locus "
                    "by VariantId. Give each variant a unique VariantId, or run without --resume.");
            }
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
    signature["CompressOutputFiles"] = params.compressOutputFiles();
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
    signature["GenotypeQualityModelVersion"]
        = params.genotypeQualityModel() ? params.genotypeQualityModel()->version : std::string();
    return signature.dump(2);
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
        return "the checkpoint's ResumeInfo record could not be parsed";
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

std::string readCheckpointRunSignature(const ResumeCheckpointPaths& paths)
{
    return readRunSignatureFromFile(paths.json, checkpointIsCompressed(paths.json));
}

ResumeLoadResult loadResumeCheckpoint(
    const ResumeCheckpointPaths& paths, const SampleParameters& sampleParams,
    const LocusDescriptionCatalog& catalog, const std::string& outputPrefix, bool demultiplexPerContig)
{
    const std::string runSignature = readRunSignatureFromFile(paths.json, checkpointIsCompressed(paths.json));

    const bool compressedJson = checkpointIsCompressed(paths.json);
    const bool compressedVcf = checkpointIsCompressed(paths.vcf);

    ResumeLoadResult result;
    std::vector<ProcessedLocus> processedLoci = readProcessedLoci(paths.processedLoci);

    std::unordered_map<std::string, std::int32_t> contigOfLocus;
    std::unordered_map<std::string, std::string> locusOfVariant;
    contigOfLocus.reserve(catalog.size());
    for (const LocusDescription& locusDescription : catalog)
    {
        contigOfLocus[locusDescription.locusId()] = locusDescription.locusContigIndex();
        for (const std::string& variantId : locusDescription.variantIds())
        {
            locusOfVariant[variantId] = locusDescription.locusId();
        }
    }

    // Every locus the checkpoint claims is done has to be in the catalog: without its contig the record
    // cannot be routed to a temp file, and it would silently vanish from the final output.
    for (const ProcessedLocus& processedLocus : processedLoci)
    {
        if (contigOfLocus.find(processedLocus.locusId) == contigOfLocus.end())
        {
            throw std::runtime_error(
                "checkpoint locus " + processedLocus.locusId + " is not in the catalog, so the checkpoint "
                "does not belong to this run");
        }
    }

    // Verification pass: which loci actually have records on disk. The checkpoint writer flushes both
    // record files before it appends to processed_loci, so every listed locus should be present in them.
    // Checking anyway catches a checkpoint that lost buffered record writes some other way (a machine
    // crash rather than a killed process), which would otherwise drop those loci from the final output
    // without genotyping them again.
    std::unordered_set<std::string> lociWithJsonRecord;
    const bool jsonHeaderFound = scanJsonRecords(
        paths.json, compressedJson,
        [&lociWithJsonRecord](const std::string& locusId, const std::string&) {
            lociWithJsonRecord.insert(locusId);
        });
    if (!jsonHeaderFound)
    {
        throw ResumeCheckpointCorruptError("no LocusResults record in " + paths.json);
    }

    // Counted per locus rather than merely present/absent: a multi-variant locus writes one VCF line per
    // variant, so losing one of several lines has to be caught as well.
    std::unordered_map<std::string, std::size_t> vcfLineCountOfLocus;
    scanLines(paths.vcf, compressedVcf, [&](const std::string& line) {
        if (line.empty() || line[0] == '#')
        {
            return;
        }
        auto locusIt = locusOfVariant.find(variantIdOfVcfLine(line));
        if (locusIt != locusOfVariant.end())
        {
            ++vcfLineCountOfLocus[locusIt->second];
        }
    });

    // A lost write can only take a suffix of each record file, so stop at the first locus whose records
    // are not all there rather than trying to keep the ones after it.
    const std::size_t verifiedCount = static_cast<std::size_t>(std::distance(
        processedLoci.begin(),
        std::find_if(processedLoci.begin(), processedLoci.end(), [&](const ProcessedLocus& processedLocus) {
            const bool jsonRecordMissing = processedLocus.jsonRecordCount > 0
                && lociWithJsonRecord.find(processedLocus.locusId) == lociWithJsonRecord.end();
            auto vcfCountIt = vcfLineCountOfLocus.find(processedLocus.locusId);
            const std::size_t vcfLinesFound = vcfCountIt == vcfLineCountOfLocus.end() ? 0 : vcfCountIt->second;
            return jsonRecordMissing || vcfLinesFound < processedLocus.vcfLineCount;
        })));
    if (verifiedCount < processedLoci.size())
    {
        spdlog::warn(
            "Resume: the checkpoint lists {} finished loci but its record files only hold the first {}; the "
            "rest will be genotyped again",
            add_commas_at_thousands(processedLoci.size()), add_commas_at_thousands(verifiedCount));
        processedLoci.resize(verifiedCount);
    }

    // Check that this run's slices can actually be appended to (see firstUnusableProcessedLocus) before
    // rewriting anything. This only bites when a --threads > 1 checkpoint is resumed at --threads 1, whose
    // single coordinate sweep cannot append to contigs finished out of order.
    //
    // Refusing rather than trimming is the point: the rewrite below would cut the checkpoint files down to
    // the reusable part, destroying the rest of an interrupted run's work for good. Telling the operator to
    // re-run with --threads > 1 keeps every finished locus, and is useless advice if the records are already
    // gone by the time they read it.
    const std::size_t usableCount
        = firstUnusableProcessedLocus(processedLoci, catalog, contigOfLocus, demultiplexPerContig);
    if (usableCount < processedLoci.size())
    {
        throw std::runtime_error(
            "--resume: this checkpoint holds " + add_commas_at_thousands(processedLoci.size())
            + " finished loci from a run that genotyped contigs in parallel, and only "
            + add_commas_at_thousands(usableCount)
            + " of them can be appended to in the single coordinate sweep --threads 1 uses. Resume with "
            "--threads > 1 to keep all of them, or delete the .unfinished files to start over.");
    }

    for (const ProcessedLocus& processedLocus : processedLoci)
    {
        result.doneLocusIds.insert(processedLocus.locusId);
    }

    // Rewrite the processed-loci list so all three checkpoint files agree exactly on the finished set: a
    // later resume must not see a locus listed here whose records this one just dropped.
    const std::string trimmedProcessedLociPath = paths.processedLoci + ".trimmed";
    {
        std::ofstream trimmed(trimmedProcessedLociPath, std::ios::out | std::ios::trunc);
        if (!trimmed)
        {
            throw std::runtime_error("Failed to open file for writing: " + trimmedProcessedLociPath);
        }
        for (const ProcessedLocus& processedLocus : processedLoci)
        {
            trimmed << formatProcessedLocusLine(
                processedLocus.locusId, processedLocus.jsonRecordCount, processedLocus.vcfLineCount);
        }
        trimmed.flush();
        if (!trimmed)
        {
            throw std::runtime_error(
                "Failed to write " + trimmedProcessedLociPath + " (" + std::strerror(errno) + ")");
        }
    }
    renameOver(trimmedProcessedLociPath, paths.processedLoci);

    SliceRebuilder rebuilder(outputPrefix, sampleParams, demultiplexPerContig);

    // Pass 1: the JSON records. Everything past the last processed locus, or written for a locus whose
    // processed_loci entry never made it to disk, is dropped and will be genotyped again.
    const std::string trimmedJsonPath = paths.json + ".trimmed";
    {
        OutputFile trimmed(trimmedJsonPath, compressedJson);
        trimmed.write(checkpointJsonHeader(sampleParams, runSignature));

        bool firstRecord = true;
        scanJsonRecords(
            paths.json, compressedJson,
            [&](const std::string& locusId, const std::string& recordText) {
                if (result.doneLocusIds.find(locusId) == result.doneLocusIds.end())
                {
                    return;
                }
                if (!firstRecord)
                {
                    trimmed.write(", ");
                }
                trimmed.write(recordText);
                firstRecord = false;
                rebuilder.addJsonRecord(contigOfLocus.at(locusId), recordText);
            });
        trimmed.close();
        result.hasExistingJsonRecords = !firstRecord;
    }
    renameOver(trimmedJsonPath, paths.json);

    // Pass 2: the VCF lines, keyed to loci through the catalog's variant ids.
    const std::string trimmedVcfPath = paths.vcf + ".trimmed";
    {
        OutputFile trimmed(trimmedVcfPath, compressedVcf);
        trimmed.write(vcfDocumentHeader(sampleParams.id()));

        scanLines(paths.vcf, compressedVcf, [&](const std::string& line) {
            if (line.empty() || line[0] == '#')
            {
                return;
            }
            const std::string variantId = variantIdOfVcfLine(line);
            auto locusIt = locusOfVariant.find(variantId);
            if (locusIt == locusOfVariant.end())
            {
                // A line whose variant is not in the catalog, or one truncated before its INFO column.
                // Either way it cannot be placed, and its locus will be genotyped again.
                return;
            }
            const std::string& locusId = locusIt->second;
            if (result.doneLocusIds.find(locusId) == result.doneLocusIds.end())
            {
                // The locus's processed_loci entry never reached disk, so its JSON record was dropped
                // above; drop its VCF lines too and let it be genotyped again.
                return;
            }
            trimmed.write(line);
            rebuilder.addVcfLine(contigOfLocus.at(locusId), line);
        });
        trimmed.close();
    }
    renameOver(trimmedVcfPath, paths.vcf);

    result.rebuiltSlices = rebuilder.finish();
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
    ResumeCheckpointPaths paths, const SampleParameters& sampleParams, const std::string& runSignature, bool append,
    bool hasExistingJsonRecords, std::size_t abortAfterLoci)
    : paths_(std::move(paths))
    , compressJson_(checkpointIsCompressed(paths_.json))
    , compressVcf_(checkpointIsCompressed(paths_.vcf))
    , firstJsonRecord_(!append || !hasExistingJsonRecords)
    , abortAfterLoci_(abortAfterLoci)
{
    if (!append)
    {
        // Start the checkpoint from scratch, discarding anything left by an earlier run that this one is
        // not resuming.
        {
            OutputFile json(paths_.json, compressJson_);
            json.write(checkpointJsonHeader(sampleParams, runSignature));
        }
        {
            OutputFile vcf(paths_.vcf, compressVcf_);
            vcf.write(vcfDocumentHeader(sampleParams.id()));
        }
        std::ofstream processedLoci(paths_.processedLoci, std::ios::out | std::ios::trunc);
        if (!processedLoci)
        {
            throw std::runtime_error("Failed to open checkpoint file for writing: " + paths_.processedLoci);
        }
    }

    thread_ = std::thread([this]() { run(); });
}

ResumeCheckpointWriter::~ResumeCheckpointWriter()
{
    try
    {
        close();
    }
    catch (...)
    {
        // Never throw from a destructor. A checkpoint write failure that reaches here has already been
        // reported by close() on the normal path.
    }
}

void ResumeCheckpointWriter::recordLocus(std::string locusId, std::string jsonText, std::string vcfText)
{
    const std::size_t kMaxQueuedLoci = 4096;

    std::unique_lock<std::mutex> lock(mutex_);
    queueNotFull_.wait(lock, [this]() { return queue_.size() < kMaxQueuedLoci || stopRequested_ || writeError_; });
    if (writeError_)
    {
        std::exception_ptr error = writeError_;
        lock.unlock();
        std::rethrow_exception(error);
    }
    queue_.push_back(QueuedLocus{ std::move(locusId), std::move(jsonText), std::move(vcfText) });
    lock.unlock();
    queueNotEmpty_.notify_one();
}

void ResumeCheckpointWriter::close()
{
    if (closed_)
    {
        return;
    }
    closed_ = true;

    {
        std::lock_guard<std::mutex> lock(mutex_);
        stopRequested_ = true;
    }
    queueNotEmpty_.notify_all();
    queueNotFull_.notify_all();
    if (thread_.joinable())
    {
        thread_.join();
    }

    if (writeError_)
    {
        std::rethrow_exception(writeError_);
    }
}

void ResumeCheckpointWriter::run()
{
    while (true)
    {
        std::vector<QueuedLocus> batch;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            queueNotEmpty_.wait(lock, [this]() { return !queue_.empty() || stopRequested_; });
            if (queue_.empty())
            {
                return;  // stop requested and nothing left to write
            }
            batch.assign(
                std::make_move_iterator(queue_.begin()), std::make_move_iterator(queue_.end()));
            queue_.clear();
        }
        queueNotFull_.notify_all();

        try
        {
            writeBatch(batch);
        }
        catch (...)
        {
            std::lock_guard<std::mutex> lock(mutex_);
            writeError_ = std::current_exception();
            queueNotFull_.notify_all();
            return;
        }
    }
}

void ResumeCheckpointWriter::writeBatch(const std::vector<QueuedLocus>& batch)
{
    std::string jsonText;
    std::string vcfText;
    std::string locusIdText;
    for (const QueuedLocus& locus : batch)
    {
        if (!locus.jsonText.empty())
        {
            if (!firstJsonRecord_)
            {
                jsonText += ", ";
            }
            jsonText += locus.jsonText;
            firstJsonRecord_ = false;
        }
        vcfText += locus.vcfText;
        locusIdText
            += formatProcessedLocusLine(locus.locusId, locus.jsonText.empty() ? 0 : 1, countLines(locus.vcfText));
    }

    // The record files are flushed before any locus id is appended, so a locus listed in processed_loci is
    // always already present in both of them. That ordering is what lets a resume trust the list.
    appendToFile(paths_.json, jsonText, compressJson_);
    appendToFile(paths_.vcf, vcfText, compressVcf_);
    appendToFile(paths_.processedLoci, locusIdText, false);

    checkpointedLocusCount_ += batch.size();
    if (abortAfterLoci_ != 0 && checkpointedLocusCount_ >= abortAfterLoci_)
    {
        // Test hook (--internal-abort-after-loci): leave the process exactly as an external kill would,
        // with the checkpoint on disk and nothing unwound.
        spdlog::warn("Aborting after {} checkpointed loci (--internal-abort-after-loci)", checkpointedLocusCount_);
        std::_Exit(137);
    }
}

}
