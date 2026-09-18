//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#include "io/StreamingOutputMerge.hh"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>

namespace ehunter
{

namespace
{

// Push a gzip compressor onto the stream iff the path ends in "gz", matching the mechanism used by
// IterativeJsonWriter/IterativeVcfWriter.
void openMergedOutput(const std::string& finalPath, std::ofstream& outFile, boost::iostreams::filtering_ostream& outStream)
{
    outFile.open(finalPath, std::ios::out | std::ios::binary);
    if (!outFile)
    {
        throw std::runtime_error("Failed to open file: " + finalPath);
    }

    if (finalPath.size() > 2 && finalPath.substr(finalPath.size() - 2) == "gz")
    {
        outStream.push(boost::iostreams::gzip_compressor());
    }
    outStream.push(outFile);
}

// Flush the merged output and fail loudly if any of the writing went wrong. Without this a full disk or a
// write error would leave a truncated final output that the caller reports as a successful run -- and, under
// --resume, would take the checkpoint files needed to recover with it.
void closeMergedOutput(
    const std::string& finalPath, std::ofstream& outFile, boost::iostreams::filtering_ostream& outStream)
{
    outStream.flush();
    outStream.reset();
    outFile.flush();
    const bool failed = !outFile;
    outFile.close();
    if (failed || !outFile)
    {
        throw std::runtime_error("Failed to write " + finalPath + " (" + std::strerror(errno) + ")");
    }
}

// Read the first `count` bytes of an already-open file, from its current position.
std::string readBytes(std::ifstream& inFile, std::size_t count)
{
    std::string buffer(count, '\0');
    inFile.read(&buffer[0], static_cast<std::streamsize>(count));
    buffer.resize(static_cast<std::size_t>(inFile.gcount()));
    inFile.clear();
    return buffer;
}

// Copy bytes [begin, end) of `inFile` to `outStream` in bounded chunks, so a multi-gigabyte region file
// never has to be held in memory at once.
void copyRange(std::ifstream& inFile, std::streamoff begin, std::streamoff end, std::ostream& outStream)
{
    const std::size_t kChunkSize = 1 << 20;
    std::vector<char> chunk(kChunkSize);
    inFile.clear();
    inFile.seekg(begin);
    std::streamoff remaining = end - begin;
    while (remaining > 0)
    {
        const std::streamsize wanted
            = static_cast<std::streamsize>(std::min<std::streamoff>(remaining, static_cast<std::streamoff>(kChunkSize)));
        inFile.read(chunk.data(), wanted);
        const std::streamsize got = inFile.gcount();
        if (got <= 0)
        {
            break;
        }
        outStream.write(chunk.data(), got);
        remaining -= got;
    }
}

}  // namespace

void mergeRegionVcfFiles(const std::string& finalPath, const std::vector<std::string>& regionTempPaths)
{
    std::ofstream outFile;
    boost::iostreams::filtering_ostream outStream;
    openMergedOutput(finalPath, outFile, outStream);

    for (size_t fileIndex = 0; fileIndex != regionTempPaths.size(); ++fileIndex)
    {
        std::ifstream inFile(regionTempPaths[fileIndex], std::ios::in | std::ios::binary);
        if (!inFile)
        {
            throw std::runtime_error("Failed to open file: " + regionTempPaths[fileIndex]);
        }

        // The first region file contributes its header and data lines; subsequent files contribute
        // only their data lines (lines not starting with '#'), dropping the duplicate header.
        const bool isFirstFile = (fileIndex == 0);
        std::string line;
        while (std::getline(inFile, line))
        {
            if (isFirstFile || line.empty() || line.front() != '#')
            {
                outStream << line << "\n";
            }
        }
    }

    closeMergedOutput(finalPath, outFile, outStream);
}

void mergeRegionJsonFiles(
    const std::string& finalPath, const std::vector<std::string>& regionTempPaths, const std::string& runInfoJson)
{
    // This merge depends on IterativeJsonWriter's exact wrapper format (single source of truth:
    // io/IterativeJsonWriter.cpp). Each region file is:
    //     "{\n  \"SampleParameters\": <json>,\n  \"LocusResults\": {" + <body> + "\n  },\n  \"RunInfo\": <json>\n}\n"
    // where <body> is the raw record text (records separated by ", "; empty for a region with no
    // records). We keep region 0's prefix (SampleParameters is identical across regions), join the
    // non-empty bodies with ", ", discard every region's own RunInfo (each reflects only that region
    // worker's local start/finish time, not the true run-wide completion), and append the caller-supplied
    // runInfoJson instead.
    static const std::string kBodyStartMarker = "\"LocusResults\": {";
    static const std::string kBodyEndMarker = "\n  },\n  \"RunInfo\"";

    std::ofstream outFile;
    boost::iostreams::filtering_ostream outStream;
    openMergedOutput(finalPath, outFile, outStream);

    // The header and the footer sit at known ends of the file, so each region file is located by reading
    // only those two ends and then streamed through in chunks. A region file holds a whole contig's records
    // (every locus, under --resume at --threads 1), which can run to gigabytes, so it must never be slurped
    // into memory.
    //
    // The windows start small and grow rather than being fixed, because the footer carries the run's whole
    // command line inside its RunInfo record: with a long one (a --locus list naming thousands of loci, say)
    // a fixed window would miss the marker and fail the run at its very last step, after all the genotyping.
    // Growth stops at kMaxMarkerWindow, well past any command line the OS will accept, so a genuinely
    // markerless file reports that instead of being read into memory whole.
    const std::size_t kInitialMarkerWindow = 64 * 1024;
    const std::size_t kMaxMarkerWindow = 64u * 1024 * 1024;

    bool wroteAnyBody = false;
    for (size_t fileIndex = 0; fileIndex != regionTempPaths.size(); ++fileIndex)
    {
        const std::string& path = regionTempPaths[fileIndex];
        std::ifstream inFile(path, std::ios::in | std::ios::binary);
        if (!inFile)
        {
            throw std::runtime_error("Failed to open file: " + path);
        }
        inFile.seekg(0, std::ios::end);
        const std::streamoff fileSize = inFile.tellg();
        inFile.seekg(0);

        std::string header;
        size_t markerPos = std::string::npos;
        for (std::size_t window = kInitialMarkerWindow; markerPos == std::string::npos; window *= 2)
        {
            inFile.clear();
            inFile.seekg(0);
            header = readBytes(inFile, static_cast<std::size_t>(std::min<std::streamoff>(fileSize, window)));
            markerPos = header.find(kBodyStartMarker);
            if (markerPos == std::string::npos
                && (static_cast<std::streamoff>(header.size()) >= fileSize || window >= kMaxMarkerWindow))
            {
                throw std::runtime_error("Missing \"LocusResults\" marker in file: " + path);
            }
        }
        const std::streamoff bodyStart = static_cast<std::streamoff>(markerPos + kBodyStartMarker.size());

        std::streamoff bodyEnd = -1;
        for (std::size_t window = kInitialMarkerWindow; bodyEnd < 0; window *= 2)
        {
            const std::streamoff footerStart
                = std::max<std::streamoff>(bodyStart, fileSize - static_cast<std::streamoff>(window));
            inFile.clear();
            inFile.seekg(footerStart);
            const std::string footer = readBytes(inFile, static_cast<std::size_t>(fileSize - footerStart));
            const size_t footerMarkerPos = footer.rfind(kBodyEndMarker);
            if (footerMarkerPos != std::string::npos)
            {
                bodyEnd = footerStart + static_cast<std::streamoff>(footerMarkerPos);
                break;
            }
            if (footerStart <= bodyStart || window >= kMaxMarkerWindow)
            {
                throw std::runtime_error("Missing closing markers in file: " + path);
            }
        }
        if (bodyEnd < bodyStart)
        {
            throw std::runtime_error("Missing closing markers in file: " + path);
        }

        // Region 0's prefix (everything up to and including the body-start marker) is canonical.
        if (fileIndex == 0)
        {
            outStream << header.substr(0, static_cast<std::size_t>(bodyStart));
        }

        if (bodyEnd > bodyStart)
        {
            if (wroteAnyBody)
            {
                outStream << ", ";
            }
            copyRange(inFile, bodyStart, bodyEnd, outStream);
            wroteAnyBody = true;
        }
    }

    outStream << "\n  },\n  \"RunInfo\": " << runInfoJson << "\n}\n";

    closeMergedOutput(finalPath, outFile, outStream);
}

}
