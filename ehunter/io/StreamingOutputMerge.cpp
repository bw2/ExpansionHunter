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

std::string readWholeFile(const std::string& path)
{
    std::ifstream inFile(path, std::ios::in | std::ios::binary);
    if (!inFile)
    {
        throw std::runtime_error("Failed to open file: " + path);
    }
    std::ostringstream buffer;
    buffer << inFile.rdbuf();
    return buffer.str();
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
}

void mergeRegionJsonFiles(const std::string& finalPath, const std::vector<std::string>& regionTempPaths)
{
    // This merge depends on IterativeJsonWriter's exact wrapper format (single source of truth:
    // io/IterativeJsonWriter.cpp). Each region file is:
    //     "{\n  \"SampleParameters\": <json>,\n  \"LocusResults\": {" + <body> + "\n  }\n}\n"
    // where <body> is the raw record text (records separated by ", "; empty for a region with no
    // records). We keep region 0's prefix (SampleParameters is identical across regions), join the
    // non-empty bodies with ", ", and append the wrapper suffix.
    static const std::string kBodyStartMarker = "\"LocusResults\": {";
    static const std::string kBodyEndMarker = "\n  }\n}";

    std::ofstream outFile;
    boost::iostreams::filtering_ostream outStream;
    openMergedOutput(finalPath, outFile, outStream);

    bool wroteAnyBody = false;
    for (size_t fileIndex = 0; fileIndex != regionTempPaths.size(); ++fileIndex)
    {
        const std::string content = readWholeFile(regionTempPaths[fileIndex]);

        const size_t markerPos = content.find(kBodyStartMarker);
        if (markerPos == std::string::npos)
        {
            throw std::runtime_error("Missing \"LocusResults\" marker in file: " + regionTempPaths[fileIndex]);
        }
        const size_t bodyStart = markerPos + kBodyStartMarker.size();

        const size_t bodyEnd = content.rfind(kBodyEndMarker);
        if (bodyEnd == std::string::npos || bodyEnd < bodyStart)
        {
            throw std::runtime_error("Missing closing markers in file: " + regionTempPaths[fileIndex]);
        }

        // Region 0's prefix (everything up to and including the body-start marker) is canonical.
        if (fileIndex == 0)
        {
            outStream << content.substr(0, bodyStart);
        }

        const std::string body = content.substr(bodyStart, bodyEnd - bodyStart);
        if (!body.empty())
        {
            if (wroteAnyBody)
            {
                outStream << ", ";
            }
            outStream << body;
            wroteAnyBody = true;
        }
    }

    outStream << "\n  }\n}\n";
}

}
