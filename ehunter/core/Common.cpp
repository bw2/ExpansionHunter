//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>,
//         Mitch Bekritsky <mbekritsky@illumina.com>, Richard Shaw
// Concept: Michael Eberle <meberle@illumina.com>
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

#include "core/Common.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <regex>

#include <sys/resource.h>

using std::string;

namespace ehunter
{

Sex decodeSampleSex(const std::string& encoding)
{
    if (encoding == "male")
    {
        return Sex::kMale;
    }
    else if (encoding == "female")
    {
        return Sex::kFemale;
    }
    else
    {
        throw std::invalid_argument(encoding + " is invalid sex; must be either male or female");
    }
}

std::ostream& operator<<(std::ostream& out, Sex sex)
{
    switch (sex)
    {
    case Sex::kFemale:
        out << "Female";
        break;
    case Sex::kMale:
        out << "Male";
        break;
    }

    return out;
}

std::ostream& operator<<(std::ostream& out, ReadType readType)
{
    switch (readType)
    {
    case ReadType::kFlanking:
        out << "FLANKING";
        break;
    case ReadType::kRepeat:
        out << "INREPEAT";
        break;
    case ReadType::kSpanning:
        out << "SPANNING";
        break;
    case ReadType::kOther:
        out << "OTHER";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, AlleleCount alleleCount)
{
    switch (alleleCount)
    {
    case AlleleCount::kOne:
        out << "One";
        break;
    case AlleleCount::kTwo:
        out << "Two";
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, NumericInterval numericInterval)
{
    out << numericInterval.start() << "-" << numericInterval.end();
    return out;
}

bool isURL(const std::string& path)
{
    static const std::regex url_regex(".*?://.*");
    return std::regex_match(path, url_regex);
}

std::time_t currentEpochSeconds()
{
    return std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
}

std::string formatLocalTimestamp(std::time_t epochSeconds)
{
    std::tm localTime{};
    localtime_r(&epochSeconds, &localTime);
    char buffer[32];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%S", &localTime);
    return std::string(buffer);
}

std::string formatRuntime(std::time_t durationSeconds)
{
    // A backward wall-clock step (NTP/manual adjustment) during a long run can make this negative;
    // clamp to 0 rather than print a misleading truncated value (e.g. "-1s" for a ~-1h duration).
    const long long total = std::max<long long>(0, static_cast<long long>(durationSeconds));
    const long long hours = total / 3600;
    const long long minutes = (total % 3600) / 60;
    const long long seconds = total % 60;

    std::ostringstream out;
    if (hours > 0)
    {
        out << hours << "h ";
    }
    if (hours > 0 || minutes > 0)
    {
        out << minutes << "m ";
    }
    out << seconds << "s";
    return out.str();
}

double peakRssMemoryMB()
{
    struct rusage usage
    {
    };
    getrusage(RUSAGE_SELF, &usage);
    // ru_maxrss is in KB on Linux but bytes on macOS/BSD.
#ifdef __APPLE__
    const double bytes = static_cast<double>(usage.ru_maxrss);
#else
    const double bytes = static_cast<double>(usage.ru_maxrss) * 1024.0;
#endif
    return std::round(bytes / (1024.0 * 1024.0) * 10.0) / 10.0;
}

}
