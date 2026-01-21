//
// ExpansionHunter
// Copyright 2016-2021 Illumina, Inc.
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
// Adapted from REViewer (Copyright 2020 Illumina, Inc., GPL-3.0 license)
//

#pragma once

#include <list>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>

#include "graphcore/Path.hh"

#include "reviewer/Aligns.hh"
#include "reviewer/Origin.hh"

namespace ehunter
{
namespace reviewer
{

/// Types of visual features in lane plots
enum class FeatureType
{
    kRect,
    kRectWithLeftBreak,
    kRectWithRightBreak,
    kLine,
    kArrows,
    kVerticalLine
};

/// A visual feature in a lane plot segment
struct Feature
{
    Feature(FeatureType type, int length, std::string fill, std::string stroke)
        : type(type)
        , length(length)
        , fill(std::move(fill))
        , stroke(std::move(stroke))
    {
    }
    FeatureType type;
    int length;
    boost::optional<std::string> label;
    std::string fill;
    std::string stroke;
};

/// A segment of aligned features
struct Segment
{
    Segment(int start, std::vector<Feature> features, double opacity)
        : start(start)
        , features(std::move(features))
        , opacity(opacity)
    {
        end = start;
        for (const auto& feature : this->features)
        {
            end += feature.length;
        }
    }

    int start;
    int end;
    std::vector<Feature> features;
    double opacity;
};

/// A lane containing segments
struct Lane
{
    Lane(int height, std::vector<Segment> segments)
        : height(height)
        , segments(std::move(segments))
    {
    }

    int height;
    std::vector<Segment> segments;
};

/// A complete lane plot (one per haplotype)
using LanePlot = std::vector<Lane>;

/// Generate visual blueprint for SVG rendering
/// @param paths Haplotype paths
/// @param fragById Map of fragments by ID
/// @param fragAssignment Fragment origin assignment
/// @param fragPathAlignsById Fragment path alignments
/// @return Vector of lane plots (one per haplotype)
std::vector<LanePlot> generateBlueprint(
    std::vector<graphtools::Path> paths,
    const FragById& fragById,
    const FragAssignment& fragAssignment,
    const FragPathAlignsById& fragPathAlignsById);

}  // namespace reviewer
}  // namespace ehunter
