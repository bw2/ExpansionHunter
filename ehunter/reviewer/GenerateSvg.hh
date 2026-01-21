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

#include <string>
#include <vector>

#include "reviewer/LanePlot.hh"

namespace ehunter
{
namespace reviewer
{

/// Generate SVG file from lane plots
/// @param lanePlots Lane plots to render
/// @param outputPath Output file path
void generateSvg(const std::vector<LanePlot>& lanePlots, const std::string& outputPath);

}  // namespace reviewer
}  // namespace ehunter
