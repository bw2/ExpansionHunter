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

#include "alignment/AlignmentFilters.hh"

#include <list>
#include <vector>

#include "graphalign/GraphAlignmentOperations.hh"
#include "graphalign/LinearAlignmentOperations.hh"
#include "graphcore/PathOperations.hh"

#include "alignment/OperationsOnAlignments.hh"

using graphtools::GraphAlignment;
using graphtools::NodeId;
using graphtools::Operation;
using graphtools::OperationType;
using graphtools::Path;
using std::list;
using std::string;
using std::vector;

namespace ehunter
{

bool checkIfLocallyPlacedReadPair(
    boost::optional<GraphAlignment> readAlignment, boost::optional<GraphAlignment> mateAlignment,
    int kMinNonRepeatAlignmentScore)
{
    int nonRepeatAlignmentScore = 0;

    if (readAlignment)
    {
        nonRepeatAlignmentScore += scoreAlignmentToNonloopNodes(*readAlignment);
    }

    if (mateAlignment)
    {
        nonRepeatAlignmentScore += scoreAlignmentToNonloopNodes(*mateAlignment);
    }

    return nonRepeatAlignmentScore >= kMinNonRepeatAlignmentScore;
}

}
