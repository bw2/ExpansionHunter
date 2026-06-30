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

#pragma once

#include <string>

namespace ehunter
{

// Result of comparing a sequence against a perfect tiling of a motif.
struct MotifPurity
{
    long matchedBases = 0; // bases of `sequence` equal to the motif base at the same (mod-length) phase
    long totalBases = 0;   // length of `sequence`
};

// Count how many bases of `sequence` match a perfect tiling of `motif` anchored at sequence[0]
// (i.e. sequence[i] is compared to motif[i % motif.length()]). Both sides are upper-cased, so
// soft-masked (lowercase) reference still matches; an 'N' counts as a mismatch. Returns {0, 0}
// when `motif` is empty. The caller's purity is matchedBases / totalBases.
//
// This is a flat, gap-free comparison: a single indel relative to the motif shifts the phase for
// all downstream bases, so it is a purity heuristic, not a true alignment. Used for the reference
// repeat region and for fast-path spanning reads; the full genotyper instead derives matched vs
// total bases from the graph alignment's match/mismatch operations.
MotifPurity motifTilingPurity(const std::string& sequence, const std::string& motif);

}
