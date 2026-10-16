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

#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "io/BamletWriter.hh"

#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "locus/LocusAnalyzer.hh"
#include "locus/LocusFindings.hh"
#include "locus/LocusSpecification.hh"
#include "sample/HtsFileSeeker.hh"
#include "sample/MateExtractor.hh"

namespace ehunter
{

// `contigInfo` (the reads file's contigs) and `typicalReadLength` are used only by --output-motif-composition: each
// thread opens its own reference reader with `contigInfo`, and `typicalReadLength` decides which motif sizes qualify.
SampleFindings htsSeekingSampleAnalysis(
    const ProgramParameters& programParams, const HeuristicParameters& heuristicParams,
    const RegionCatalog& regionCatalog, BamletWriterPtr alignmentWriter, const ReferenceContigInfo& contigInfo,
    int typicalReadLength);

// Genotypes one locus from the reads seeking mode collects for it: every primary paired read in the locus's target
// and off-target read-extraction regions, then the mates of far-apart pairs. Seeking mode genotypes every locus with
// it, and optimized-streaming mode calls it for loci that use OfftargetRegions, so their records match seeking mode's. `locusAnalyzers` holds exactly the
// locus's analyzer (AnalyzerFinder routes reads by index into this vector). Under --output-motif-composition the
// returned findings also carry the locus's MotifComposition, computed from the reads as they were before graph
// alignment (which reverse-complements reads in place), with `reference` and `typicalReadLength` as
// addMotifComposition uses them.
LocusFindings genotypeLocusWithSeekingModeReads(
    const ProgramParameters& programParams, std::vector<std::unique_ptr<locus::LocusAnalyzer>>& locusAnalyzers,
    htshelpers::HtsFileSeeker& htsFileSeeker, htshelpers::MateExtractor& mateExtractor,
    graphtools::AlignerSelector& alignerSelector, Reference& reference, int typicalReadLength);

}
