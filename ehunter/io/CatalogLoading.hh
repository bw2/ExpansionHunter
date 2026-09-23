//
// Expansion Hunter
// Copyright 2016-2019 Illumina, Inc.
// All rights reserved.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
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

#pragma once

#include <string>
#include <vector>

#include "thirdparty/json/json.hpp"

#include "core/Common.hh"
#include "core/Parameters.hh"
#include "core/Reference.hh"
#include "locus/LocusSpecification.hh"

namespace ehunter
{

LocusDescriptionCatalog loadLocusDescriptions(
	const ProgramParameters& params, const Reference& reference);
RegionCatalog convertLocusDescriptionsToLocusSpecs(
    LocusDescriptionCatalog& locusDescriptionCatalog, const HeuristicParameters& heuristicParams, Reference& reference);

// The optional "KnownMotifs" field of a catalog locus record: motifs to treat as known at this locus in
// --output-motif-composition, as written in the catalog (empty if the field is absent). Entries the motif composition
// calculation cannot use are skipped there, with a warning (see selectCatalogKnownMotifs). Throws if the field is
// not an array of strings.
std::vector<std::string> decodeKnownMotifs(const nlohmann::json& locusJson, const std::string& locusId);

}
