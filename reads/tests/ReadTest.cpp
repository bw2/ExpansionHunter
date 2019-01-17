//
// Expansion Hunter
// Copyright (c) 2016 Illumina, Inc.
//
// Author: Egor Dolzhenko <edolzhenko@illumina.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "reads/Read.hh"

#include "gtest/gtest.h"

using namespace ehunter;

TEST(ReadInitialization, TypicalCoreInfo_CoreInfoAddedToRead)
{
    ReadId readId("frag1", MateNumber::kSecondMate);
    Read read(readId, "ATTC");
    EXPECT_EQ("frag1", read.fragmentId());
}
