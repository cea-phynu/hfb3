//==============================================================================
// HFB3
// Copyright CEA, DAM F-91297 Arpajon, France
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//==============================================================================

/** \file
 *  \brief Test suite for the IOhfb3 class.
 */

#include "gtest/gtest.h"
#include "io_hfb3.h"
#include "datatree.h"

//==============================================================================

/*
TEST(IOhfb3, IOhfb3)
{
}
*/

//==============================================================================

TEST(IOhfb3, fromContent)
{
  std::string content =
R"toto(#!HFB3!#

[action]
basisOptimization F
saveResultFiles   True
)toto";

  DataTree dt = IOhfb3().fromContent(content);

  bool value = true;
  dt.get(value, "action/basisOptimization", true);

  ASSERT_EQ(value, false);

  value = false;
  dt.get(value, "action/saveResultFiles", true);

  ASSERT_EQ(value, true);
}


