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
 *  \brief Test suite for the Action class.
 */

#include "gtest/gtest.h"
#include "action.h"

//==============================================================================

TEST(Action, Action)
{
  DataTree dataTree("misc/data/test_2ct_256Fm.hfb3");
  Action action(dataTree);
  ASSERT_TRUE(dataTree == action.dataTree);
}

//==============================================================================

TEST(Action, Action_)
{
  std::string filename = "misc/data/test_2ct_256Fm.hfb3";
  Action action(filename);
  DataTree dataTree(filename);
  ASSERT_TRUE(dataTree == action.dataTree);
}
