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
 *  \brief Test suite for the Solveur_ws class.
 */

#include "gtest/gtest.h"
#include "solver_ws.h"

//==============================================================================

/*
TEST(SolverWS, SolverWS)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, calcWS)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, calcWSWithFields)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, getError)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, updateCurrentDef)
{
}
*/

//==============================================================================

TEST(SolverWS, init)
{
    auto dataTree = DataTree::getDefault();

    dataTree.set("constraints/q10t", 10.0);
    dataTree.set("constraints/q20t", 800.0);
    dataTree.set("constraints/q30t", 4000.0);

    SolverWS solver(dataTree);

    solver.init();

#ifdef WS_IGNORE_Q10_CONSTRAINT
    ASSERT_EQ(solver.dim, 2);
#else
    ASSERT_EQ(solver.dim, 3);
#endif
}


//==============================================================================

/*
TEST(SolverWS, nextIter)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, finalize)
{
}
*/

//==============================================================================

/*
TEST(SolverWS, calc)
{
}
*/
