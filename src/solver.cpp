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

#include "solver.h"

/** \file
 *  \brief Methods of the Solver class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > Solver::validKeys =
  {
    { "solver/forceNonEmptyKappa", "Force non-empty kappa matrices at the start of HFB iterations", "True", "B" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

std::vector<std::string> Solver::statusStr =
  {
    "Iterating", "Starting", "Converged", "Diverged", "Maxiter", "End of scheme"
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

Solver::Solver(const std::string &filename) : Solver(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

Solver::Solver(const DataTree &_dataTree) : Solver(_dataTree, State(_dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

Solver::Solver(const DataTree &_dataTree, State _state) : state(_state), dataTree(_dataTree)
{
  DBG_ENTER;

  bestState = state;
  dataTree.get(forceNonEmptyKappa, "solver/forceNonEmptyKappa", true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the optimization loop.
 */

void Solver::init(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the next iteration.
 */

bool Solver::nextIter(void)
{
  DBG_ENTER;

  DBG_RETURN(true);
}

