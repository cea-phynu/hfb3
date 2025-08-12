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
 *  \brief Headers for the SolverAlternator class.
 */

#include "solver_alternator.h"
#include "solver_hfb_gradient.h"
#include "solver_hfb_broyden.h"
#include "tools.h"
#include <cmath>
#include <cstdio>
#include <string>

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct> SolverAlternator::validKeys =
  {
    { "solver/alternator/maxIter", "Maximum number of iterations when alternating solvers",  "200", "I" },
    { "solver/alternator/scheme" , "Scheme for the alternating solvers"                   , "[GB]", "S" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

SolverAlternator::SolverAlternator(const std::string &filename) : SolverAlternator(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree object.
 */

SolverAlternator::SolverAlternator(const DataTree &_dataTree) : SolverAlternator(_dataTree, State(_dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree and a State.
 */

SolverAlternator::SolverAlternator(const DataTree &_dataTree, State _state) :
  Solver(_dataTree, _state),
  solverWS(_dataTree, _state),
  solverHFBGradient(_dataTree, _state),
  solverHFBBroyden(_dataTree, _state)
{
  DBG_ENTER;

  dataTree.get(scheme               , "solver/alternator/scheme"               , true);
  dataTree.get(maxIter              , "solver/alternator/maxIter"              , true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the basis optimization loop.
 */

void SolverAlternator::init()
{
  DBG_ENTER;

  startTime = Tools::clock();

  solverHFBBroyden.discrete = Discrete(&state.basis, Mesh::regular(-10.0, 0.0, -25.0, 10.0, 0.0, 25.0, 101, 1, 251));
  solverHFBGradient.discrete = Discrete(&state.basis, Mesh::regular(-10.0, 0.0, -25.0, 10.0, 0.0, 25.0, 101, 1, 251));

  Tools::mesg("SolAlt", info(false));

  schemeList = {};
  finalId = -1;

  for (auto c: scheme)
  {
    if      (c == 'W') schemeList.push_back(WOODSSAXON);
    else if (c == 'B') { schemeList.push_back(BROYDEN); finalId = schemeList.size() - 1; }
    else if (c == 'G') { schemeList.push_back(GRADIENT); finalId = schemeList.size() - 1; }
    else if (c == '[') schemeList.push_back(DO);
    else if (c == ']') schemeList.push_back(ENDDO);
    else
    {
      ERROR("Unknown scheme letter: '%c'", c);
    }
  }

  schemeId = 0;
  solverInit = true;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Calculate the next iteration.
 */

bool SolverAlternator::nextIter()
{
  DBG_ENTER;

  // DEBUG("scheme: %d", schemeId);
  status = Solver::ITERATING;

  switch (schemeList.at(schemeId))
  {
    case WOODSSAXON:
    {
      if (solverInit)
      {
        if (!state.empty())
        {
          Tools::mesg("SolAlt", "skipping the WS solver, since an initial state is given");
          solverInit = true;
          schemeId++;
          break;
        }

        Tools::mesg("SolAlt", "Starting the WS solver");

        solverWS.init();
        solverInit = false;
      }
      bool wsIterating = solverWS.nextIter();

      if (!wsIterating)
      {
        state = solverWS.state;
        solverInit = true;
        schemeId++;
      }

      // do not include WS iterations in the iteration total
      // nbIter++;

      break;
    }
    case BROYDEN:
    {
      if (solverInit)
      {
        Tools::mesg("SolAlt", "Starting the Broyden HFB solver");

        solverHFBBroyden.state = state;
        solverHFBBroyden.init();

        solverInit = false;
      }
      bool broydenIterating = solverHFBBroyden.nextIter();

      if (!broydenIterating)
      {
        state = solverHFBBroyden.state;

        // only converge if final solver
        if ((schemeId == finalId) && (solverHFBBroyden.status == Solver::CONVERGED))
        {
          status = Solver::CONVERGED;
          break;
        }

        solverInit = true;
        schemeId++;
      }
      nbIter++;
      break;
    }
    case GRADIENT:
    {
      if (solverInit)
      {
        if (!state.empty())
        {
          Tools::mesg("SolAlt", "skipping the Gradient HFB solver, since an initial state is given");
          solverInit = true;
          schemeId++;
          break;
        }

        Tools::mesg("SolAlt", "Starting the Gradient HFB solver");

        solverHFBGradient.state = state;
        solverHFBGradient.init();
        solverInit = false;
      }
      bool gradientIterating = solverHFBGradient.nextIter();

      if (!gradientIterating)
      {
        state = solverHFBGradient.state;

        // only converge if final solver
        if ((schemeId == finalId) && (solverHFBGradient.status == Solver::CONVERGED))
        {
          status = Solver::CONVERGED;
          break;
        }

        solverInit = true;
        schemeId++;
      }
      nbIter++;
      break;
    }
    case DO:
    {
      schemeLoopIndex = schemeId;
      schemeId++;
      break;
    }
    case ENDDO:
    {
      ASSERT(schemeLoopIndex >= 0, "Missing loop start '[' in scheme.");
      schemeId = schemeLoopIndex;
      schemeId++;
      break;
    }
  } // end of switch

  //============================================================================

  if (status == Solver::CONVERGED)
  {
    state.converged = true;
    DBG_RETURN(false);
  }

  if (nbIter >= maxIter)
  {
    status = Solver::MAXITER;
    state.converged = false;
    DBG_RETURN(false);
  }

  if (schemeId >= scheme.size())
  {
    status = Solver::SCHEME_END;
    state.converged = false;
    DBG_RETURN(false);
  }

  // INFO(PF("%d/%d ", nbIter, maxIter) + Solver::statusStr[status]);

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

void SolverAlternator::finalize(void)
{
  DBG_ENTER;

  Tools::mesg("SolAlt", "End of optimization");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string SolverAlternator::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"state ", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"status", Solver::statusStr[status]},
    }, true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"SolverAlternator", ""},
      {"state ", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"inter.", solverHFBBroyden.interaction.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"target", Tools::infoStr(cvgTarget)},
      {"status", Solver::statusStr[status]},
      {"scheme", scheme},
      {"momen.", solverHFBBroyden.multipoleOperators.info(true)},
    }, false);
  }

  DBG_RETURN(result);
}
