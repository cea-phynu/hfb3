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
 *  \brief Headers for the SolverBasis class.
 */

#include "solver_basis.h"
#include "gradientwalk.h"
#include "solver_hfb_gradient.h"
#include "solver_hfb_broyden.h"
#include "solver_alternator.h"
#include "basis.h"
#include "plot.h"
#include "tools.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > SolverBasis::validKeys =
  {
    { "solver/basis/sphericalBasis", "Force the use of a spherical basis (b_r = b_z)"       , "false", "B" },
    { "solver/basis/cvgTarget"     , "Convergence target value"                             , "1e-4" , "D" },
    { "solver/basis/maxIter"       , "Maximum number of iterations"                         , "40"   , "I" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

SolverBasis::SolverBasis(const std::string &filename) : SolverBasis(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree object.
 */

SolverBasis::SolverBasis(const DataTree &_dataTree) : SolverBasis(_dataTree, State(_dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree and a State.
 */

SolverBasis::SolverBasis(const DataTree &_dataTree, const State &_state) :
  Solver(_dataTree, _state)
{
  DBG_ENTER;

  dataTree.get(cvgTarget, "solver/basis/cvgTarget", true);
  dataTree.get(maxIter,   "solver/basis/maxIter",   true);

  // number of basis parameters to optimize ? 1ct: [b_r, b_z], 2ct: [d_0, b_r, b_z]
  dim = state.basis.dMax == 1 ? 2 : 3;

  // br = bz for a spherical basis
  if (sphericalBasis) dim = 1;

  ASSERT( !sphericalBasis || (state.basis.dMax == 1), "Using a spherical 2ct basis !?");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the basis optimization loop.
 */

void SolverBasis::init()
{
  DBG_ENTER;

  startTime = Tools::clock();

  arma::vec previousParams;

  if (dim == 1) previousParams = {state.basis.b_r};

  if (dim == 2) previousParams = {state.basis.b_r, state.basis.b_z};

  if (dim == 3) previousParams = {state.basis.b_r, state.basis.b_z, state.basis.d_0};

  // TODO: remove hardcoded values

  vmin  = arma::ones(dim) * 1.00;

  vmax  = arma::ones(dim) * 3.50;

  vstep = arma::ones(dim) * 0.05;

  if (dim == 3)
  {
    vmin(2)  =  5.0;
    vmax(2)  = 15.0;
    vstep(2) =  0.5;
  }

  minimizer.xMin           = vmin;
  minimizer.xMax           = vmax;
  minimizer.xStep          = vstep;
  minimizer.initialCoords  = previousParams;

  minimizer.init();

  Tools::mesg("SolBas", info(false));

  DBG_LEAVE;

}

//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Calculate the next iteration.
 */

bool SolverBasis::nextIter()
{
  DBG_ENTER;

  status = Solver::ITERATING;

  arma::vec next;
  next = minimizer.getEval();

  if (next.empty())
  {
    Tools::mesg("SolBas", PF_GREEN("no points left to eval => minimum found"));

    status = Solver::CONVERGED;

    finalize();

    DBG_RETURN(false);
  }

  bool valid;
  valid = calcHFB(next);
  minimizer.addEval(next, -value, valid);
  Tools::mesg("SolBas", getHistTable());

  status = Solver::ITERATING;

  nbIter++;

  if (nbIter >= maxIter)
  {
    Tools::mesg("SolBas", PF_RED("maximum number of iterations reached (%d)", nbIter));

    status = Solver::MAXITER;

    finalize();
    DBG_RETURN(false);
  }

  double convergence = minimizer.getConvergence();

  if (convergence < cvgTarget)
  {
    Tools::mesg("SolBas", PF_GREEN("target value reached %e < %e", convergence, cvgTarget));

    status = Solver::CONVERGED;

    finalize();
    DBG_RETURN(false);
  }

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

void SolverBasis::finalize(void)
{
  DBG_ENTER;

  state = bestState;

  Tools::mesg("SolBas", "End of optimization");
  // Tools::mesg("SolBas", getHistTable());

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the next iteration.
 */

bool SolverBasis::calcHFB(const arma::vec &basisParams, const std::string &label)
{
  DBG_ENTER;

  //=======================
  //== basis preparation ==
  //=======================

  double d_0           = state.basis.d_0;
  double b_r           = state.basis.b_r;
  double b_z           = state.basis.b_z;
  INT    nOscil        = state.basis.nOscil;
  INT    n_zMaxImposed = state.basis.n_zMaxImposed;
  double g_q           = state.basis.g_q;

  if (dim == 1)
  {
    b_r = basisParams(0);
    b_z = b_r;

  }
  else if (dim == 2)
  {
    b_r = basisParams(0);
    b_z = basisParams(1);
  }
  else if (dim == 3)
  {
    b_r = basisParams(0);
    b_z = basisParams(1);
    d_0 = basisParams(2);
  }
  else
  {
    ERROR("wrong number of basis parameters");
  }

  Basis newBasis = Basis(d_0, b_r, b_z, nOscil, n_zMaxImposed, g_q);

  Tools::mesg("SolBas", PF_GREEN("HFB calculation with basis parameters " + newBasis.getNiceInfo()));

  state = bestState.convertTo(newBasis);

  // state.getDataTree().save(PF("run%04d.hfb3", iter));

  //======================
  //== HFB calculation ==
  //======================

  DataTree customDataTree = dataTree;
  double cvgTarget = 1e-6;

  dataTree.get(cvgTarget, "solver/broyden/cvgTarget");
  customDataTree.set("solver/broyden/cvgTarget", cvgTarget * 100.0);

  SolverAlternator solverAlternator(customDataTree, state);
  solverAlternator.init();
  while (solverAlternator.nextIter());
  state = solverAlternator.state;

  value = state.converged ? state.totalEnergy: 1e99;

  histCoords = arma::resize(histCoords, dim, nbIter + 1);
  histValues = arma::resize(histValues, nbIter + 1,   1);
  histIters  = arma::resize(histIters, nbIter + 1,   1);

  histCoords.col(nbIter) = basisParams;
  histValues(nbIter) = value;
  histIters(nbIter) = solverAlternator.nbIter;

  newBest = false;


  if ((value < bestValue) && state.converged)
  {
    bestValue = value;
    newBest = true;
    bestState = state;
  }

  // Tools::mesg("SolBas", PF_MAGENTA("end of HFB: (basis: " + Tools::vecToStr(basisParams) + PF(", Ehfb: %.3e)", value)));

  if (useBokeh)
  {
    // Plot::save(PF("plot%04d.html", iter));

    Plot::clear("Ene HFB [MeV]");
    Plot::clear("Convergence");

    Plot::slot(4);
    Plot::curve("#iter basis", "Best Ene HFB [MeV]", nbIter, value);
  }

  if (solverAlternator.status == Solver::CONVERGED)
  {
    Tools::mesg("SolBas",
                PF("#it: %03d ", nbIter) +
                Tools::vecToStr(basisParams) + " " +
                PF("ene: ") + PF_GREEN("%9.3f", value) + " " +
                label + (newBest ? PF_GREEN(" *") : ""));
  }
  else
  {
    Tools::mesg("SolBas",
                PF("#it: %03d ", nbIter) +
                Tools::vecToStr(basisParams) + " " +
                PF_RED("NOT CONVERGED ") + " " +
                label);
  }

  bestState.calculationLength = Tools::clock() - startTime;

#ifdef DEBUG_META
  DBG_RETURN(!std::isnan(value));
#else
  DBG_RETURN(state.converged);
#endif

}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string SolverBasis::getHistTable(void) const
{
  DBG_ENTER;

  std::string result;
  result += TABLE_CENTER + TABLE_BLUE + "Iter" + TABLE_TD + TABLE_YELLOW
            + "Basis params. " + ((dim == 3) ? "[br, bz, d0]" : "[br, bz]")+ TABLE_TD
            + "Ehfb"             + TABLE_TD
            + "#iter"            + TABLE_TD
            + TABLE_TR;

  UINT idMin = histValues.index_min();

  for (UINT i = 0; i < histValues.n_rows; i++)
  {
    std::string eneStr = PF("%.3f", histValues(i));
    std::string iterStr = PF("%d", histIters(i));

    if (histValues(i) > 1e98) eneStr = TABLE_RED + "none";

    if (i == idMin)
      result += TABLE_BLUE + TABLE_RIGHT + PF("* %d", i) + TABLE_TD;
    else
      result += TABLE_GREEN + TABLE_RIGHT + PF("%d", i) + TABLE_TD;

    result += Tools::vecToStr(histCoords.col(i)) + TABLE_TD;
    result += TABLE_RIGHT + eneStr + TABLE_TD;
    result += TABLE_RIGHT + iterStr + TABLE_TR;
  }

  result += TABLE_TD + TABLE_TD + TABLE_TD + TABLE_RED + PF("tot: %d", arma::accu(histIters)) + TABLE_TR;

  DBG_RETURN(Tools::printTable(result));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string SolverBasis::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"state.", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"status", Solver::statusStr[status]},
    }, true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"SolverBasis", ""},
      {"state.", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"dim.  ", Tools::infoStr(dim)},
      {"status", Solver::statusStr[status]},
      {"vmin  ", Tools::vecToStr(vmin)},
      {"vmax  ", Tools::vecToStr(vmax)},
      {"vstep ", Tools::vecToStr(vstep)},
    }, false);
  }

  DBG_RETURN(result);
}
