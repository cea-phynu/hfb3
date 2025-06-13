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

#include "mixing.h"
#include "tools.h"

/** \file
 *  \brief Methods of the Mixing class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > Mixing::validKeys =
  {
    { "solver/broyden/mixing/size"        , "Number of configurations", "8"  , "I" },
    { "solver/broyden/mixing/linearFactor", "Linear mixing factor"    , "0.1", "D" },
    { "solver/broyden/mixing/mixingfactor", "General mixing factor"   , "0.6", "D" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

Mixing::Mixing(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================


/** The constructor from a DataTree.
 */

Mixing::Mixing(const DataTree &dataTree)
{
  DBG_ENTER;

  INT _mixingSize = 8;
  double _linearFactor = 0.1;
  double _mixingFactor = 0.6;

  dataTree.get(_mixingSize  , "solver/broyden/mixing/size"        , true);
  dataTree.get(_linearFactor, "solver/broyden/mixing/linearFactor", true);
  dataTree.get(_mixingFactor, "solver/broyden/mixing/mixingfactor", true);

  (*this) = Mixing(_mixingSize, _linearFactor, _mixingFactor);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================


/** The constructor from the parameter values.
 */

Mixing::Mixing(INT _mixingSize, double _linearFactor, double _mixingFactor) :
  mixingSize(_mixingSize),
  linearFactor(_linearFactor),
  mixingFactor(_mixingFactor)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Store a new vector and do the mixing.
 */

bool Mixing::newVec(const arma::vec &vec)
{
  DBG_ENTER;

  if (F.empty())
  {
    UINT size = vec.n_rows;
    F  = arma::zeros(size);
    dF = arma::zeros(size, mixingSize);
    dx = arma::zeros(size, mixingSize);
    x  = vec;
  }

  arma::vec y = vec;

  // We don't permute the columns like in the reference: instead, we overwrite
  // them in an (anti) circular way
  INT l = (mixingSize - 1 - m) % mixingSize;

  if (l < 0) l += mixingSize;

  dF.col(l) = -F;
  F         = y - x;
  dF.col(l) += F;

  if ((m < mixingSize) || (mixingSize == 1))
  {
    if (firstLinear)
    {
      // Tools::mesg("Mixing", "using linear method (factor: " + PF("%5.3f", linearFactor) + ")");
      firstLinear = false;
    }

    // Tools::info("y", y);
    // Tools::info("x", x);

    y = (1 - linearFactor) * x + linearFactor * y;
  }
  else
  {
    if (firstBroyden)
    {
      // Tools::mesg("Mixing", "using Broyden method (factor: " + PF("%5.3f, mixingSize: %1d", mixingFactor, mixingSize) + ")");
      firstBroyden = false;
    }

    dFdF = dF.t() * dF;
    dFF = dF.t() * F;
    dFdF += w0 * w0 * arma::diagmat(dFdF);

    bool retry = false;

    try
    {
      gamma = arma::solve(dFdF, dFF, arma::solve_opts::no_approx);
    }
    catch (const std::runtime_error &e)
    {
      retry = true;
    }

    if (retry)
    {
      try
      {
        gamma = arma::solve(dFdF, dFF, arma::solve_opts::fast);
      }
      catch (const std::runtime_error &e)
      {
        Tools::warning("no solution found in Mixing::newSol()");
        DBG_RETURN(false);
      }
    }

    //    gamma.print("poids");
    y = x + mixingFactor * F - dx * arma::shift(gamma, 1) - mixingFactor * dF * gamma;
  }

  dx.col(l)     = y - x;
  x             = y;
  m += 1;

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the mixed vector.
 */

const arma::vec & Mixing::getVec()
{
  DBG_ENTER;

  DBG_RETURN(x);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Mixing::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"m.size", Tools::infoStr(mixingSize)},
      {"lfact.", Tools::infoStr(linearFactor)},
      {"mfact.", Tools::infoStr(mixingFactor)},
    }, true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"Mixing", ""},
      {"m.size", Tools::infoStr(mixingSize)},
      {"lfact.", Tools::infoStr(linearFactor)},
      {"mfact.", Tools::infoStr(mixingFactor)},
    }, false);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check if linear mode is still active.
 */

bool Mixing::isLinearMode(void)
{
  DBG_ENTER;

  DBG_RETURN(m < mixingSize);
}
