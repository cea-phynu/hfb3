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

#include "gradientwalk.h"
#include "tools.h"

/** \file
 *  \brief Methods of the GradientWalk class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

void GradientWalk::init()
{
  DBG_ENTER;

  start = true;
  restart = true;

  dim  = xMin.n_rows;
  nEval = 0;
  convergence = 1.0;
  bestEval = -1;

  ASSERT((dim > 0), "unsupported dimension");

  // Initial candidate Direction (simple gradient walk method)
  candidateDirection = arma::zeros<IVEC >(dim);
  candidateDirection(0) = 1;

  // Initial Direction (simple gradient walk method)
  direction = arma::zeros<IVEC >(dim);
  direction(0) = 0;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

// Add a new evaluated point
void GradientWalk::addEval(const arma::vec &_xEval, const double _yEval, bool valid)
{
  DBG_ENTER;

  // INFO(Tools::vecToStr(_xEval) + " -> " + PF("%f", _yEval));
  UINT sizeOld = xEval.n_cols;

  // Resizing (x,y,index)Eval
  xEval.resize(dim, sizeOld + 1);
  //INFO("xEval: " + Tools::infoStr(xEval));

  yEval.resize(sizeOld + 1);
  //INFO("yEval: " + Tools::infoStr(yEval));

  // Copy data
  xEval.col(sizeOld) = _xEval;

  if (valid)
  {
    yEval(sizeOld)     = _yEval;
  }
  else
  {
    if (bestEval > -1) yEval(sizeOld)     = yEval(bestEval) - 10.0; // dummy energy
    else yEval(sizeOld) = 0.0;
  }

  restart = false;
  nEval += 1;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

// If starting coords does not converge, use random points
arma::vec GradientWalk::startEval()
{
  DBG_ENTER;

  if (start)
  {
    start = false;
    DBG_RETURN(initialCoords);
  }

  arma::vec span = xMax - xMin;
  arma::vec rng = arma::vec().randu(dim);
  arma::vec res = arma::round(rng % span / xStep) % xStep + xMin;
  DBG_RETURN(res);
}

//==============================================================================
//==============================================================================
//==============================================================================

// Get points to evaluate
arma::vec GradientWalk::getEval(void)
{
  DBG_ENTER;

  if (restart)
  {
    DBG_RETURN(startEval());
  }

  arma::vec nextPoint;

  bestEval = yEval.index_max();

  if (bestEval == nEval - 1)
  {
    // new bestEval
    if (direction != candidateDirection)
    {
      direction = candidateDirection;
      // Tools::mesg("GradWa", PF_RED("New direction: " + Tools::ivecToStr(direction)));
    }

    if (nEval > 1)
    {
      UINT previousBestEval = yEval.head(nEval - 1).index_max();
      convergence = yEval(bestEval) - yEval(previousBestEval);
      // INFO("convergence: %e", convergence);
    }
  }

  if (!candidateDirection.empty())
  {
    arma::vec candidateCoords = xEval.col(bestEval) + xStep % direction;

    if (Tools::findColInMat(candidateCoords, xEval) == -1)
    {
      nextPoint = candidateCoords;
      DBG_RETURN(nextPoint);
    }
  }

  for (INT i = 0; i < dim; i++)
  {
    IVEC simpleDirection = arma::zeros<IVEC >(dim);
    simpleDirection(i) = 1;

    arma::vec candidateCoords = xEval.col(bestEval) + xStep % simpleDirection;

    if (Tools::findColInMat(candidateCoords, xEval) == -1)
    {
      nextPoint = candidateCoords;
      candidateDirection = simpleDirection;
      break;
    }

    simpleDirection(i) = -1;

    candidateCoords = xEval.col(bestEval) + xStep % simpleDirection;

    if (Tools::findColInMat(candidateCoords, xEval) == -1)
    {
      nextPoint = candidateCoords;
      candidateDirection = simpleDirection;
      break;
    }
  }

  if (nextPoint.empty())
  {
    convergence = 0.0;
  }

  DBG_RETURN(nextPoint);
}

//==============================================================================
//==============================================================================
//==============================================================================

double GradientWalk::getConvergence()
{
  return convergence;
}

//==============================================================================
//==============================================================================
//==============================================================================

const arma::vec GradientWalk::getOptimum(void)
{
  DBG_ENTER;

  if (bestEval != -1)
  {
    DBG_RETURN(xEval.col(bestEval));
  }

  DBG_RETURN(arma::vec());
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string GradientWalk::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"dim.  ", Tools::infoStr(dim)},
      {"nEval ", Tools::infoStr(nEval)},
    }, true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"GaussianProcess", ""},
      {"dim.  ", Tools::infoStr(dim)},
      {"nEval ", Tools::infoStr(nEval)},
      {"xMin  ", Tools::vecToStr(xMin)},
      {"xMax  ", Tools::vecToStr(xMax)},
    }, false);
  }

  DBG_RETURN(result);
}

