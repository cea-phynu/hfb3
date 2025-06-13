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

#ifndef GRADIENTWALK_H
#define GRADIENTWALK_H

#include "global.h"
#include "tools.h"

/** \file
 *  \brief Headers for the GradientWalk class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Simple implementation of a gradient walk minimization method.
 */

class GradientWalk
{
public:

  const std::string info(bool isShort = USE_SHORT_INFO) const;

  void init(void);
  void addEval(const arma::vec &_xEval, const double _yEval, bool valid = true);
  arma::vec getEval(void);
  double getConvergence(void);
  const arma::vec getOptimum(void);

  /// space
  arma::vec xMin;
  arma::vec xMax;
  arma::vec xStep;
  arma::vec initialCoords;

private:
  arma::mat xEval;
  arma::vec yEval;

  // direction
  IVEC candidateDirection;
  IVEC direction;

  // start
  bool start;
  bool restart;

  UINT dim;
  UINT bestEval;
  UINT nEval;

  double convergence;

  arma::vec startEval();
};

#endif
