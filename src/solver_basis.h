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

#ifndef SOLVER_BASIS_H
#define SOLVER_BASIS_H

/** \file
 *  \brief Headers for the SolverBasis class.
 */

#include "global.h"
#include "gradientwalk.h"
#include "solver.h"
#include "interaction.h"
#include "discrete.h"
#include "mixing.h"
#include "multi.h"
#include "multipole_operators.h"
#include "solver_hfb_gradient.h"
#include "solver_hfb_broyden.h"

class DataTree;

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description, optional default value and type for each input or output.
#define SOLVER_BASIS_VALID_KEYS \
{ "solver/basis/sphericalBasis", "Force the use of a spherical basis (b_r = b_z)"       , "false", "B" }, \
{ "solver/basis/cvgTarget"     , "Convergence target value"                             , "1e-4" , "D" }, \
{ "solver/basis/maxIter"       , "Maximum number of iterations"                         , "40"   , "I" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief The basis parameters optimizer.
 *
 *  This class represents the optimizer of the basis parameters.
 */

class SolverBasis : public Solver
{
public:

  explicit SolverBasis(const std::string &filename);                   // #TEST#
  explicit SolverBasis(const DataTree &, const State &);               // #TEST#
  explicit SolverBasis(const DataTree &);                              // #TEST#

  //============================================================================

  void init(void);                                                     // #TEST#
  bool nextIter(void);                                                 // #TEST#
  void finalize(void);                                                 // #TEST#
  bool calcHFB(const arma::vec &, const std::string &label = "");      // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string getHistTable(void) const;                          // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// Number of parameters to consider.
  INT dim = 0;

  /// Lowest value reached so far.
  double bestValue = 1e10;

  /// Determinant of the covariant matrix.
  double determinant = 1e99;

  /// Density of the covariant matrix.
  double density = 1e99;

  /// Use a spherical basis (b_r = b_z) ?
  bool sphericalBasis = false;

  /// Total binding energy.
  double energy = 1e99;

  /// Best total binding energy.
  double bestEnergy = 1e99;

  /// Convergence value.
  double convergence = 1e99;

  /// Convergence target value.
  double cvgTarget = 1e99;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  /// An instance of the minimizer
  GradientWalk minimizer;

  /// Did the last calculation improve the best State ?
  bool newBest = false;

  /// History of the coordinates tried.
  arma::mat previousCoords;

  /// History of the values obtained.
  arma::vec previousEnergy;

  /// History of the values obtained.
  IVEC previousIter;

  /// Minimal values of the basis parameters.
  arma::vec vmin;

  /// Maximal values of the basis parameters.
  arma::vec vmax;

  /// Mesh discretization steps.
  arma::vec vstep;
};

#endif // SOLVER_BASIS_H
