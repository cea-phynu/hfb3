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

#ifndef SOLVER_HFB_GRADIENT_H
#define SOLVER_HFB_GRADIENT_H

/** \file
 *  \brief Headers for the SolverHFBGradient class.
 */

#include "discrete.h"
#include "interaction.h"
#include "global.h"
#include "multi.h"
#include "solver.h"
#include "multipole_operators.h"
#include "geometrical_operators.h"
#include "fragments.h"

class DataTree;

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description, optional default value and type for each input or output.
#define SOLVER_HFB_GRADIENT_VALID_KEYS \
{ "solver/gradient/cvgTarget"      , "Convergence target value"                                       , "1e-1", "D" }, \
{ "solver/gradient/maxIter"        , "Maximum number of iterations"                                   , "50"  , "I" }, \
{ "solver/gradient/cvgTargetLambda", "Convergence target value for lambda-iterations"                 , "1e-5", "D" }, \
{ "solver/gradient/maxIterLambda"  , "Maximum number of iterations for lambda-iterations"             , "10"  , "I" }, \
{ "solver/gradient/randomSeed"     , "Seed used for the generation of initial random U and V matrices", "1337", "I" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief A Gradient solver for the HFB equations.
 *
 *  This class implements a solver for the HFB equations using the Gradient method.
 */

class SolverHFBGradient : public Solver
{
public:
  explicit SolverHFBGradient(const std::string &filename);             // #TEST#
  explicit SolverHFBGradient(const DataTree &);                        // #TEST#
  explicit SolverHFBGradient(const DataTree &, State);                 // #TEST#

  //============================================================================

  void init(void);                                                     // #TEST#
  bool nextIter(void);                                                 // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// A Interaction object.
  Interaction interaction;

  /// A Discrete object.
  Discrete discrete;

  /// A MultipoleOperators object.
  MultipoleOperators multipoleOperators;

  /// A Fragments object.
  Fragments fragments;

  /// A GeometricalOperators object.
  GeometricalOperators geometricalOperators;

  //============================================================================

  /// Total binding energy.
  double energy = 1e99;

  /// Best energy found.
  double bestEnergy = 1e99;

  /// Convergence value.
  double convergence = 1e99;

  /// Convergence target value.
  double cvgTarget = 1e99;

  /// Target value for Lambda variation.
  double cvgTargetLambda = 1e-09;

  /// Maximum lambda-iterations.
  INT maxIterLambda = 20;

  /// Activate to enable the calculation of fragment properties at each HFB iteration.
  bool fragInLoop = false;

  //============================================================================
  //============================================================================
  //============================================================================

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

private:

  /// Local number of iterations.
  int localNbIter = 0;

  bool gradIter(Multi<arma::mat> oG, arma::vec mu, arma::vec nu);
  bool calcUV(Multi<arma::mat> &_U, Multi<arma::mat> &_V, Multi<arma::mat> &Uj, Multi<arma::mat> &Vj);
  Multi<arma::mat> getGradient();
  bool getCholeskyDecomposition(Multi<arma::mat> &_G, Multi<arma::mat> &L);

  void finalize(Multi<arma::mat> &G, double coeff = 1.);

  void calcRhoKappa(Multi<arma::mat> &U, Multi<arma::mat> &V, Multi<arma::mat> &r, Multi<arma::mat> &k);
  static void convFromSPtoQP(arma::mat &sp11, arma::mat &sp20, arma::mat &qp11, arma::mat &qp20, arma::mat &U, arma::mat &V);
  static void convFromQPtoSP(arma::mat &sp11, arma::mat &sp20, arma::mat &qp11, arma::mat &qp20, arma::mat &U, arma::mat &V);
  void calcConstraintsQP(Multi<arma::mat> &U, Multi<arma::mat> &V, Multi<arma::mat> &, arma::vec &);
  arma::mat calcConstraints(UINT iso, UVEC &idxHO, arma::mat &M);

  void initRandomWF(void);

  /// Seed for the generation of random initial U and V matrices.
  INT randomSeed = 1337;

  /// U matrix
  Multi<arma::mat> U;

  /// V matrix
  Multi<arma::mat> V;

  /// rho matrix
  Multi<arma::mat> rho;

  /// kappa matrix
  Multi<arma::mat> kappa;

  /// G matrix
  Multi<arma::mat> G;

  /// Lambda iteration number.
  UINT lambdaIter = 0;

  /// The list of disabled interaction
  IVEC disabledInteraction;

  /// The length of one HFB iteration.
  double iterLength = 0.0;

  /// Parameters for heavy ball
  arma::vec mu, nu;

  //============================================================================
  //============================================================================
  //============================================================================
};

#endif     // SOLVER_HFB_GRADIENT_H
