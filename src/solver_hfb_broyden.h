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

#ifndef SOLVER_HFB_BROYDEN_H
#define SOLVER_HFB_BROYDEN_H

/** \file
 *  \brief Headers for the SolverHFBBroyden class.
 */

#include "global.h"
#include "solver.h"
#include "interaction.h"
#include "discrete.h"
#include "mixing.h"
#include "multi.h"
#include "multipole_operators.h"
#include "geometrical_operators.h"

class DataTree;

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define SOLVER_HFB_BROYDEN_VALID_KEYS \
{ "solver/broyden/cvgTarget"           , "Convergence target value."                                         , "1e-6" , "D" }, \
{ "solver/broyden/maxIter"             , "Maximum number of iterations"                                      , "200"  , "I" }, \
{ "solver/broyden/emptyKappaProtection", "Force non-empty kappa matrices"                                    , "False", "B" }, \
{ "solver/broyden/cvgTargetLambda"     , "Convergence target value for lambda-iterations"                    , "1e-5" , "D" }, \
{ "solver/broyden/maxIterLambda"       , "Maximum number of iterations for lambda-iterations"                , "10"   , "I" }, \
{ "solver/broyden/earlyLambdaMixing"   , "Mix Lagrange multipliers in the linear phase of the Broyden mixing", "True" , "B" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief A Broyden-mixing solver for the HFB equations.
 *
 *  This class implements a solver for the HFB equations using the Broyden mixing method.
 */

class SolverHFBBroyden : public Solver
{
public:

  explicit SolverHFBBroyden(const std::string &filename);              // #TEST#
  explicit SolverHFBBroyden(const DataTree &);                         // #TEST#
  explicit SolverHFBBroyden(const DataTree &, State);                  // #TEST#

  //============================================================================

  void calc(void);                                                     // #TEST#
  void init(void);                                                     // #TEST#
  bool nextIter(void);                                                 // #TEST#

  bool singleHFBiter(void);                                            // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  bool mixStates(void);                                                // #TEST#
  void bokehPlot(void);                                                // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// A Mixing instance.
  Mixing mixing;

  /// A Interaction instance.
  Interaction interaction;

  /// A Discrete instance.
  Discrete discrete;

  /// A MultipoleOperators instance.
  MultipoleOperators multipoleOperators;

  /// A GeometricalOperators instance.
  GeometricalOperators geometricalOperators;

  /// The last State instance.
  State lastState;

  //============================================================================

  /// Total binding energy.
  double energy = 1e99;

  /// Best total binding energy.
  double bestEnergy = 1e99;

  /// Convergence value.
  double convergence = 1e99;

  /// Convergence target value.
  double cvgTarget = 1e99;

  /// Target value for Lambda variation.
  double cvgTargetLambda = 1e-09;

  /// Maximum lambda-iterations.
  INT maxIterLambda = 20;

  /// Protection against empty kappa matrices.
  bool emptyKappaProtection = true;

  /// Total U matrix.
  Multi<arma::mat> U;

  /// Total V matrix.
  Multi<arma::mat> V;

  /// Particle energies (iso).
  Multi<arma::vec> partEne;

  /// Quasi-particle energies (iso).
  Multi<arma::vec> qpEne;

  /// Activate to enable the calculation of fragment properties at each HFB iteration.
  bool fragInLoop = false;

  /// Activate to mix the constraint Lagrange multipliers during the mixing linear mode.
  bool earlyLambdaMixing = true;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  /// Next Rho matrices.
  Multi<arma::mat> rhoNext;

  /// Next Kappa matrices.
  Multi<arma::mat> kappaNext;

  /// Lambda iteration number.
  IVEC lambdaIter = arma::zeros<IVEC >(2);

  /// Step for lambda.
  arma::vec lambdaStep = arma::ones(2);

  /// The list of disabled interaction.
  IVEC disabledInteraction;

  ///  The HF states organized in (iso, omega) blocks.
  Multi<arma::mat> hfStatesOmega;

  //============================================================================
  //============================================================================
  //============================================================================

  bool adjustConstraints(void);
};

#endif // SOLVER_HFB_BROYDEN_H
