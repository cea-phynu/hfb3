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

#ifndef SOLVER_WS_H
#define SOLVER_WS_H

/** \file
 *  \brief Headers for the SolverWS class.
 */

#include "global.h"
#include "solver.h"
#include "gradientwalk.h"
#include "interaction.h"
#include "state.h"
#include "discrete.h"
#include "multi.h"
#include "multipole_operators.h"
#include "state.h"

class DataTree;

/** \brief The solver of the WS equation.
 *
 *  This class represents the core solver of the WS equation.
 */

class SolverWS : public Solver
{
public:

  explicit SolverWS(const std::string filename);                       // #TEST#
  explicit SolverWS(const DataTree &);                                 // #TEST#
  explicit SolverWS(const DataTree &_dataTree, State _state);          // #TEST#

  //============================================================================

  void calcWS();                                                       // #TEST#
  void calcWSWithInteraction(Interaction &interactionWS);

  double getError(double scalingFactor = 5.0);                         // #TEST#
  std::string niceStr(Multi<double> &def);
  void updateCurrentDef(const arma::vec &v);                           // #TEST#

  void init(void);                                                     // #TEST#
  bool nextIter(void);                                                  // #TEST#
  void finalize(void);                                                 // #TEST#
  State calc();                                                        // #TEST#
  void bokehPlot(void);

  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// The interaction.
  Interaction wsInteraction;

  /// A MultipoleOperators object.
  MultipoleOperators multipoleOperators;

  /// Woods-Saxon individual energies.
  arma::mat indivEnerWS;

  /// Woods-Saxon individual states.
  arma::cube indivStateWS;

  /// Woods-Saxon individual occupations.
  arma::mat indivOccupWS;

  /// Persistent Discrete instance.
  Discrete discrete;

  /// The best deformation parameters found by a Minimization instance.
  Multi<double> bestDef;

  /// The best error found by a Minimization instance until now.
  double bestError = 1e99;

  /// Best State found by a Minimizer instance.
  State bestState;

  /// Number of parameters to consider.
  UINT dim = 0;

  /// Initial constraint for $\langle \beta_{20}\rangle$.
  double beta20tInit = 0.1;

  /// Initial constraint for $\langle Q_{30}\rangle$.
  double q30tInit = 1000.0;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  /// Original constraints.
  std::map<std::string, Constraint> originalConstraints;

  /// An instance of the current deformation for calcWS.
  Multi<double> currentDef;

  /// An instance of GradientWalk.
  GradientWalk gradientWalk;

  /// Minimal values of the basis parameters.
  arma::vec vmin;

  /// Maximal values of the basis parameters.
  arma::vec vmax;

  /// Mesh discretization steps.
  arma::vec vstep;

  /// Using a temporary constraint ?
  bool temporaryConstraint = false;
};

#endif // SOLVER_WS_H
