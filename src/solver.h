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

#ifndef SOLVER_H
#define SOLVER_H

/** \file
 *  \brief Headers for the Solver class.
 */

#include <list>
#include "global.h"
#include "generic.h"
#include "state.h"
#include "datatree.h"

/** \brief Virtual class for the solvers.
 *
 * This class is a virtual class for classes implementing a solver.
 */

class Solver : public Generic
{
public:

  /// Return codes.
  enum {ITERATING, STARTING, CONVERGED, DIVERGED, MAXITER, SCHEME_END};

  explicit Solver(const std::string &filename);                        // #TEST#
  explicit Solver(const DataTree &);                                   // #TEST#
  explicit Solver(const DataTree &, State);                            // #TEST#

  virtual void init(void);
  virtual bool nextIter(void);

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// List of status types.
  static std::vector<std::string> statusStr;

  /// The current State.
  State state;

  /// A copy of the dataTree.
  DataTree dataTree;

  /// The best State so far.
  State bestState;

  /// The maximum number of iterations.
  INT maxIter = 200;

  /// The current iteration index.
  INT nbIter = 0;

  /// Status of the solver.
  INT status = 0;

  /// Starting time of the calculation [s].
  double startTime = 0.0;

  /// Length of the calculation [s].
  double calculationLength = 0.0;

  /// The target convergence value.
  double cvgTarget = 1e-3;

  /// The current value to be minimized.
  double value = 999.0;

  /// Must plot local densities ?
  bool plotDensities = false;

  /// Protection against empty kappa matrices.
  bool forceNonEmptyKappa = true;
};

#endif // SOLVER_H
