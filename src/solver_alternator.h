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

#ifndef SOLVER_ALTERNATOR_H
#define SOLVER_ALTERNATOR_H

/** \file
 *  \brief Headers for the SolverAlternator class.
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
#include "solver_ws.h"

class DataTree;

/** \brief The basis parameters optimizer.
 *
 *  This class represents the optimizer of the basis parameters.
 */

class SolverAlternator : public Solver
{
public:

  /// Scheme item types.
  enum { WOODSSAXON, BROYDEN, GRADIENT, DO, ENDDO };

  explicit SolverAlternator(const std::string &filename);              // #TEST#
  explicit SolverAlternator(const DataTree &);                         // #TEST#
  explicit SolverAlternator(const DataTree &, State);                  // #TEST#

  //============================================================================

  void init(void);                                                     // #TEST#
  bool nextIter(void);                                                 // #TEST#
  void finalize(void);                                                 // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// Scheme for the alternating solvers.
  std::string scheme = "";

  /// Scheme for the alternating solvers (vector<int> version).
  std::vector<INT> schemeList = {};

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// Current scheme index.
  INT schemeId = 0;

  /// Optional scheme loop index.
  INT schemeLoopIndex = -1;

  /// Flag for solver initialization.
  bool solverInit = true;

  /// Id of the final solver.
  INT finalId = -1;

  //============================================================================
  //============================================================================
  //============================================================================

  /// The Woods-Saxon solver.
  SolverWS solverWS;

  /// The Gradient HFB solver.
  SolverHFBGradient solverHFBGradient;

  // The Broyden HFB solver.
  SolverHFBBroyden solverHFBBroyden;

  //============================================================================
  //============================================================================
  //============================================================================

private:

};

#endif // SOLVER_ALTERNATOR_H
