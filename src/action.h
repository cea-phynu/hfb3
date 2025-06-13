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

#ifndef ACTION_H
#define ACTION_H

#include "global.h"
#include "generic.h"
#include "datatree.h"
#include "state.h"

/** \file
 *  \brief Headers for the Action class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief A class implementing various high-level use cases.
 *
 * This class implements some of the main high-level use-cases of the solver.
 * These use-cases can be used as entry points or possible workflow examples.
 */

class Action : public Generic
{
  public:

  //============================================================================
  //============================================================================
  //============================================================================

  /// A structure to store a blocking try.
  struct BlockingTry
  {
    INT id;           ///< Id of the try.
    std::string name; ///< Name of the try.
    std::string neut; ///< Blocked neutron state.
    std::string prot; ///< Blocked proton state.
    double ene;       ///< Energy obtained.
    INT nbIter;       ///< Number of iterations performed.
    bool converged;   ///< Did the solver converge ?
  };

  //============================================================================
  //============================================================================
  //============================================================================

  // Constructors
  explicit Action(const DataTree &);
  explicit Action(const std::string &);

  // Entry point
  void run(void);

  // Actions
  void calcHFBblocking(void);
  void calcHFBquiet(void);
  void calcHFB(void);
  void calcWS(void);
  void calcWSHFB(void);
  void calcEnergies(void);
  void calcObservables(void);
  void calcBasis(void);
  void calcMultipoleExpansion(void);

  // Misc
  const std::string getBlockingHistTable(void) const;
  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// The main DataTree instance.
  DataTree dataTree;

  /// The final State
  State state;

  /// The blocking tries.
  std::list<BlockingTry> blockingTrials;

  /// Optimize the basis parameters ?
  bool basisOptimization = false;

  /// Save result files ?
  bool saveResultFiles = true;

  /// Name of the job.
  std::string jobName = "unnamed";

  /// Action to be done.
  std::string action = "ws_hfb";

  /// Number of blocking tries for each isospin value.
  int nbBlockingTrials = 3;
};

#endif // ACTION_H
