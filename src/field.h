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

#ifndef FIELD_H
#define FIELD_H

/** \file
 *  \brief Headers for the Field class.
 */

#include "global.h"
#include "generic.h"
#include "basis.h"
#include "multi.h"
#include "state.h"
#include "wspot.h"

class State;

/** \brief Parent class for the derived Field* classes.
 *
 * This class must be derived from.
 */

class Field : public Generic
{
public :

  /// The parameters for each field instance
  typedef std::map<std::string, double> Parameters;

  enum
  {
    DIRECT,
    EXCHANGE,
    PAIRING
  };

  explicit Field(Parameters _parameters, State *);                     // #TEST#
  virtual void calcField(void);                                        // #TEST#
  virtual void calcEnergy(void);                                       // #TEST#
  virtual void setDef(const Multi<double> &_def);                      // #TEST#
  virtual WSPot getWSPot(void);                                        // #TEST#
  virtual void clear(void);                                            // #TEST#
  double getEnergy(INT iso, INT type);                                 // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string ramInfo() const;
  const std::string timeInfo() const;

  //============================================================================

  /// A pointer to a State instance.
  State *state_ptr;

  /// A reference to a State instance.
  State &state;

  /// The field parameters.
  Parameters parameters;

  /// A reference to a Basis instance.
  Basis &basis;

  /// A reference to the rho matrices.
  Multi<arma::mat> &rho;

  /// A reference to the kappa matrices.
  Multi<arma::mat> &kappa;

  //============================================================================

  /// Number of Gauss-Laguerre nodes.
  INT nGLA = 0;

  /// Number of Gauss-Hermite nodes.
  INT nGHE = 0;

  /// Number of Gauss-Legendre nodes.
  INT nGLE = 0;

  /// The calculating length.
  double calculatingLength = -2.0;

  /// Size of objects
  std::map<std::string, UINT> ramTable;

  /// Time taken
  std::map<std::string, double> timeTable;

  /// Should this field contribute to the total hamiltonian ?
  bool contributeToEnergy = true;

  /// Contributions to the total hamiltonian and/or delta.
  Multi<arma::mat> field;

  /// Contributions to the total energy.
  Multi<double> energy;

  /// The name of the field.
  std::string name = "";

  /// The short name of the field.
  std::string shortName = "";

  /// Field is enabled ?
  IVEC enabled = {1, 1, 1};

  /// Arbitrary factor for the energy.
  double energyFactor = 1.0;

  /// Field must be cleared by Interaction (at every HFB iteration).
  bool mustBeCleared = true;

  /// Warning message.
  std::string warningStr;

  /// The type names.
  static const std::vector<std::string> typeStr;
};

#endif // FIELD_H
