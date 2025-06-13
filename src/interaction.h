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

#ifndef INTERACTION_H
#define INTERACTION_H

/** \file
 *  \brief Headers for the Interaction class.
 */

#include <memory>
#include "global.h"
#include "generic.h"
#include "field.h"

class State;

/** \brief Store the Field instances.
 *
 * This class stores the instances of the Field* classes.
 */

class Interaction : public Generic
{
public :

  Interaction(const std::string &_interactionName = "D1S", State *_state = NULL); // #TEST#
  Interaction(DataTree dataTree, State *_state = NULL);                // #TEST#
  void clear(void);                                                    // #TEST#
  void calcFields(bool ignoreMissingRhoKappa = false);                 // #TEST#
  void calcEnergies(void);                                             // #TEST#
  const std::map<std::string, double> getEnergyContributions(void) const; // #TEST#
  const arma::mat getHamiltonianContributions(INT iso, INT type);      // #TEST#
  static const DataTree getInteractionDataTree(const std::string& interactionName);

  std::shared_ptr<Field> operator()(std::string name, INT id = 0);

  void addField(const std::string &fieldName);                         // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string getNiceInfo(void) const;                           // #TEST#
  std::vector<Field::Parameters> getParametersFromDataTree(DataTree &dataTree, const std::string &fieldName);

  //============================================================================
  //============================================================================
  //============================================================================

  /// The interaction name.
  std::string interactionName = "noname";

  /// A DataTree instance
  DataTree dataTree;

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// A pointer to a State instance.
  State *state;

  /// Total energies.
  arma::vec totalEnergy;

  /// List of Field instances.
  std::vector<std::shared_ptr<Field>> fieldsList;

  /// Table of energy contributions.
  std::map<std::string, double> energyContributions;

  /// Total calculation length.
  double calcLength = 0.0;
};

#endif // INTERACTION_H
