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

#ifndef FIELD_KINETIC_H
#define FIELD_KINETIC_H

/** \file
 *  \brief Headers for the FieldKinetic class.
 */

#include "global.h"
#include "field.h"

/** \brief This class calculates the Kinetic field.
 *
 * This class is derived from the Field class and calculates the Kinetic field.
 */

class FieldKinetic : public Field
{
public:
  explicit FieldKinetic(Parameters _parameters, State *_state);
  void calcField(void) override;                                      // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;
};

#endif // FIELD_KINETIC_H
