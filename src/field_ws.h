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

#ifndef FIELD_WS_H
#define FIELD_WS_H

/** \file
 *  \brief Headers for the FieldWS class.
 */

#include "global.h"
#include "field.h"
#include "wspot.h"

/** \brief This class calculates the Woods-Saxon field.
 *
 * This class is derived from the Field class and calculates the Woods-Saxon field.
 */

class FieldWS : public Field
{
public:

  FieldWS(Field::Parameters fp, State *_state);     // #TEST#
  void calcField(void) override;                                      // #TEST#
  void setDef(const Multi<double> &_def) override;                    // #TEST#
  WSPot getWSPot(void) override;                                      // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  WSPot wspot;

private:

  /// The deformation parameters.
  Multi<double> def;
};

#endif // FIELD_WOODS_SAXON_H
