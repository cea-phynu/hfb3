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

#ifndef FIELD_DENSITY_H
#define FIELD_DENSITY_H

/** \file
 *  \brief Headers for the FieldDensity class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "discrete.h"
#include "state.h"

/** \brief This class calculates the Density field.
 *
 * This class is derived from the Field class and calculates the Density field.
 */

class FieldDensity : public Field
{
public:
  FieldDensity(Field::Parameters fp, State *_state);      // #TEST#
  void calcField(void) override;                                      // #TEST#

private:
  void calcPreZ(void);
  void calcPreR(void);

  /// 4 Discrete instances.
  Discrete discrete;
  Discrete discrete0;
  Discrete discrete1;
  Discrete discrete2;

  //Pre-calculations for quadratures.
  Multi<arma::vec> preZ ; //preZ(n_zalpha,n_zbeta,d_alpha+2*d_beta)(z)
  Multi<arma::vec> preR ; //preR(m,n_alpha,n_beta)(r)
};

#endif // FIELD_DENSITY_H
