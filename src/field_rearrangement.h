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

#ifndef FIELD_REARRANGEMENT_H
#define FIELD_REARRANGEMENT_H

/** \file
 *  \brief Headers for the FieldRearrangement class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "discrete.h"
#include "quadratures.h"

/** \brief This class calculates the Rearrangement field.
 *
 * This class is derived from the Field class and calculates the Rearrangement field.
 */

class FieldRearrangement : public Field
{
public:
  explicit FieldRearrangement(Field::Parameters fp, State *_state);
  void calcField(void) override;                                      // #TEST#

private:
  void calcPreZ(void);
  void calcPreR(void);

  /// 4 Discrete instances.
  Discrete discrete;
  Discrete discrete0;

  //Pre-calculations for quadratures.
  Multi<arma::vec> preZ ; //preZ(n_zalpha,n_zbeta,d_alpha+2*d_beta)(z)
  Multi<arma::vec> preR ; //preR(m,n_alpha,n_beta)(r)
};

#endif // FIELD_REARRANGEMENT_H
