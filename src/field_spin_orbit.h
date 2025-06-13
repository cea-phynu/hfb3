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

#ifndef FIELD_SPIN_ORBIT_H
#define FIELD_SPIN_ORBIT_H

/** \file
 *  \brief Headers for the FieldSpinOrbit class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"

/** \brief This class calculates the SpinOrbit field.
 *
 * This class is derived from the Field class and calculates the SpinOrbit field.
 */

class FieldSpinOrbit : public Field
{
public:
  explicit FieldSpinOrbit(Field::Parameters fp, State *_state);
  void calcField(void) override;

  //Previous field, used for testing only (quadrature).
  void calcLegacy(void);

  //Interaction calculations with quadratures. For testing only.
  Multi<arma::mat> fieldLegacy;

  void calcXz(void);
  void calcXr(void);
  void calcYr(void);
  void calcYz(void);

  //============================================================================
  //============================================================================
  //============================================================================

  //Xz for spin-orbit term. Xz(index,n_za,d_alpha+2*d_gamma)(n_zbeta,n_zdelta,d_alpha+2*d_beta).
  Multi<arma::cube> Xz;

  //Xr for spin-orbit term. Xr(index,n_a,mp)(n_beta,n_delta).
  Multi<arma::mat> Xr;

  //Yr for spin-orbit term. Yr(index,m,n_alpha,n_gamma)(n_a).
  Multi<arma::vec> Yr;

  //Yz for spin-orbit term. Yz(index,n_zalpha,n_zgamma,d_alpha+2*d_gamma)(n_za).
  Multi<arma::vec> Yz;

};

#endif // FIELD_SPIN_ORBIT_H
