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

#ifndef FIELD_COULOMB_EXACT_H
#define FIELD_COULOMB_EXACT_H

/** \file
 *  \brief Headers for the FieldCoulombExact class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "multi.h"

/** \brief This class calculates the Coulomb Exchange Exact field (and some other interaction).
 *
 * This class is derived from the Field class and calculates the Coulomb Exchange Exact field (and some other interaction).
 */

class FieldCoulombExact : public Field
{
public:

  explicit FieldCoulombExact(Field::Parameters fp, State *_state);
  void calcField(void) override;

private:
  void calcDirectField(void);
  void calcIzx(void);
  void calcIrx(void);

  INT nQuad;

  /// Ir exchange. Irx(x,+1 or 0,m,n_alpha,n_gamma)(n_delta + n_beta*nMax(m)).
  FMulti<arma::vec> Irx;

  /// Iz exchange. Izx(x,d_alpha+2*d_delta,n_zalpha,n_zdelta)(n_zgamma,n_zbeta,d_gamma+2*d_beta).
  FMulti<arma::cube> Izx;

  Axis axis_exact;

  void calcIcoul(void);
  void calcJcoul(void);

  //Multi for Coulomb Integral.
  FMulti<arma::vec> Icoul; //Icoul(d_alpha+2*d_gamma,d_delta+2*d_beta,n_d)(n_zd)

  //Multi for Jcouloumb
  FMulti<arma::cube> Jcoul; //Jcoul(n_za,d_alpha+2*d_gamma,mp)(vectorized(n_delta,n_zdelta,d_delta),vectorized(n_beta,n_zbeta,d_beta),n_a)
};

#endif // FIELD_COULOMB_EXACT_H
