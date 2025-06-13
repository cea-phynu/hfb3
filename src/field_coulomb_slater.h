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

#ifndef FIELD_COULOMB_SLATER_H
#define FIELD_COULOMB_SLATER_H

/** \file
 *  \brief Headers for the FieldCoulombSlater class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "multi.h"
#include "discrete.h"

/** \brief This class calculates the Coulomb interaction with the Slater approximation.
 *
 * This class calculates the Coulomb interaction with the Slater approximation.
 */

class FieldCoulombSlater : public Field
{
public:
  explicit FieldCoulombSlater(Field::Parameters fp, State *_stateds);
  void calcField(void) override;
  void calcEnergy(void) override;

private:
  void calcIcoul(void);
  void calcJcoul(void);

  //Multi for Coulomb Integral.
  Multi<arma::vec> Icoul; //Icoul(d_alpha+2*d_gamma,d_delta+2*d_beta,n_d)(n_zd)

  //Multi for Jcouloumb
  Multi<arma::cube>Jcoul; //Jcoul(n_za,d_alpha+2*d_gamma,mp)(vectorized(n_delta,n_zdelta,d_delta),vectorized(n_beta,n_zbeta,d_beta),n_a)

  //2 Discrete instances.
  Discrete discrete;
  Discrete discrete0;
};

#endif // FIELD_COULOMB_SLATER_H
