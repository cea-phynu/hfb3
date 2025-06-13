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

#ifndef FIELD_REARRANGEMENT_FR_H
#define FIELD_REARRANGEMENT_FR_H

/** \file
 *  \brief Headers for the FieldRearrangementFR class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "discrete.h"

/** \brief This class calculates the RearrangementFR field.
 *
 * This class is derived from the Field class and calculates the RearrangementFR field.
 */

class FieldRearrangementFR : public Field
{
public:
  FieldRearrangementFR(Field::Parameters fp, State *_state);      // #TEST#
  void calcField(void) override;                                      // #TEST#

  void calcPreZ(void);
  void calcPreR(void);
  void calcIz(void);
  void calcIr(void);
  void calcWaveZ(void);
  void calcIndexZ(void);
  void calcWaveMatZ(void);
  void calcRVals(void);

  //Z objects.
  Multi<arma::cube> IzDirect; //IzDirect(z)(n_zdelta,n_zbeta,d_delta+2*d_beta)
  Multi<arma::mat> IzExchange; //IzExchange(d_alpha,n_zalpha,index_z)(n_zgamma,d_gamma)
  Multi<arma::vec> zValsExchange;
  Multi<arma::cube> waveZ;
  Multi<arma::mat> waveMatZ;
  Multi<arma::vec> preZ;

  //R objects.
  Multi<arma::mat> IrDirect; //IrDirect(r,mp)(n_delta,n_beta)
  Multi<arma::vec> IrExchange; //IrDirect(r,m,n_gamma)(VECTOR(mp,n_beta))
  Multi<arma::vec> IrExchange_3; //IrDirect(r,m,n_gamma)(VECTOR(mp,n_beta))
  Multi<arma::vec> IrExchange_1_p1; //IrDirect(r,m,n_gamma)(VECTOR(mp,n_beta))
  Multi<arma::vec> rVals;
  Multi<arma::vec> rValsExchange;
  Multi<arma::vec> preR;

  //2 Discrete instances.
  Discrete discrete;
  Discrete discrete0;

  IVEC index_M_0;
  IVEC index_M_1_p1;
  IVEC index_M_1_p2;
  INT index_mn_0;
  INT index_mn_1_p1;
  INT index_mn_1_p2;

  INT fullIndex_pp;
  INT fullIndex_mm;
  INT fullIndex_mp;
  IVEC index_Z_pp;
  IVEC index_Z_mm;
  IVEC index_Z_mp;
};

#endif // FIELD_REARRANGEMENT_FR_H
