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

#ifndef FIELD_DENSITY_FR_H
#define FIELD_DENSITY_FR_H

/** \file
 *  \brief Headers for the FieldDensityFR class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"
#include "discrete.h"

/** \brief This class calculates the DensityFR field.
 *
 * This class is derived from the Field class and calculates the DensityF2 field.
 */

class FieldDensityFR : public Field
{
public:
  FieldDensityFR(Field::Parameters fp, State *_state);    // #TEST#
  void calcField(void) override;                                      // #TEST#

  void calcIz(void);
  void calcPreZ(void);
  void calcWaveZ(void);
  void calcWaveMatZ(void);
  void calcIr(void);
  void calcPreR(void);
  void calcRVals(void);

  //z quantities.
  FMulti<arma::cube> IzDirect;
  FMulti<arma::mat>  IzExchange;
  Multi<arma::vec>  newIzDirect;
  Multi<arma::vec>  zValsExchange;
  Multi<arma::cube> waveZ;
  Multi<arma::mat>  waveMatZ;
  Multi<arma::vec>  preZ;

  //R quantities.
  Multi<arma::mat> IrDirect;
  Multi<arma::vec> IrExchange;
  Multi<arma::vec> IrExchange_3;
  Multi<arma::vec> IrExchange_1_p1;
  Multi<arma::vec> IrExchange_1_p2;
  Multi<arma::vec> newIrDirect;
  FMulti<arma::vec> rVals;
  FMulti<arma::vec> rValsExchange;
  FMulti<arma::vec> rValsExchange2;
  Multi<arma::vec> preR;
  IVEC index_M_0;
  IVEC index_M_1_p1;
  IVEC index_M_1_p2;
  INT index_mn_0;
  INT index_mn_1_p1;
  INT index_mn_1_p2;


  //Discrete Instances.
  Discrete discrete;
  Discrete discrete0;

private:
  // big allocations

  FMulti<arma::cube> RmatDirect;
  FMulti<arma::mat> RmatExchange_0;
  FMulti<arma::mat> Kmat_0;
  FMulti<arma::mat> RmatExchange_1_p1;
  FMulti<arma::mat> RmatExchange_1_p2;
  FMulti<arma::mat> Kmat_1_p1;
  FMulti<arma::mat> Kmat_1_p2;

  FMulti<arma::cube> Rzze_0;
  FMulti<arma::cube> Kzze_0;

  FMulti<arma::cube> Rzze_1_p1;
  FMulti<arma::cube> Kzze_1_p1;
  FMulti<arma::cube> Rzze_1_p2;
  FMulti<arma::cube> Kzze_1_p2;

  FMulti<arma::mat> Rfinal_0;
  FMulti<arma::mat> Kfinal_0;
  FMulti<arma::mat> Rfinal_3;
  FMulti<arma::mat> Kfinal_3;

  FMulti<arma::mat> Rfinal_1_p1;
  FMulti<arma::mat> Kfinal_1_p1;
  FMulti<arma::mat> Rfinal_1_p2;
  FMulti<arma::mat> Kfinal_1_p2;

  FMulti<arma::vec> Rz;
  FMulti<arma::vec> Rz2;

  FMulti<arma::vec> Rz_0;
  FMulti<arma::vec> Kz_0;
  FMulti<arma::vec> Rz_3;
  FMulti<arma::vec> Kz_3;

  FMulti<arma::vec> Rz_1_p1;
  FMulti<arma::vec> Kz_1_p1;
  FMulti<arma::vec> Rz_1_p2;
  FMulti<arma::vec> Kz_1_p2;


  FMulti<arma::cube> Gamma;

  FMulti<arma::cube> GammaE_0;
  FMulti<arma::cube> Delta_0;
  FMulti<arma::cube> GammaE_3;
  FMulti<arma::cube> Delta_3;


  FMulti<arma::cube> GammaE_1;
  FMulti<arma::cube> Delta_1;

};

#endif // FIELD_DENSITY_FR_H
