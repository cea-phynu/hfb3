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

#ifndef FIELD_CDM2_H
#define FIELD_CDM2_H

/** \file
 *  \brief Headers for the FieldCDM2 class.
 */

#include "global.h"
#include "generic.h"
#include "field.h"

/** \brief This class calculates the CDM2 field.
 *
 * This class is derived from the Field class and calculates the CDM2 field.
 */

class FieldCDM2 : public Field
{
public:

  explicit FieldCDM2(Field::Parameters _parameters, State *_state);
  void calcField(void) override;                              // #TEST#

private:

  void calcP(void);
  void calcPblocks(void);
  void calcIz(void);
  void calcIr(void);

  //============================================================================
  //============================================================================
  //============================================================================

  //Cube standing for <alpha_z|nabla_0|beta_z>.
  arma::cube Iz_0;

  //Cube standing for <alpha_z|beta_z>.
  arma::cube Iz_mp;

  //Multi standing for <alpha_r|nabla_+|beta_r>.
  Multi<arma::mat> Ir_plus;

  //Multi for 0 and +- blocks calculations.
  Multi<arma::mat> P;

  //Multi for blocks calculations
  Multi<arma::mat> Pblocks;
  Multi<arma::mat> SizeInit;

  //Matrix for time even/odd blocks.
  arma::mat p;
};

#endif // FIELD_CDM2_H
