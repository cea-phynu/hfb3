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

#ifndef FIELD_CENTRAL_H
#define FIELD_CENTRAL_H

/** \file
 *  \brief Headers for the FieldCentral class.
 */

#include "global.h"
#include "field.h"
#include "multi.h"

/** \brief This class calculates the Central Direct field.
 *
 * This class is derived from the Field class and calculates the Central interaction.
 */

class FieldCentral : public Field
{
public:

  explicit FieldCentral(Parameters _parameters, State *_state);
  void calcField(void) override;

private:

  void calcIr(void);
  void calcIz(void);
  void calcJr(void);
  void calcJz(void);

  //============================================================================
  //============================================================================
  //============================================================================

  /// Ir direct. IrDirect(ig,m,n_alpha,n_gamma)(n_delta + n_beta*nMax(m)).
  Multi<arma::vec> IrDirect;

  /// Ir exchange. IrExchange(ig,m,n_alpha,n_gamma)(n_delta + n_beta*nMax(m)).
  Multi<arma::vec> IrExchange;

  /// Iz direct. IzDirect(ig,d_alpha,d_gamma,n_zalpha,n_zgamma)(n_zdelta,n_zbeta,d_delta+2*d_beta).
  FMulti<arma::cube> IzDirect;

  /// Iz exchange. IzExchange(ig,d_alpha,d_delta,n_zalpha,n_zdelta)(n_zgamma,n_zbeta,d_gamma+2*d_beta).
  FMulti<arma::cube> IzExchange;

  /// Jr. Jr(ig,m_alpha,n_alpha,m_beta,n_beta)(n_b).
  Multi<arma::vec> Jr;

  /// Jz. Jz(ig,d_alpha+2*d_beta,d_gamma+2*d_delta,n_zbeta,n_zdelta)(n_z).
  Multi<arma::vec> Jz;

  /// Enumeration for the local field types.
  enum
  {
    MEAN_CENT_DIRECT,
    MEAN_CENT_EXCHANGE_PLUS,
    MEAN_CENT_EXCHANGE_MINUS,
    PAIR_CENT_DIRECT,
    PAIR_CENT_EXCHANGE
  };

  /// Central field parameters.
  double w, b, h, m, p;
};

#endif // FIELD_CENTRAL_H
