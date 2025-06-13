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

#include "geometrical_operators.h"
#include "global.h"
#include "quadratures.h"
#include "tools.h"
#include "geometry.h"
#include "discrete.h"

/** \file
 *  \brief Methods of the GeometricalOperators class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

GeometricalOperators::GeometricalOperators(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

GeometricalOperators::GeometricalOperators(State _state) : state(_state), multipoleOperators(_state)
{
  DBG_ENTER;

  Basis &basis = state.basis;
  basis.calcWDN();

  Discrete discrete(&state.basis, Mesh::regular(0, 0, -20, 20, 0, 20, 201, 1, 401));
  arma::mat densn = discrete.getLocalXZ(state.rho(NEUTRON), true);
  arma::mat densp = discrete.getLocalXZ(state.rho(PROTON ), true);
  arma::mat denst = densn + densp;
  Geometry geomt(discrete.mesh, denst, state.sys);

  izNeck = geomt.izNeck;

  if (izNeck != -1)
  {
    zNeck  = geomt.neckPos;
  }
  else
  {
    zNeck = 0.0;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the geometrical operator \f$\hat{Q}_{N}\f$.
 *
 *  For a definition of the \f$\hat{Q}_{N}\f$ multipole moment operators, cf. \ref
 * qneck.
 */

void GeometricalOperators::calcQNeckMatrix(void)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  calcQNk();  // Calculate the qNeck matrix elements in the 2ct HO basis.

  if (!qneck0.empty())
  {
    DBG_LEAVE;
  }

  matQNk = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  matS0 = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  for (UINT ia = 0; ia < basis.HOqn.nb; ia++)
  {
    INT ma   = basis.HOqn(0, ia);
    INT n_za = basis.HOqn(2, ia);
    INT da   = basis.HOqn(3, ia);
    INT sa   = basis.HOqn(4, ia);

    for (UINT ib = 0; ib < basis.HOqn.nb; ib++)
    {
      INT mb   = basis.HOqn(0, ib);
      INT n_zb = basis.HOqn(2, ib);
      INT db   = basis.HOqn(3, ib);
      INT sb   = basis.HOqn(4, ib);

      if (ma != mb) continue;

      if (sa == sb) matS0(ia, ib) = 1.0;

      matQNk(ia, ib) = matQN(n_za, da)(n_zb, db);

    }     // ib
  }       // ia

  qneck0(0) = matS0 % matQNk % multipoleOperators.matR2k(0);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

// Calculate the \f$\langle n_z^a, d_a|\hat{Q}_{N}|n_z^b, d_b\rangle\f$ values.

void GeometricalOperators::calcQNk()
{
  DBG_ENTER;

  Basis &basis = state.basis;
  basis.calcTalmanz();

  double aN = 1.0;
  double fac = 1.0 / sqrt(pow(aN, 2) + pow(basis.b_z, 2));

  double sum0, dab, fac_exp;
  INT    n_zcMin;

  for (INT n_za = 0; n_za < basis.n_zGlobalMax; n_za++)
  {
    for (INT da = 0; da < basis.dMax; da++)
    {
      matQN(n_za, da) = arma::zeros(basis.n_zGlobalMax, basis.dMax);
      arma::mat &mat    = matQN(n_za, da);

      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax; n_zb++)
      {
        for (INT db = 0; db < basis.dMax; db++)
        {
          dab     = (1 - da - db) * basis.d_0 / 2.0;
          fac_exp = exp(- pow( ( zNeck - dab ) * fac, 2 ) );
          arma::vec &talvec = basis.talmanz(n_za, da, n_zb, db);

          n_zcMin = 0;
          if (da == db) n_zcMin = abs(n_za - n_zb);

          sum0     = 0.0;

          for (INT n_zc = n_zcMin; n_zc < n_za + n_zb + 1; n_zc++) sum0 += talvec[n_zc] * mnz[n_zc]
                * pow(basis.b_z / 2.0 * fac, n_zc)
                * basis.hermite(n_zc, (zNeck - dab) * fac);

          mat(n_zb, db) = sum0 * fac * aN * fac_exp;
        }     // db
      }       // n_zb
    }         // da
  }           // n_za

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the geometrical operator \f$\langle
 * m,n_a,a,d_a|\hat{Q}_{N}|m,n_b,b,d_b\rangle\f$.
 *
 *  For a definition of \f$\langle m,n_a,a,d_a|\hat{Q}_{N}|m,n_b,b,d_b\rangle\f$
 *  cf. \ref qneck.
 */

void GeometricalOperators::calcQneck(const arma::mat &rhon, const arma::mat &rhop)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  if (rhon.empty() || rhop.empty())
  {
    DBG_LEAVE;
  }

  calcQNeckMatrix();
  basis.calcWDN();

  qneck = arma::zeros(1, 3);

  arma::mat rn    = 2 * rhon; // The coefficient 2 takes into account time reversal negative states.
  arma::mat rp    = 2 * rhop;
  arma::mat ra    = rn + rp;

  qneck(0, NEUTRON) = arma::accu(qneck0(0) % rn);
  qneck(0, PROTON)  = arma::accu(qneck0(0) % rp);
  qneck(0, TOTAL)   = arma::accu(qneck0(0) % ra);

  DBG_LEAVE;
}

