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

#include "multipole_operators.h"
#include "global.h"
#include "quadratures.h"
#include "tools.h"

/** \file
 *  \brief Methods of the MultipoleOperators class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

MultipoleOperators::MultipoleOperators(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

MultipoleOperators::MultipoleOperators(State _state) : state(_state)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  if ((general.compatibility != General::COMPAT_NONE) && (!general.compatWarningDisplayed))
  {
    general.compatWarningDisplayed = true;
    Tools::warning("Using non-standard Multipole Moment Operators ! (general.compatibility mode set)");
  }

  // dependencies
  basis.calcWDN();
  calcZk();
  calcRk();
  calcQlmObs(state.rho);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the multipole moment operators \f$\hat{Q}_{\lambda,0}\f$.
 *
 *  For a definition of the \f$\hat{Q}_{\lambda,0}\f$ multipole moment operators, cf. \ref
 * multipol_cyl.
 */

void MultipoleOperators::calcQlmHO(void)
{
  DBG_ENTER;

  calcRkn();

  Basis &basis = state.basis;

  if (!qlmHO.empty())
  {
    DBG_LEAVE;
  }

  for (INT k = 0; k < 7; k++)
  {
    matZk(k) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  }

  for (INT k = 0; k < 4; k++)
  {
    matR2k(k) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  }

  matS0 = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  arma::mat matDeltaOmega0 = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  arma::mat matDeltaOmega1 = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  arma::mat matDeltaOmega2 = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  for (UINT ia = 0; ia < basis.HOqn.nb; ia++)
  {
    INT ma   = basis.HOqn(0, ia);
    INT na   = basis.HOqn(1, ia);
    INT n_za = basis.HOqn(2, ia);
    INT da   = basis.HOqn(3, ia);
    INT sa   = basis.HOqn(4, ia);

    for (UINT ib = 0; ib < basis.HOqn.nb; ib++)
    {
      INT mb   = basis.HOqn(0, ib);
      INT nb   = basis.HOqn(1, ib);
      INT n_zb = basis.HOqn(2, ib);
      INT db   = basis.HOqn(3, ib);
      INT sb   = basis.HOqn(4, ib);

      if (ma == mb    ) matDeltaOmega0(ia, ib) = 1.0;
      if (ma == mb + 1) matDeltaOmega1(ia, ib) = 1.0;
      if (ma == mb + 2) matDeltaOmega2(ia, ib) = 1.0;

      if (sa == sb) matS0(ia, ib) = 1.0;

      for (INT k = 0; k < 7; k++)
      {
        // TODO: optimize me
        matZk(k)(ia, ib) = matZ(k, n_za, da)(n_zb, db);
      }

      if (ma == mb)
      {
        for (INT k = 0; k < 4; k++)
        {
          matR2k(k)(ia, ib) = matR(k, ma, na)(nb);
        }
      }
    }     // ib
  }       // ia

  switch(general.compatibility)
  {
    case General::COMPAT_NONE:
      // see https://en.wikipedia.org/wiki/Table_of_spherical_harmonics for a list of some YLM functions
      // see https://en.wikipedia.org/wiki/List_of_common_coordinate_transformations for the common coordinate transformations

      // QLM = rho * r^l * YLM (see Ring, Schuck)
      qlmHO(0,  0) = 1. / 2.  * sqrt( 1. / PI) * matS0 % matDeltaOmega0 % ( 1 * matZk(0) % matR2k(0));
      qlmHO(1,  0) = 1. / 2.  * sqrt( 3. / PI) * matS0 % matDeltaOmega0 % ( 1 * matZk(1) % matR2k(0));
      qlmHO(2,  0) = 1. / 4.  * sqrt( 5. / PI) * matS0 % matDeltaOmega0 % ( 2 * matZk(2) % matR2k(0) -   1 * matR2k(1) % matZk(0));
      qlmHO(3,  0) = 1. / 4.  * sqrt( 7. / PI) * matS0 % matDeltaOmega0 % ( 2 * matZk(3) % matR2k(0) -   3 * matR2k(1) % matZk(1));
      qlmHO(4,  0) = 1. / 16. * sqrt( 9. / PI) * matS0 % matDeltaOmega0 % ( 8 * matZk(4) % matR2k(0) -  24 * matR2k(1) % matZk(2) +  3 * matR2k(2) % matZk(0));
      qlmHO(5,  0) = 1. / 16. * sqrt(11. / PI) * matS0 % matDeltaOmega0 % ( 8 * matZk(5) % matR2k(0) -  40 * matR2k(1) % matZk(3) + 15 * matR2k(2) % matZk(1));
      qlmHO(6,  0) = 1. / 32. * sqrt(13. / PI) * matS0 % matDeltaOmega0 % (16 * matZk(6) % matR2k(0) - 120 * matR2k(1) % matZk(4) + 90 * matR2k(2) % matZk(2) - 5 * matR2k(3) % matZk(0));

      break;
    case General::COMPAT_BERGER:
      // QLM = sqrt( (4 * PI) / (2 * L + 1) ) * rho * r^l * YLM (see, e.g., Messiah, Albert (1999) "Quantum mechanics : two volumes bound as one")
      // Berger added a factor 2 in the Q20 definition !!!
      qlmHO(0, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(0) % matR2k(0));
      qlmHO(1, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(1) % matR2k(0));
      qlmHO(2, 0) = matS0 % matDeltaOmega0 % (2.0 * matZk(2) % matR2k(0) - 1.0          * matR2k(1) % matZk(0));
      qlmHO(3, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(3) % matR2k(0) - 3.0   / 2.0  * matR2k(1) % matZk(1));
      qlmHO(4, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(4) % matR2k(0) - 24.0  / 8.0  * matR2k(1) % matZk(2) + 3.0  / 8.0  * matR2k(2) % matZk(0));
      qlmHO(5, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(5) % matR2k(0) - 40.0  / 8.0  * matR2k(1) % matZk(3) + 15.0 / 8.0  * matR2k(2) % matZk(1));
      qlmHO(6, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(6) % matR2k(0) - 120.0 / 16.0 * matR2k(1) % matZk(4) + 90.0 / 16.0 * matR2k(2) % matZk(2) - 5.0 / 16.0 * matR2k(3) % matZk(0));
      break;
    case General::COMPAT_ROBLEDO:
      // QLM = sqrt( (4 * PI) / (2 * L + 1) ) * rho * r^l * YLM (see, e.g., Messiah, Albert (1999) "Quantum mechanics : two volumes bound as one" or wikipedia)
      qlmHO(0, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(0) % matR2k(0));
      qlmHO(1, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(1) % matR2k(0));
      qlmHO(2, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(2) % matR2k(0) - 0.5          * matR2k(1) % matZk(0));
      qlmHO(3, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(3) % matR2k(0) - 3.0   / 2.0  * matR2k(1) % matZk(1));
      qlmHO(4, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(4) % matR2k(0) - 24.0  / 8.0  * matR2k(1) % matZk(2) + 3.0  / 8.0  * matR2k(2) % matZk(0));
      qlmHO(5, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(5) % matR2k(0) - 40.0  / 8.0  * matR2k(1) % matZk(3) + 15.0 / 8.0  * matR2k(2) % matZk(1));
      qlmHO(6, 0) = matS0 % matDeltaOmega0 % (1.0 * matZk(6) % matR2k(0) - 120.0 / 16.0 * matR2k(1) % matZk(4) + 90.0 / 16.0 * matR2k(2) % matZk(2) - 5.0 / 16.0 * matR2k(3) % matZk(0));
      break;
    case General::COMPAT_HFBTHO:
      /* HFBTHO uses the definitions of General::COMPAT_ROBLEDO for Q00,Q10,
      *  the definition of General::COMPAT_BERGER for Q20,
      *  and the definitions General::COMPAT_NONE for Q30 and higher.
      *  No, it doesn't make any sense.
      */
      qlmHO(0, 0) =                             matS0 % matDeltaOmega0 % ( 1 * matZk(0) % matR2k(0));
      qlmHO(1, 0) =                             matS0 % matDeltaOmega0 % ( 1 * matZk(1) % matR2k(0));
      qlmHO(2, 0) =                             matS0 % matDeltaOmega0 % ( 2 * matZk(2) % matR2k(0) -   1 * matR2k(1) % matZk(0));
      qlmHO(3, 0) = 1. / 4.  * sqrt( 7. / PI) * matS0 % matDeltaOmega0 % ( 2 * matZk(3) % matR2k(0) -   3 * matR2k(1) % matZk(1));
      qlmHO(4, 0) = 1. / 16. * sqrt( 9. / PI) * matS0 % matDeltaOmega0 % ( 8 * matZk(4) % matR2k(0) -  24 * matR2k(1) % matZk(2) +  3 * matR2k(2) % matZk(0));
      qlmHO(5, 0) = 1. / 16. * sqrt(11. / PI) * matS0 % matDeltaOmega0 % ( 8 * matZk(5) % matR2k(0) -  40 * matR2k(1) % matZk(3) + 15 * matR2k(2) % matZk(1));
      qlmHO(6, 0) = 1. / 32. * sqrt(13. / PI) * matS0 % matDeltaOmega0 % (16 * matZk(6) % matR2k(0) - 120 * matR2k(1) % matZk(4) + 90 * matR2k(2) % matZk(2) - 5 * matR2k(3) % matZk(0));
      break;
    default:
      ERROR("Unknown general.compatibility value.");
  }

  // TODO: check that this is still true for non-default compatibility values
  qlmHO(2, -2) =  1. / 4. * sqrt(15. / (2.0 * PI)) * (2.0 * PI) * matS0 % matDeltaOmega2.t() % matZk(0) % matRkn(3);
  qlmHO(2, -1) =  1. / 2. * sqrt(15. / (2.0 * PI)) * (2.0 * PI) * matS0 % matDeltaOmega1.t() % matZk(1) % matRkn(2);
  qlmHO(2,  1) = -1. / 2. * sqrt(15. / (2.0 * PI)) * (2.0 * PI) * matS0 % matDeltaOmega1     % matZk(1) % matRkn(2);
  qlmHO(2,  2) =  1. / 4. * sqrt(15. / (2.0 * PI)) * (2.0 * PI) * matS0 % matDeltaOmega2     % matZk(0) % matRkn(3);


  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\langle m,n_a|\hat{r}^{k}|m,n_b\rangle\f$ values.
 */

void MultipoleOperators::calcRkn(void)
{
  Multi<MAT> rkn;

  Axis axis(Axis::GAUSS_LAGUERRE, 60);

  for (INT k = 2; k < 4; k++)
  {
    VEC crkw = pow(state.basis.b_r, k - 1) / (2.0 * PI) * arma::pow(arma::sqrt(axis.p), k - 1) % axis.w;

    for (INT ma = 0; ma < state.basis.mMax; ma++)
    {
      for (INT mb = ma; mb < state.basis.mMax; mb++)
      {
        rkn(k, ma, mb) = arma::zeros(state.basis.nMax(ma), state.basis.nMax(mb));

        for (INT na = 0; na < state.basis.nMax(ma); na++)
        {
          VEC Ra = state.basis.rPartNormReduced(axis.p, ma, na);

          for (INT nb = 0; nb < state.basis.nMax(mb); nb++)
          {
            VEC Rb = state.basis.rPartNormReduced(axis.p, mb, nb);

            VEC func = Ra % Rb % crkw;
            rkn(k, ma, mb)(na, nb) = arma::accu(func);
          }
        }
        rkn(k, mb, ma) = rkn(k, ma, mb).t();
      }
    }

    matRkn(k) = arma::zeros(state.basis.HOqn.nb, state.basis.HOqn.nb);

    for (UINT ia = 0; ia < state.basis.HOqn.nb; ia++)
    {
      INT ma   = state.basis.HOqn(0, ia);
      INT na   = state.basis.HOqn(1, ia);

      for (UINT ib = ia; ib < state.basis.HOqn.nb; ib++)
      {
        INT mb   = state.basis.HOqn(0, ib);
        INT nb   = state.basis.HOqn(1, ib);

        matRkn(k)(ia, ib) = rkn(k, ma, mb)(na, nb);
        matRkn(k)(ib, ia) = matRkn(k)(ia, ib);
      }     // ib
    }       // ia
  }         // k
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the equivalent \f$\beta\f$ value from a \f$q_{20}\f$ value.
 */

double MultipoleOperators::getBetaFromQ20(double parts, double A, double vq20)
{
  DBG_ENTER;

  switch(general.compatibility)
  {
    case General::COMPAT_NONE:
      DBG_RETURN(            4.0 * PI / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2)) * vq20);
      break;
    case General::COMPAT_BERGER:
      DBG_RETURN(      sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2)) * vq20);
      break;
    case General::COMPAT_ROBLEDO:
      DBG_RETURN(2.0 * sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2)) * vq20);
      break;
    case General::COMPAT_HFBTHO:
      DBG_RETURN(      sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2)) * vq20);
      break;
    default:
      ERROR("Unknown general.compatibility value.");
  }

  // should not happen
  DBG_RETURN(-1.0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the equivalent \f$q_{20}\f$ value from a \f$\beta\f$ value.
 */

double MultipoleOperators::getQ20FromBeta(double parts, double A, double vBeta)
{
  DBG_ENTER;

  switch(general.compatibility)
  {
    case General::COMPAT_NONE:
      DBG_RETURN(vBeta / (            4.0 * PI / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2))));
      break;
    case General::COMPAT_BERGER:
      DBG_RETURN(vBeta / (      sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2))));
      break;
    case General::COMPAT_ROBLEDO:
      DBG_RETURN(vBeta / (2.0 * sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2))));
      break;
    case General::COMPAT_HFBTHO:
      DBG_RETURN(vBeta / (      sqrt(5.0 * PI) / (3.0 * parts * pow( R_0 * pow(A, 1.0 / 3.0), 2))));
      break;
    default:
      ERROR("Unknown general.compatibility value.");
  }

  // should not happen
  DBG_RETURN(-1.0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\beta\f$ values.
 */

void MultipoleOperators::calcBeta(void)
{
  DBG_ENTER;

  beta = arma::zeros(3);

  double Z = state.sys.nProt;
  double N = state.sys.nNeut;
  double A = Z + N;

  beta(NEUTRON) = getBetaFromQ20(N, A, qlmObs(2, 0, NEUTRON));
  beta(PROTON ) = getBetaFromQ20(Z, A, qlmObs(2, 0, PROTON ));
  beta(TOTAL  ) = getBetaFromQ20(A, A, qlmObs(2, 0, TOTAL  ));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\langle m,n_a|\hat{r}^{2l}|m,n_b\rangle\f$ values.
 * TODO: use A-28 in Berger's PhD + Talmanr transformation ?
 */

void MultipoleOperators::calcRk(void)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  // dependencies
  basis.calcTalmanr();

  if (!matR.empty())
  {
    DBG_LEAVE;
  }

  INT lmax = 4;

  for (INT l = 0; l < lmax; l++)
  {
    for (INT ma = 0; ma < basis.mMax; ma++)
    {
      for (INT na = 0; na < basis.nMax(ma); na++)
      {
        matR(l, ma, na) = arma::zeros(basis.nMax(ma));
        arma::vec &vec  = matR(l, ma, na);

        for (INT nb = 0; nb < basis.nMax(ma); nb++)
        {
          double     factor = pow(fact[l], 2) * pow(basis.b_r, 2 * l);
          double     sum    = 0.0;
          arma::vec &vec2   = basis.talmanr(ma, na, ma, nb);

          for (INT nc = abs(na - nb); nc <= MIN(l, na + nb + abs(ma)); nc++)
          {
            sum += pow(-1.0, nc) / (fact[nc] * fact[l - nc]) * vec2(nc);
          }

          vec(nb) = sum * factor;
        }     // nb
      }       // na
    }         // ma
  }           // l

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\langle n_z^a, d_a|\hat{z}^{k}|n_z^b, d_b\rangle\f$ values.
 * TODO: use A-27 in Berger's PhD + Talmanz transformation ?
 */

void MultipoleOperators::calcZk(void)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  // dependencies
  basis.calcTalmanz();

  if (!matZ.empty())
  {
    DBG_LEAVE;
  }

  INT    kmax = 7; // <- WARNING hard-coded value
  double sum0, sum1, factor, dab;
  INT    n_zcMin;
  double tempv;

  for (INT k = 0; k < kmax; k++)
  {
    for (INT n_za = 0; n_za < basis.n_zGlobalMax; n_za++)
    {
      for (INT da = 0; da < basis.dMax; da++)
      {
        matZ(k, n_za, da) = arma::zeros(basis.n_zGlobalMax, basis.dMax);
        arma::mat &mat    = matZ(k, n_za, da);

        for (INT n_zb = 0; n_zb < basis.n_zGlobalMax; n_zb++)
        {
          for (INT db = 0; db < basis.dMax; db++)
          {
            dab     = (1.0 - da - db) * basis.d_0 / 2.0;
            factor  = fact[k];
            sum0    = 0.0;
            n_zcMin = 0;

            if (da == db) n_zcMin = abs(n_za - n_zb);

            arma::vec &talvec = basis.talmanz(n_za, da, n_zb, db);

            for (INT n_zc = n_zcMin; n_zc < n_za + n_zb + 1; n_zc++)
            {
              sum1 = 0.0;

              for (INT l = n_zc; l <= k; l += 2)
              {
                tempv = pow(basis.b_z / 2.0, l) * pow(dab, k - l) /
                        (fact[(l - n_zc) / 2] * fact[k - l]);
                sum1 += tempv;
              }

              tempv = sum1 * talvec(n_zc) * mnz[n_zc];
              sum0 += tempv;
            }

            mat(n_zb, db) = sum0 * factor;
          }     // db
        }       // n_zb
      }         // da
    }           // n_za
  }             // k

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the multipole moments \f$\langle
 * m,n_a,a,d_a|\hat{Q}_{\lambda,0}|m,n_b,b,d_b\rangle\f$.
 *
 *  For a definition of the \f$\langle m,n_a,a,d_a|\hat{Q}_{\lambda,0}|m,n_b,b,d_b\rangle\f$
 * multipole moments, cf. \ref multipol_cyl.
 */

void MultipoleOperators::calcQlmObs(const Multi<arma::mat> &rho)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  if (rho(NEUTRON).empty() || rho(PROTON).empty())
  {
    DBG_LEAVE;
  }

  // dependencies
  calcQlmHO();
  basis.calcWDN();

  qlmObs.clear();

  arma::mat rn    = 2 * rho(NEUTRON);     // coeff 2 to take into account time reversal negative states
  arma::mat rp    = 2 * rho(PROTON );     // same
  arma::mat rt    = rn + rp;
  // arma::mat s     = basis.ORtoHO * basis.ORtoHO.t();

  qlmObs(0, 0, NEUTRON) = arma::accu(qlmHO(0, 0) % rn);
  qlmObs(0, 0, PROTON ) = arma::accu(qlmHO(0, 0) % rp);
  qlmObs(0, 0, TOTAL  ) = arma::accu(qlmHO(0, 0) % rt);
  qlmObs(1, 0, NEUTRON) = arma::accu(qlmHO(1, 0) % rn);
  qlmObs(1, 0, PROTON ) = arma::accu(qlmHO(1, 0) % rp);
  qlmObs(1, 0, TOTAL  ) = arma::accu(qlmHO(1, 0) % rt);
  qlmObs(2, 0, NEUTRON) = arma::accu(qlmHO(2, 0) % rn);
  qlmObs(2, 0, PROTON ) = arma::accu(qlmHO(2, 0) % rp);
  qlmObs(2, 0, TOTAL  ) = arma::accu(qlmHO(2, 0) % rt);
  qlmObs(2, 1, NEUTRON) = arma::accu(qlmHO(2, 1) % rn);
  qlmObs(2, 1, PROTON ) = arma::accu(qlmHO(2, 1) % rp);
  qlmObs(2, 1, TOTAL  ) = arma::accu(qlmHO(2, 1) % rt);
  qlmObs(2, 2, NEUTRON) = arma::accu(qlmHO(2, 2) % rn);
  qlmObs(2, 2, PROTON ) = arma::accu(qlmHO(2, 2) % rp);
  qlmObs(2, 2, TOTAL  ) = arma::accu(qlmHO(2, 2) % rt);
  qlmObs(3, 0, NEUTRON) = arma::accu(qlmHO(3, 0) % rn);
  qlmObs(3, 0, PROTON ) = arma::accu(qlmHO(3, 0) % rp);
  qlmObs(3, 0, TOTAL  ) = arma::accu(qlmHO(3, 0) % rt);
  qlmObs(4, 0, NEUTRON) = arma::accu(qlmHO(4, 0) % rn);
  qlmObs(4, 0, PROTON ) = arma::accu(qlmHO(4, 0) % rp);
  qlmObs(4, 0, TOTAL  ) = arma::accu(qlmHO(4, 0) % rt);
  qlmObs(5, 0, NEUTRON) = arma::accu(qlmHO(5, 0) % rn);
  qlmObs(5, 0, PROTON ) = arma::accu(qlmHO(5, 0) % rp);
  qlmObs(5, 0, TOTAL  ) = arma::accu(qlmHO(5, 0) % rt);
  qlmObs(6, 0, NEUTRON) = arma::accu(qlmHO(6, 0) % rn);
  qlmObs(6, 0, PROTON ) = arma::accu(qlmHO(6, 0) % rp);
  qlmObs(6, 0, TOTAL  ) = arma::accu(qlmHO(6, 0) % rt);

  calcBeta();
  calcNPart(rho);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the number of particles from Q00 moments.
 */

void MultipoleOperators::calcNPart(const Multi<arma::mat> &rho)
{
  DBG_ENTER;

  if (rho(NEUTRON).empty() || rho(PROTON).empty())
  {
    DBG_LEAVE;
  }

  nPart = arma::zeros(3);

  double correctionFactor = 1. / 2.  * sqrt( 1. / PI);

  if (general.compatibility == General::COMPAT_BERGER) correctionFactor = 1.0;
  if (general.compatibility == General::COMPAT_HFBTHO) correctionFactor = 1.0;

  for (INT iso = 0; iso < 3; iso++)
  {
    nPart(iso) = qlmObs(0, 0, iso) / correctionFactor;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice table with the multipole moment mean values.
 */

const std::string MultipoleOperators::getNiceInfo(void)
{
  DBG_ENTER;

  // ASSERT(state.checkSolution(), "state does not contain a solution.");

  // color constrained values
  Multi<std::string> colorStr;

  for (auto &c : state.constraints) colorStr(c.second.lambda, c.second.mu, c.second.iso) = TABLE_BLUE;

  std::list<std::list<std::string> > tempList;

  tempList.push_back(
  {
    "Mass    [n]",
    PF("%.6f", nPart(NEUTRON)),
    PF("%.6f", nPart(PROTON )),
    PF("%.6f", nPart(TOTAL  ))
  });

  std::map< std::string, std::vector<double> > tempMap;

  for (INT lambda = 0; lambda < 7; lambda++)
  {
    INT mu = 0;

    tempList.push_back(
    {
      PF("<Q%1d%d> [fm%1d]", lambda, mu, lambda),
      PF("%s%.6f", colorStr(lambda, mu, NEUTRON).c_str(), qlmObs(lambda, mu, NEUTRON)),
      PF("%s%.6f", colorStr(lambda, mu, PROTON ).c_str(), qlmObs(lambda, mu, PROTON )),
      PF("%s%.6f", colorStr(lambda, mu, TOTAL  ).c_str(), qlmObs(lambda, mu, TOTAL  ))
    });
  }

  std::string result =
    Tools::valueTable("Deformations", {"Neutron", "Proton", "Neut+Prot"}, {"", "", ""}, tempList);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string MultipoleOperators::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  std::vector<std::pair<std::string, std::string>> list;

  if (isShort)
  {
    if (qlmObs.empty())
    {
      DBG_RETURN(PF_RED("empty"));
    }

    auto keys = qlmObs.getKeys();

    for (auto key: keys)
    {
      int lambda = key[0];
      int mu     = key[1];
      int iso    = key[2];
      if (iso == TOTAL) list.push_back({PF("Q%01d%01d", lambda, mu), PF("%.6e", qlmObs(lambda, mu, TOTAL))});
    }

    list.push_back({"beta", PF("%.6e", beta(TOTAL))});
  }
  else
  {
    list.push_back({"MultipoleOperators", ""});

    if (qlmObs.empty())
    {
      list.push_back({"qlm", PF_RED("empty")});
    }
    else
    {
      for (auto key: qlmObs.getKeys())
      {
        int lambda = key[0];
        int mu     = key[1];
        int iso    = key[2];
        list.push_back({PF("<Q%01d%01d> ", lambda, mu), PF_YELLOW(Tools::strIsospin(iso)) + ":" + PF_BLUE("%13.6e", qlmObs(lambda, mu, iso))});
      }

      list.push_back({"beta  ", "(" + PF_YELLOW("neut.") + ":" + PF_BLUE("%13.6e", beta(NEUTRON)) +
                      ", " + PF_YELLOW("prot.") + ":" +
                      PF_BLUE("%13.6e", beta(PROTON)) + ", " + PF_YELLOW("total") +
                      ":" + PF_BLUE("%13.6e", beta(TOTAL)) + ")"});
    }
  }

  result += Tools::treeStr(list, isShort);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Save the results in a DataTree instance.
 */

DataTree MultipoleOperators::getDataTree(void)
{
  DBG_ENTER;

  DataTree dt;

  dt.set("multipoleOperators/qlmObs", qlmObs);
  dt.set("multipoleOperators/beta"  , beta);

  DBG_RETURN(dt);
}
