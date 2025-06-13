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

  if (general.compatibility != General::COMPAT_NONE)
  {
    Tools::warning("Using non-standard Multipole Moment Operators ! (general.compatibility mode set)");
  }

  // dependencies
  basis.calcWDN();
  calcZk();
  calcR2l();
  calcQlm(state.rho);

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

void MultipoleOperators::calcQl0Matrices(void)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  if (!ql0.empty())
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

      if (ma != mb) continue;

      if (sa == sb) matS0(ia, ib) = 1.0;

      for (INT k = 0; k < 7; k++)
      {
        // TODO: optimize me
        matZk(k)(ia, ib) = matZ(k, n_za, da)(n_zb, db);
      }

      for (INT k = 0; k < 4; k++)
      {
        // TODO: optimize me
        matR2k(k)(ia, ib) = matR(k, ma, na)(nb);
      }
    }     // ib
  }       // ia

  switch(general.compatibility)
  {
    case General::COMPAT_NONE:
      // QLM = rho * r^l * YLM (see Ring, Schuck)
      ql0(0) = 1. / 2.  * sqrt( 1. / PI) * matS0 % ( 1 * matZk(0) % matR2k(0));
      ql0(1) = 1. / 2.  * sqrt( 3. / PI) * matS0 % ( 1 * matZk(1) % matR2k(0));
      ql0(2) = 1. / 4.  * sqrt( 5. / PI) * matS0 % ( 2 * matZk(2) % matR2k(0) -   1 * matR2k(1) % matZk(0));
      ql0(3) = 1. / 4.  * sqrt( 7. / PI) * matS0 % ( 2 * matZk(3) % matR2k(0) -   3 * matR2k(1) % matZk(1));
      ql0(4) = 1. / 16. * sqrt( 9. / PI) * matS0 % ( 8 * matZk(4) % matR2k(0) -  24 * matR2k(1) % matZk(2) +  3 * matR2k(2) % matZk(0));
      ql0(5) = 1. / 16. * sqrt(11. / PI) * matS0 % ( 8 * matZk(5) % matR2k(0) -  40 * matR2k(1) % matZk(3) + 15 * matR2k(2) % matZk(1));
      ql0(6) = 1. / 32. * sqrt(13. / PI) * matS0 % (16 * matZk(6) % matR2k(0) - 120 * matR2k(1) % matZk(4) + 90 * matR2k(2) % matZk(2) - 5 * matR2k(3) % matZk(0));
      break;
    case General::COMPAT_BERGER:
      /*
      *    QLM = sqrt( (4 * PI) / (2 * L + 1) ) * rho * r^l * YLM (see, e.g., Messiah, Albert (1999) "Quantum mechanics : two volumes bound as one")
      *    !!! Berger added a factor 2 in the Q20 definition !!!
      */
      ql0(0) = matS0 % (1.0 * matZk(0) % matR2k(0));
      ql0(1) = matS0 % (1.0 * matZk(1) % matR2k(0));
      ql0(2) = matS0 % (2.0 * matZk(2) % matR2k(0) - 1.0          * matR2k(1) % matZk(0));
      ql0(3) = matS0 % (1.0 * matZk(3) % matR2k(0) - 3.0   / 2.0  * matR2k(1) % matZk(1));
      ql0(4) = matS0 % (1.0 * matZk(4) % matR2k(0) - 24.0  / 8.0  * matR2k(1) % matZk(2) + 3.0  / 8.0  * matR2k(2) % matZk(0));
      ql0(5) = matS0 % (1.0 * matZk(5) % matR2k(0) - 40.0  / 8.0  * matR2k(1) % matZk(3) + 15.0 / 8.0  * matR2k(2) % matZk(1));
      ql0(6) = matS0 % (1.0 * matZk(6) % matR2k(0) - 120.0 / 16.0 * matR2k(1) % matZk(4) + 90.0 / 16.0 * matR2k(2) % matZk(2) - 5.0 / 16.0 * matR2k(3) % matZk(0));
      break;
    case General::COMPAT_ROBLEDO:
      // QLM = sqrt( (4 * PI) / (2 * L + 1) ) * rho * r^l * YLM (see, e.g., Messiah, Albert (1999) "Quantum mechanics : two volumes bound as one" or wikipedia)
      ql0(0) = matS0 % (1.0 * matZk(0) % matR2k(0));
      ql0(1) = matS0 % (1.0 * matZk(1) % matR2k(0));
      ql0(2) = matS0 % (1.0 * matZk(2) % matR2k(0) - 0.5          * matR2k(1) % matZk(0));
      ql0(3) = matS0 % (1.0 * matZk(3) % matR2k(0) - 3.0   / 2.0  * matR2k(1) % matZk(1));
      ql0(4) = matS0 % (1.0 * matZk(4) % matR2k(0) - 24.0  / 8.0  * matR2k(1) % matZk(2) + 3.0  / 8.0  * matR2k(2) % matZk(0));
      ql0(5) = matS0 % (1.0 * matZk(5) % matR2k(0) - 40.0  / 8.0  * matR2k(1) % matZk(3) + 15.0 / 8.0  * matR2k(2) % matZk(1));
      ql0(6) = matS0 % (1.0 * matZk(6) % matR2k(0) - 120.0 / 16.0 * matR2k(1) % matZk(4) + 90.0 / 16.0 * matR2k(2) % matZk(2) - 5.0 / 16.0 * matR2k(3) % matZk(0));
      break;
    case General::COMPAT_HFBTHO:
      /* HFBTHO uses the definitions of General::COMPAT_ROBLEDO for Q00,Q10,
      *  the definition of General::COMPAT_BERGER for Q20,
      *  and the definitions General::COMPAT_NONE for Q30 and higher.
      *  No, it doesn't make any sense.
      */
      ql0(0) =                             matS0 % ( 1 * matZk(0) % matR2k(0));
      ql0(1) =                             matS0 % ( 1 * matZk(1) % matR2k(0));
      ql0(2) =                             matS0 % ( 2 * matZk(2) % matR2k(0) -   1 * matR2k(1) % matZk(0));
      ql0(3) = 1. / 4.  * sqrt( 7. / PI) * matS0 % ( 2 * matZk(3) % matR2k(0) -   3 * matR2k(1) % matZk(1));
      ql0(4) = 1. / 16. * sqrt( 9. / PI) * matS0 % ( 8 * matZk(4) % matR2k(0) -  24 * matR2k(1) % matZk(2) +  3 * matR2k(2) % matZk(0));
      ql0(5) = 1. / 16. * sqrt(11. / PI) * matS0 % ( 8 * matZk(5) % matR2k(0) -  40 * matR2k(1) % matZk(3) + 15 * matR2k(2) % matZk(1));
      ql0(6) = 1. / 32. * sqrt(13. / PI) * matS0 % (16 * matZk(6) % matR2k(0) - 120 * matR2k(1) % matZk(4) + 90 * matR2k(2) % matZk(2) - 5 * matR2k(3) % matZk(0));
      break;
    default:
      ERROR("Unknown general.compatibility value.");
  }

  DBG_LEAVE;
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

  beta(NEUTRON) = getBetaFromQ20(N, A, qlm(2, NEUTRON));
  beta(PROTON ) = getBetaFromQ20(Z, A, qlm(2, PROTON ));
  beta(TOTAL  ) = getBetaFromQ20(A, A, qlm(2, TOTAL  ));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\langle m,n_a|\hat{r}^{2l}|m,n_b\rangle\f$ values.
 */

void MultipoleOperators::calcR2l(void)
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

  INT    kmax = 7;
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

void MultipoleOperators::calcQlm(const Multi<arma::mat> &rho)
{
  DBG_ENTER;

  Basis &basis = state.basis;

  if (rho(NEUTRON).empty() || rho(PROTON).empty())
  {
    DBG_LEAVE;
  }

  // dependencies
  calcQl0Matrices();
  basis.calcWDN();
  qlm             = arma::zeros(7, 3);
  arma::mat rn    = 2 * rho(NEUTRON);     // coeff 2 to take into account time reversal negative states
  arma::mat rp    = 2 * rho(PROTON );     // same
  arma::mat ra    = rn + rp;
  // arma::mat s     = basis.ORtoHO * basis.ORtoHO.t();

  qlm(0, NEUTRON) = arma::accu(ql0(0) % rn);     // tr( M * r * M.t() ) ?
  qlm(0, PROTON)  = arma::accu(ql0(0) % rp);
  qlm(0, TOTAL)   = arma::accu(ql0(0) % ra);
  qlm(1, NEUTRON) = arma::accu(ql0(1) % rn);     // tr( Mt M A Mt M r  ) ?
  qlm(1, PROTON)  = arma::accu(ql0(1) % rp);
  qlm(1, TOTAL)   = arma::accu(ql0(1) % ra);
  qlm(2, NEUTRON) = arma::accu(ql0(2) % rn);
  qlm(2, PROTON)  = arma::accu(ql0(2) % rp);
  qlm(2, TOTAL)   = arma::accu(ql0(2) % ra);
  qlm(3, NEUTRON) = arma::accu(ql0(3) % rn);
  qlm(3, PROTON)  = arma::accu(ql0(3) % rp);
  qlm(3, TOTAL)   = arma::accu(ql0(3) % ra);
  qlm(4, NEUTRON) = arma::accu(ql0(4) % rn);
  qlm(4, PROTON)  = arma::accu(ql0(4) % rp);
  qlm(4, TOTAL)   = arma::accu(ql0(4) % ra);
  qlm(5, NEUTRON) = arma::accu(ql0(5) % rn);
  qlm(5, PROTON)  = arma::accu(ql0(5) % rp);
  qlm(5, TOTAL)   = arma::accu(ql0(5) % ra);
  qlm(6, NEUTRON) = arma::accu(ql0(6) % rn);
  qlm(6, PROTON)  = arma::accu(ql0(6) % rp);
  qlm(6, TOTAL)   = arma::accu(ql0(6) % ra);

  calcBeta();
  calcNPart(rho);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice table with the multipole moment mean values.
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
    nPart(iso) = qlm(0, iso) / correctionFactor;
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

  /* ASSERT(state.checkSolution(), "state does not contain a solution."); */

  // color constrained values
  Multi<std::string> colorStr;

  for (auto &c : state.constraints) colorStr(c.second.lm, c.second.iso) = TABLE_BLUE;

  std::list<std::list<std::string> > tempList;

  tempList.push_back(
  {
    "Mass    [n]",
    PF("%.6f", nPart(NEUTRON)),
    PF("%.6f", nPart(PROTON )),
    PF("%.6f", nPart(TOTAL  ))
  });

  for (INT lambda = 0; lambda < 7; lambda++)
  {
    tempList.push_back(
    {
      PF("<Q%1d0> [fm%1d]", lambda, lambda),
      colorStr(lambda, NEUTRON) + PF("%.6f", qlm(lambda, NEUTRON)) + TABLE_NORM,
      colorStr(lambda, PROTON ) + PF("%.6f", qlm(lambda, PROTON )) + TABLE_NORM,
      colorStr(lambda, TOTAL  ) + PF("%.6f", qlm(lambda, TOTAL  )) + TABLE_NORM
    });
  }

  tempList.push_back({"beta2    []", PF("%.6f", beta(NEUTRON)), PF("%.6f", beta(PROTON)),
                      PF("%.6f", beta(TOTAL))});

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
    if (qlm.empty())
    {
      DBG_RETURN(PF_RED("empty"));
    }

    for (UINT l = 0; l < 4; l++)
    {
      list.push_back({PF("Q%1d0", l), PF("%.6e", qlm(l, TOTAL))});
    }

    list.push_back({"beta", PF("%.6e", beta(TOTAL))});
  }
  else
  {
    list.push_back({"MultipoleOperators", ""});

    if (qlm.empty())
    {
      list.push_back({"qlm", PF_RED("empty")});
    }
    else
    {
      for (UINT l = 0; l < qlm.n_rows; l++)
      {
        list.push_back({PF("<Q%1d0> ", l),
                        "(" + PF_YELLOW("neut.") + ":" + PF_BLUE("%13.6e", qlm(l, NEUTRON)) + ", " +
                              PF_YELLOW("prot.") + ":" + PF_BLUE("%13.6e", qlm(l, PROTON) ) + ", " +
                              PF_YELLOW("total") + ":" + PF_BLUE("%13.6e", qlm(l, TOTAL)  ) + ")"});
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
