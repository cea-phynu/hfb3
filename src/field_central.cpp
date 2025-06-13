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

#include "field_central.h"
#include "tools.h"
#include "basis.h"
#include "interaction.h"
#include "multi.h"

/** \file
 *  \brief Methods of the FieldCentral class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldCentral::FieldCentral(Parameters _parameters, State *_state) :
  Field(_parameters, _state)
{
  DBG_ENTER;

  name = "central";
  shortName = "Centr.";

  IzDirect = FMulti<arma::cube>({SRANGE, SRANGE, NZRANGE, NZRANGE});

  IzExchange = FMulti<arma::cube>({SRANGE, SRANGE, NZRANGE, NZRANGE});

  w = parameters["w"];
  b = parameters["b"];
  h = parameters["h"];
  m = parameters["m"];
  p = parameters["p"];

  // INFO("central field parameters: w=%f b=%f h=%f m=%f p=%f", w, b, h, m, p);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldCentral::calcField(void)
{
  DBG_ENTER;

  auto startTime = Tools::clock();

  if (!field(NEUTRON, DIRECT  ).empty() && !field(PROTON , DIRECT  ).empty() &&
      !field(NEUTRON, EXCHANGE).empty() && !field(PROTON , EXCHANGE).empty() &&
      !field(NEUTRON, PAIRING ).empty() && !field(PROTON , PAIRING ).empty()) DBG_LEAVE;


  // Tools::info("rhon  ", state.rho(NEUTRON));
  // Tools::info("rhop  ", state.rho(PROTON ));
  // Tools::info("kappan", state.kappa(NEUTRON));
  // Tools::info("kappap", state.kappa(PROTON ));


  // dependencies
  calcIr();
  calcIz();

  //Usefull quantities
  INT M = basis.mMax;
  INT dMax = basis.dMax;
  // INT N_p = basis.nGlobalMax;
  INT N_z = basis.n_zGlobalMax;
  INT DMax = (dMax - 1) * 3 + 1;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Kappa and Rho Multis with explicit quantum numbers.
  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});

  Multi<arma::mat> Rho;   //Rho  (Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)
  Multi<arma::mat> Kappa; //Kappa(Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)

  //NOTE: Fmulti is slower on sparse matrix ?

  // FMulti<arma::mat> Rho(
  //   {0,0,0,0  ,0  ,0,0},
  //   {2,M,M,N_p,N_p,4,4}
  // );
  // FMulti<arma::mat> Kappa(
  //   {0,0,0,0  ,0  ,0,0},
  //   {2,M,M,N_p,N_p,4,4}
  // );

  for (UINT i = 0; i < nbBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m = bras(0, 0);
    INT n = bras(1, 0);
    INT d = bras(3, 0);
    INT s = bras(4, 0);

    for (UINT j = 0; j < nbBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT mp = kets(0, 0);
      INT np = kets(1, 0);
      INT dp = kets(3, 0);
      INT sp = kets(4, 0);

      if (mp - sp != m - s) continue; //Symétrie block omega

      {
        Rho(NEUTRON, m, mp, n, np, d + 2 * dp, s + 2 * sp) = rho(NEUTRON  ).submat(bras.filter, kets.filter);
        Rho(PROTON , m, mp, n, np, d + 2 * dp, s + 2 * sp) = rho(PROTON    ).submat(bras.filter, kets.filter);
        Kappa(NEUTRON, m, mp, n, np, d + 2 * dp, s + 2 * sp) = kappa(NEUTRON).submat(bras.filter, kets.filter);
        Kappa(PROTON , m, mp, n, np, d + 2 * dp, s + 2 * sp) = kappa(PROTON  ).submat(bras.filter, kets.filter);
      }//if
    }//j
  }//i



  ////////////////////////////////////////////////////////////////////////////////////////////////////////One builds Rmat and Kmat(pairing).
  //Usefull quantities for Rmat and Kmat.

  arma::wall_clock timer;
  timer.tic();

  double fac_direct_iso;
  double fac_direct_noniso;
  double fac_E_iso_pp;
  double fac_E_noniso_pp;
  double fac_E_iso_mm;
  double fac_E_noniso_mm;
  double fac_P_pp;
  double fac_P_mm;

  fac_direct_iso    = 2 * w + b - (2 * h + m);
  fac_direct_noniso = 2 * w + b;

  fac_E_iso_pp    = -w - b + h + m;
  fac_E_iso_mm    = -b + m;
  fac_E_noniso_pp =  h + m;
  fac_E_noniso_mm =  m;

  fac_P_pp    =  w  - h;
  fac_P_mm    =  -b + m;

  // Multi<arma::cube> RmatDirect;     //Rmat for the mean central direct term. RmatDirect(ig,Isospin,mp,n_beta,n_delta)(n_zbeta,n_zdelta,d_beta+2*d_delta)
  // Multi<arma::cube> RmatExchange_0; //Rmat for the mean central exchange term when m_alpha = m_gamma.    RmatExchange_0(ig,Isospin,m_beta,n_beta,n_delta)(n_zbeta,n_zdelta,d_beta+2*d_delta).
  // Multi<arma::cube> RmatExchange_1; //Rmat for the mean central exchange term when m_alpha = m_gamma +1. RmatExchange_0(ig,Isospin,m_beta,n_beta,n_delta)(n_zbeta,n_zdelta,d_beta+2*d_delta).
  // Multi<arma::cube> Kmat_0;         //Kmat for the central pairing term when m_alpha = m_gamma.    Kmat_0(ig,Isospin,m_beta,n_beta,n_delta)(n_zbeta,n_zdelta,d_beta+2*d_delta).
  // Multi<arma::cube> Kmat_1;         //Kmat for the central pairing term when m_alpha = m_gamma +1. Kmat_1(ig,Isospin,m_beta,n_beta,n_delta)(n_zbeta,n_zdelta,d_beta+2*d_delta).

  FMulti<arma::cube> RmatDirect(    {SRANGE, MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> RmatExchange_0({SRANGE, MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> RmatExchange_1({SRANGE, MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Kmat_0(        {SRANGE, MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Kmat_1(        {SRANGE, MRANGE, NRANGE, NRANGE});

  //calculates RmatDirect, RmatExchange_0 and Kmat_0.
  INT mp = 0; //mp=0 requires a special treatment to avoid negative omega blocks.
  {
    for (int n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
    {
      for (int n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        RmatDirect(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatDirect(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        RmatExchange_0(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatExchange_0(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        Kmat_0(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        Kmat_0(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        for (int d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (int d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            RmatDirect(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_direct_iso * (Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0))
                + fac_direct_noniso * (Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0));
            RmatDirect(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_direct_noniso * (Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0))
                + fac_direct_iso * (Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0));

            RmatExchange_0(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_noniso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_iso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_noniso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);

            RmatExchange_0(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_noniso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_iso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_noniso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);

            Kmat_0(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_mm * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                                                                                  + fac_P_pp * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);

            Kmat_0(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_pp * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_P_mm * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);

          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp

  for (int mp = 1 ; mp < M ; mp++)
  {
    for (int n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
    {
      for (int n_delta = 0 ; n_delta < basis.nMax(mp); n_delta++)
      {
        RmatDirect(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatDirect(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        RmatExchange_0(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatExchange_0(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatExchange_0(NEUTRON, -mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        RmatExchange_0(PROTON , -mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        Kmat_0(NEUTRON, mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        Kmat_0(PROTON , mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        Kmat_0(NEUTRON, -mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);
        Kmat_0(PROTON , -mp, n_beta, n_delta) = arma::zeros(basis.n_zMax(mp, n_beta), basis.n_zMax(mp, n_delta), DMax);

        for (int d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (int d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            RmatDirect(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_direct_iso * (Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3))
                + fac_direct_noniso * (Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3));
            RmatDirect(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_direct_noniso * (Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3))
                + fac_direct_iso * (Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3));

            RmatExchange_0(NEUTRON, -mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3) + fac_E_iso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_noniso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3) + fac_E_noniso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);
            RmatExchange_0(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_iso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3)
                + fac_E_noniso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_noniso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3);

            RmatExchange_0(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_iso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3)
                + fac_E_noniso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0) + fac_E_noniso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3);
            RmatExchange_0(PROTON , -mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_E_iso_pp * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3) + fac_E_iso_mm * Rho(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_E_noniso_pp * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3) + fac_E_noniso_mm * Rho(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);

            Kmat_0(NEUTRON, -mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_pp * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3)
                + fac_P_mm * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);
            Kmat_0(NEUTRON, mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_pp * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_P_mm * Kappa(NEUTRON, mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3);

            Kmat_0(PROTON , mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_pp * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0)
                + fac_P_mm * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3);
            Kmat_0(PROTON , -mp, n_beta, n_delta).slice(d_beta + 2 * d_delta) = fac_P_pp * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 3)
                + fac_P_mm * Kappa(PROTON , mp, mp, n_beta, n_delta, d_beta + 2 * d_delta, 0);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp

  //One has to separate mp <= 0 and mp >0 because the structures of Rmat and Kmat differ at that point.
  for (int mp = -M + 2 ; mp < 1 ; mp++) //mp <= 0.
  {
    for (int n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
    {
      for (int n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
      {
        RmatExchange_1(NEUTRON, mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);
        RmatExchange_1(PROTON , mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);

        Kmat_1(NEUTRON, mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);
        Kmat_1(PROTON , mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);

        for (int d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (int d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            INT D = 2 * d_beta + d_delta;
            RmatExchange_1(NEUTRON, mp, n_delta, n_beta).slice(D) = (w - h) * Rho(NEUTRON, abs(mp), abs(mp - 1), n_delta, n_beta, D, 2)
                - h * Rho(PROTON , abs(mp), abs(mp - 1), n_delta, n_beta, D, 2);
            RmatExchange_1(PROTON , mp, n_delta, n_beta).slice(D) = (w - h) * Rho(PROTON , abs(mp), abs(mp - 1), n_delta, n_beta, D, 2)
                - h * Rho(NEUTRON, abs(mp), abs(mp - 1), n_delta, n_beta, D, 2);

            Kmat_1(NEUTRON, mp, n_delta, n_beta).slice(D) = -(w - h + b - m) * Kappa(NEUTRON, abs(mp), abs(mp - 1), n_delta, n_beta, D, 2);
            Kmat_1(PROTON , mp, n_delta, n_beta).slice(D) = -(w - h + b - m) * Kappa(PROTON , abs(mp), abs(mp - 1), n_delta, n_beta, D, 2);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp <= 0

  for (int mp = 1 ; mp < M ; mp++) //mp > 0.
  {
    for (int n_beta = 0 ; n_beta < basis.nMax(mp - 1) ; n_beta++)
    {
      for (int n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        RmatExchange_1(NEUTRON, mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);
        RmatExchange_1(PROTON , mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);

        Kmat_1(NEUTRON, mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);
        Kmat_1(PROTON , mp, n_delta, n_beta) = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);

        for (int d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (int d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            INT D = 2 * d_beta + d_delta;
            RmatExchange_1(NEUTRON, mp, n_delta, n_beta).slice(D) = (-w + h) * Rho(NEUTRON, mp, mp - 1, n_delta, n_beta, D, 1)
                + h * Rho(PROTON , mp, mp - 1, n_delta, n_beta, D, 1);
            RmatExchange_1(PROTON , mp, n_delta, n_beta).slice(D) = (-w + h) * Rho(PROTON , mp, mp - 1, n_delta, n_beta, D, 1)
                + h * Rho(NEUTRON, mp, mp - 1, n_delta, n_beta, D, 1);

            Kmat_1(NEUTRON, mp, n_delta, n_beta).slice(D) = (w - h + b - m) * Kappa(NEUTRON, mp, mp - 1, n_delta, n_beta, D, 1);
            Kmat_1(PROTON , mp, n_delta, n_beta).slice(D) = (w - h + b - m) * Kappa(PROTON , mp, mp - 1, n_delta, n_beta, D, 1);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp > 0

  ramTable["Rmat"] =
    RmatDirect.size()
  + RmatExchange_0.size()
  + RmatExchange_1.size()
  + Kmat_0.size()
  + Kmat_1.size();

  timeTable["Rmat"] = timer.toc();
  timer.tic();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Calculates IzR and IzK merging Iz with Rmat and Kmat.
  // Multi<arma::vec>  IzRDirect; //RmatDirect with IzDirect. IzRDirect(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_beta + n_delta*basis.nMax(mp) + indexM(mp));
  // Multi<arma::vec>  IzRExchange_0; //RmatExchange_1 with IzExchange. IzRExchange_1(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index3M(mp+M-2));
  // Multi<arma::vec>  IzRExchange_1; //RmatExchange_1 with IzExchange. IzRExchange_1(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index3M(mp+M-2));
  // Multi<arma::vec>  IzRExchange_3; //RmatExchange_0 with IzExchange in a reverse order. IzRExchange_3(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index2M(-mp+M-1));
  //
  // Multi<arma::vec>  IzK_0; //Kmat_0 with IzExchange. IzK_0(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_beta + n_delta*basis.nMax(mp) + indexM(mp));
  // Multi<arma::vec>  IzK_1; //Kmat_1 with IzExchange. Izk_1(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index3M(mp+M-2));
  // Multi<arma::vec>  IzK_3; //Kmat_0 with IzExchange in a reverse order. IzK_3(ISO,ig,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index2M(-mp+M-1));

  FMulti<arma::vec>  IzRDirect(    {SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzRExchange_0({SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzRExchange_1({SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzRExchange_3({SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_0(        {SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_1(        {SRANGE, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_3(        {SRANGE, DDRANGE, NZRANGE, NZRANGE});


  INT sizeInit = 0;
  INT sizeInit2 = 0;
  INT sizeInit3 = 0;
  IVEC indexM(M);
  IVEC index2M(2 * M);
  IVEC index3M(2 * M - 1);

  for (INT m = 0 ; m < M ; m++)
  {
    indexM(m) = sizeInit;
    sizeInit += basis.nMax(m) * basis.nMax(m);
  }

  for (INT m = -M + 1 ; m < M ; m++)
  {
    index2M(m + M - 1) = sizeInit2;
    sizeInit2 += basis.nMax(abs(m)) * basis.nMax(abs(m));
  }

  for (INT m = -M + 2 ; m < M ; m++)
  {
    index3M(m + M - 2) = sizeInit3;
    sizeInit3 += basis.nMax(abs(m)) * basis.nMax(abs(m - 1));
  }

  arma::cube IzcubeD;

  //Calculates IzRDirect.
  for (int n_zalpha = 0; n_zalpha < N_z ; n_zalpha++)
  {
    for (int n_zgamma = 0; n_zgamma < n_zalpha + 1 ; n_zgamma++) //One uses the symmetry gamma <-> alpha.
    {
      for (int d_alpha = 0; d_alpha < dMax ; d_alpha++)
      {
        for (int d_gamma = 0; d_gamma < dMax ; d_gamma++)
        {
          INT D = d_alpha + 2 * d_gamma;
          INT D_bar = 2 * d_alpha + d_gamma;
          arma::vec &IzRDN = IzRDirect(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzRDP = IzRDirect(PROTON , D, n_zalpha, n_zgamma);
          IzRDN = arma::zeros(sizeInit);
          IzRDP = arma::zeros(sizeInit);

          for (int mp = 0 ; mp < M ; mp++)
          {
            for (int n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
            {
              for (int n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++) //One uses the symmetry beta <-> delta.
              {
                INT N = n_beta + n_delta * basis.nMax(mp) + indexM(mp);
                INT N_sym = n_delta + n_beta * basis.nMax(mp) + indexM(mp);

                IzcubeD = IzDirect(d_alpha, d_gamma, n_zalpha, n_zgamma).subcube(
                  0, 0, 0,
                  basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_beta) - 1, DMax - 1
                );
                IzRDN(N) = ACCU(IzcubeD % RmatDirect(NEUTRON, mp, n_delta, n_beta));
                IzRDN(N_sym) = IzRDN(N);
                IzRDP(N) = ACCU(IzcubeD % RmatDirect(PROTON , mp, n_delta, n_beta));
                IzRDP(N_sym) = IzRDP(N);
              }//np
            }//n
          }//mp

          IzRDirect(NEUTRON, D_bar, n_zgamma, n_zalpha) = IzRDN;
          IzRDirect(PROTON , D_bar, n_zgamma, n_zalpha) = IzRDP;
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  arma::cube IzcubeE_0;
  //Calculates IzRExchange_0, IzRExchange_3, IzK_0 and IzK_3.
  for (int n_zalpha = 0; n_zalpha < N_z ; n_zalpha++)
  {
    for (int n_zgamma = 0; n_zgamma < n_zalpha + 1 ; n_zgamma++) //One uses the symmetry gamma <-> alpha.
    {
      for (int d_alpha = 0; d_alpha < dMax ; d_alpha++)
      {
        for (int d_gamma = 0; d_gamma < dMax ; d_gamma++)
        {
          INT D = d_alpha + 2 * d_gamma;
          INT D_sym = d_gamma + 2 * d_alpha;

          arma::vec &IzREN_0 = IzRExchange_0(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzREP_0 = IzRExchange_0(PROTON , D, n_zalpha, n_zgamma);
          arma::vec &IzREN_0_sym = IzRExchange_0(NEUTRON, D_sym, n_zgamma, n_zalpha);
          arma::vec &IzREP_0_sym = IzRExchange_0(PROTON , D_sym, n_zgamma, n_zalpha);
          IzREN_0 = arma::zeros(sizeInit2);
          IzREP_0 = arma::zeros(sizeInit2);
          IzREN_0_sym = arma::zeros(sizeInit2);
          IzREP_0_sym = arma::zeros(sizeInit2);

          arma::vec &IzREN_3 = IzRExchange_3(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzREP_3 = IzRExchange_3(PROTON , D, n_zalpha, n_zgamma);
          arma::vec &IzREN_3_sym = IzRExchange_3(NEUTRON, D_sym, n_zgamma, n_zalpha);
          arma::vec &IzREP_3_sym = IzRExchange_3(PROTON , D_sym, n_zgamma, n_zalpha);
          IzREN_3 = arma::zeros(sizeInit2);
          IzREP_3 = arma::zeros(sizeInit2);
          IzREN_3_sym = arma::zeros(sizeInit2);
          IzREP_3_sym = arma::zeros(sizeInit2);


          arma::vec &IzKN_0 = IzK_0(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzKP_0 = IzK_0(PROTON , D, n_zalpha, n_zgamma);
          arma::vec &IzKN_0_sym = IzK_0(NEUTRON, D_sym, n_zgamma, n_zalpha);
          arma::vec &IzKP_0_sym = IzK_0(PROTON , D_sym, n_zgamma, n_zalpha);
          IzKN_0 = arma::zeros(sizeInit2);
          IzKP_0 = arma::zeros(sizeInit2);
          IzKN_0_sym = arma::zeros(sizeInit2);
          IzKP_0_sym = arma::zeros(sizeInit2);

          arma::vec &IzKN_3 = IzK_3(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzKP_3 = IzK_3(PROTON , D, n_zalpha, n_zgamma);
          arma::vec &IzKN_3_sym = IzK_3(NEUTRON, D_sym, n_zgamma, n_zalpha);
          arma::vec &IzKP_3_sym = IzK_3(PROTON , D_sym, n_zgamma, n_zalpha);
          IzKN_3 = arma::zeros(sizeInit2);
          IzKP_3 = arma::zeros(sizeInit2);
          IzKN_3_sym = arma::zeros(sizeInit2);
          IzKP_3_sym = arma::zeros(sizeInit2);

          for (int mp = -M + 1 ; mp < M ; mp++)
          {
            for (int n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
            {
              for (int n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
              {
                IzcubeE_0 = IzExchange(d_alpha, d_gamma, n_zalpha, n_zgamma).subcube(0, 0, 0, basis.n_zMax(abs(mp), n_delta) - 1, basis.n_zMax(abs(mp), n_beta) - 1, DMax - 1);
                INT N = n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1);
                INT N_sym = n_beta + n_delta * basis.nMax(abs(mp)) + index2M(mp + M - 1);
                IzREN_0(N) = ACCU(IzcubeE_0 % RmatExchange_0(NEUTRON, mp, n_delta, n_beta));
                IzREP_0(N) = ACCU(IzcubeE_0 % RmatExchange_0(PROTON , mp, n_delta, n_beta));
                IzREN_0_sym(N_sym) = IzREN_0(N);
                IzREP_0_sym(N_sym) = IzREP_0(N);

                INT Np = n_delta + n_beta * basis.nMax(abs(mp)) + index2M(-mp + M - 1);
                INT Np_sym = n_beta + n_delta * basis.nMax(abs(mp)) + index2M(-mp + M - 1);
                IzREN_3(Np) = IzREN_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                IzREP_3(Np) = IzREP_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                IzREN_3_sym(Np_sym) = IzREN_3(Np);
                IzREP_3_sym(Np_sym) = IzREP_3(Np);

                IzKN_0(N) = ACCU(IzcubeE_0 % Kmat_0(NEUTRON, mp, n_delta, n_beta));
                IzKP_0(N) = ACCU(IzcubeE_0 % Kmat_0(PROTON , mp, n_delta, n_beta));
                IzKN_0_sym(N_sym) = IzKN_0(N);
                IzKP_0_sym(N_sym) = IzKP_0(N);

                IzKN_3(Np) = IzKN_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                IzKP_3(Np) = IzKP_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                IzKN_3_sym(Np_sym) = IzKN_3(Np);
                IzKP_3_sym(Np_sym) = IzKP_3(Np);

              }//n_beta
            }//n_delta
          }//mp
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  arma::cube IzcubeE_1;
  //Calculates IzRExchange_1 and IzK_1.
  for (int n_zalpha = 0; n_zalpha < N_z ; n_zalpha++)
  {
    for (int n_zgamma = 0; n_zgamma < N_z ; n_zgamma++)
    {
      for (int d_alpha = 0; d_alpha < dMax ; d_alpha++)
      {
        for (int d_gamma = 0; d_gamma < dMax ; d_gamma++)
        {
          INT D = d_alpha + 2 * d_gamma;

          arma::vec &IzKN_1 = IzK_1(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzKP_1 = IzK_1(PROTON , D, n_zalpha, n_zgamma);
          IzKN_1 = arma::zeros(sizeInit3);
          IzKP_1 = arma::zeros(sizeInit3);

          arma::vec &IzREN_1 = IzRExchange_1(NEUTRON, D, n_zalpha, n_zgamma);
          arma::vec &IzREP_1 = IzRExchange_1(PROTON , D, n_zalpha, n_zgamma);
          IzREN_1 = arma::zeros(sizeInit3);
          IzREP_1 = arma::zeros(sizeInit3);

          for (int mp = -M + 2 ; mp < M ; mp++)
          {
            for (int n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
            {
              for (int n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
              {
                INT N = n_delta + n_beta * basis.nMax(abs(mp)) + index3M(mp + M - 2);
                IzcubeE_1 = IzExchange(d_alpha, d_gamma, n_zalpha, n_zgamma).subcube(0, 0, 0, basis.n_zMax(abs(mp), n_delta) - 1, basis.n_zMax(abs(mp - 1), n_beta) - 1, DMax - 1);

                IzKN_1(N) = ACCU(IzcubeE_1 % Kmat_1(NEUTRON, mp, n_delta, n_beta));
                IzKP_1(N) = ACCU(IzcubeE_1 % Kmat_1(PROTON , mp, n_delta, n_beta));

                IzREN_1(N) = ACCU(IzcubeE_1 % RmatExchange_1(NEUTRON, mp, n_delta, n_beta));
                IzREP_1(N) = ACCU(IzcubeE_1 % RmatExchange_1(PROTON , mp, n_delta, n_beta));
              }//np
            }//n
          }//mp
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  ramTable["IzR"] =
    IzRDirect.size()
  + IzRExchange_0.size()
  + IzRExchange_1.size()
  + IzRExchange_3.size()
  + IzK_0.size()
  + IzK_1.size()
  + IzK_3.size();
  timeTable["IzR"] = timer.toc();
  timer.tic();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Gamma_Direct, Gamma_Exchange and Delta from merging Ir with IzR and IzK.
  // Multi<arma::cube> Gamma_Direct; //GammaDirect(Iso,m,n_alpha,n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Gamma_Exchange; //GammaDirect(Iso,m,n_alpha,n_gamma,s_alpha+2*s_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Delta; ////Delta(Iso,m,n_alpha,n_gamma,s_alpha+2*s_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).

  FMulti<arma::cube> Gamma_Direct(  {SRANGE, MPRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Gamma_Exchange({SRANGE, MPRANGE, NRANGE, NRANGE, {0, 4}});
  FMulti<arma::cube> Delta(         {SRANGE, MPRANGE, NRANGE, NRANGE, {0, 4}});

  for (INT m = 0 ; m < M ; m++)
  {
    for (INT n_alpha = 0; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0; n_gamma < n_alpha + 1 ; n_gamma++) //One uses the symmetry gamma <-> alpha.
      {
        arma::cube &GammaCubeDN = Gamma_Direct(NEUTRON, m, n_alpha, n_gamma);
        arma::cube &GammaCubeDP = Gamma_Direct(PROTON , m, n_alpha, n_gamma);
        arma::cube &GammaCubeDN_sym = Gamma_Direct(NEUTRON, m, n_gamma, n_alpha);
        arma::cube &GammaCubeDP_sym = Gamma_Direct(PROTON , m, n_gamma, n_alpha);
        GammaCubeDN = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeDP = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeDN_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        GammaCubeDP_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

        arma::cube &GammaCubeE_0N = Gamma_Exchange(NEUTRON, m, n_alpha, n_gamma, 0);
        arma::cube &GammaCubeE_0P = Gamma_Exchange(PROTON , m, n_alpha, n_gamma, 0);
        arma::cube &GammaCubeE_3N = Gamma_Exchange(NEUTRON, m, n_alpha, n_gamma, 3);
        arma::cube &GammaCubeE_3P = Gamma_Exchange(PROTON , m, n_alpha, n_gamma, 3);
        arma::cube &GammaCubeE_0N_sym = Gamma_Exchange(NEUTRON, m, n_gamma, n_alpha, 0);
        arma::cube &GammaCubeE_0P_sym = Gamma_Exchange(PROTON , m, n_gamma, n_alpha, 0);
        arma::cube &GammaCubeE_3N_sym = Gamma_Exchange(NEUTRON, m, n_gamma, n_alpha, 3);
        arma::cube &GammaCubeE_3P_sym = Gamma_Exchange(PROTON , m, n_gamma, n_alpha, 3);
        GammaCubeE_0N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeE_0P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeE_3N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeE_3P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        GammaCubeE_0N_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        GammaCubeE_0P_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        GammaCubeE_3N_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        GammaCubeE_3P_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);



        arma::cube &DeltaCube_0N     = Delta(NEUTRON, m, n_alpha, n_gamma, 0);
        arma::cube &DeltaCube_0P     = Delta(PROTON , m, n_alpha, n_gamma, 0);
        arma::cube &DeltaCube_3N     = Delta(NEUTRON, m, n_alpha, n_gamma, 3);
        arma::cube &DeltaCube_3P     = Delta(PROTON , m, n_alpha, n_gamma, 3);
        arma::cube &DeltaCube_0N_sym = Delta(NEUTRON, m, n_gamma, n_alpha, 0);
        arma::cube &DeltaCube_0P_sym = Delta(PROTON , m, n_gamma, n_alpha, 0);
        arma::cube &DeltaCube_3N_sym = Delta(NEUTRON, m, n_gamma, n_alpha, 3);
        arma::cube &DeltaCube_3P_sym = Delta(PROTON , m, n_gamma, n_alpha, 3);

        DeltaCube_0N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        DeltaCube_0P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        DeltaCube_3N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        DeltaCube_3P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        DeltaCube_0N_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        DeltaCube_0P_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        DeltaCube_3N_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        DeltaCube_3P_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);


        arma::vec &IrvecD = IrDirect(m, n_alpha, n_gamma);
        arma::vec &IrvecE_0 = IrExchange(0, m, n_alpha, n_gamma);

        for (int d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (int d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            INT D = d_alpha + 2 * d_gamma;

            for (int n_zalpha = 0; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
            {
              for (int n_zgamma = 0; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
              {
                GammaCubeDN (n_zalpha, n_zgamma, D) += ACCU(IrvecD % IzRDirect(NEUTRON, D, n_zalpha, n_zgamma));
                GammaCubeDP (n_zalpha, n_zgamma, D) += ACCU(IrvecD % IzRDirect(PROTON , D, n_zalpha, n_zgamma));
                GammaCubeE_0N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzRExchange_0(NEUTRON, D, n_zalpha, n_zgamma));
                GammaCubeE_0P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzRExchange_0(PROTON , D, n_zalpha, n_zgamma));
                GammaCubeE_3N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzRExchange_3(NEUTRON, D, n_zalpha, n_zgamma));
                GammaCubeE_3P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzRExchange_3(PROTON , D, n_zalpha, n_zgamma));

                DeltaCube_0N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzK_0(NEUTRON, D, n_zalpha, n_zgamma));
                DeltaCube_0P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzK_0(PROTON , D, n_zalpha, n_zgamma));
                DeltaCube_3N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzK_3(NEUTRON, D, n_zalpha, n_zgamma));
                DeltaCube_3P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_0 % IzK_3(PROTON , D, n_zalpha, n_zgamma));
              }//n_zgamma
            }//n_zalpha
          }//d_gamma
        }//d_alpha

        for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            INT D = d_alpha + 2 * d_gamma;
            INT D_bar = d_gamma + 2 * d_alpha;
            GammaCubeDN_sym.slice(D_bar) = (GammaCubeDN.slice(D)).t();
            GammaCubeDP_sym.slice(D_bar) = (GammaCubeDP.slice(D)).t();

            GammaCubeE_0N_sym.slice(D_bar) = (GammaCubeE_0N.slice(D)).t();
            GammaCubeE_0P_sym.slice(D_bar) = (GammaCubeE_0P.slice(D)).t();
            GammaCubeE_3N_sym.slice(D_bar) = (GammaCubeE_3N.slice(D)).t();
            GammaCubeE_3P_sym.slice(D_bar) = (GammaCubeE_3P.slice(D)).t();

            DeltaCube_0N_sym.slice(D_bar)  = (DeltaCube_0N.slice(D)).t();
            DeltaCube_0P_sym.slice(D_bar)  = (DeltaCube_0P.slice(D)).t();
            DeltaCube_3N_sym.slice(D_bar)  = (DeltaCube_3N.slice(D)).t();
            DeltaCube_3P_sym.slice(D_bar)  = (DeltaCube_3P.slice(D)).t();
          }//d_gamma
        }//d_alpha

        // if ((m == 0) && (n_gamma == 0) && (n_alpha == 0))
        // {
        //   Tools::info("cube 0N ", DeltaCube_0N);
        // }

      }//n_gamma
    }//n_alpha
  }//m

  for (INT m = 1 ; m < M ; m++)
  {
    for (INT n_alpha = 0; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        arma::cube &GammaCubeE_1N = Gamma_Exchange(NEUTRON, m, n_alpha, n_gamma, 1);
        arma::cube &GammaCubeE_1P = Gamma_Exchange(PROTON , m, n_alpha, n_gamma, 1);
        GammaCubeE_1N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);
        GammaCubeE_1P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);

        arma::cube &DeltaCube_1N = Delta(NEUTRON, m, n_alpha, n_gamma, 1);
        arma::cube &DeltaCube_1P = Delta(PROTON , m, n_alpha, n_gamma, 1);
        DeltaCube_1N = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);
        DeltaCube_1P = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);

        arma::vec &IrvecE_1 = IrExchange(1, m, n_alpha, n_gamma);

        for (int d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (int d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            for (int n_zalpha = 0; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
            {
              for (int n_zgamma = 0; n_zgamma < basis.n_zMax(m - 1, n_gamma) ; n_zgamma++)
              {
                INT D = d_alpha + 2 * d_gamma;
                GammaCubeE_1N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_1 % IzRExchange_1(NEUTRON, D, n_zalpha, n_zgamma));
                GammaCubeE_1P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_1 % IzRExchange_1(PROTON , D, n_zalpha, n_zgamma));

                DeltaCube_1N (n_zalpha, n_zgamma, D) += ACCU(IrvecE_1 % IzK_1(NEUTRON, D, n_zalpha, n_zgamma));
                DeltaCube_1P (n_zalpha, n_zgamma, D) += ACCU(IrvecE_1 % IzK_1(PROTON , D, n_zalpha, n_zgamma));
              }//n_zgamma
            }//n_zalpha
          }//d_gamma
        }//d_alpha
      }//n_gamma
    }//n_alpha
  }//m

  ramTable["Gamma"] =
    Gamma_Direct.size()
  + Gamma_Exchange.size()
  + Delta.size();
  timeTable["Gamma"] = timer.toc();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Builds the interaction.

  field(NEUTRON, DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(NEUTRON, EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(NEUTRON, PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  UINT nbMBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});

  for (UINT i = 0; i < nbMBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT n_alpha =   bras(1, 0);
    INT d_alpha =   bras(3, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT n_gamma = kets(1, 0);
      INT d_gamma = kets(3, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma && s_alpha == s_gamma)//Mean Direct part.
      {
        arma::mat &GmatN = Gamma_Direct(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        arma::mat &GmatP = Gamma_Direct(PROTON , m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = GmatN;
        field(PROTON , DIRECT).submat(bras.filter, kets.filter) = GmatP;
      }//if

      if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0) //Mean Exchange part ++ and Pairing part ++.
      {
        arma::mat &GmatN = Gamma_Exchange(NEUTRON, m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);
        arma::mat &GmatP = Gamma_Exchange(PROTON , m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GmatN;
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GmatP;

        arma::mat &DmatN = Delta(NEUTRON, m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);
        arma::mat &DmatP = Delta(PROTON , m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) = DmatN;
        field(PROTON , PAIRING).submat(bras.filter, kets.filter) = DmatP;
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean exchange part -- and Pairing part --.
      {
        arma::mat &GmatN = Gamma_Exchange(NEUTRON, m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
        arma::mat &GmatP = Gamma_Exchange(PROTON , m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GmatN;
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GmatP;

        arma::mat &DmatN = Delta(NEUTRON, m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
        arma::mat &DmatP = Delta(PROTON , m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) = DmatN;
        field(PROTON , PAIRING).submat(bras.filter, kets.filter) = DmatP;
      }//if

      if (m_alpha == m_gamma + 1 && m_alpha > 0 && s_alpha == 1 && s_gamma == 0) //Mean exchange part -+ and Pairing part -+.
      {
        arma::mat &GmatN = Gamma_Exchange(NEUTRON, m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
        arma::mat &GmatP = Gamma_Exchange(PROTON , m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GmatN;
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GmatP;

        arma::mat &DmatN = Delta(NEUTRON, m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
        arma::mat &DmatP = Delta(PROTON , m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);

        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) = DmatN;
        field(PROTON , PAIRING).submat(bras.filter, kets.filter) = DmatP;
      }//if

      if (m_alpha == m_gamma - 1 && m_alpha < M - 1 && s_alpha == 0 && s_gamma == 1) //Mean exchange part +- and Pairing part +-.
      {
        arma::mat &GmatN = Gamma_Exchange(NEUTRON, m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);
        arma::mat &GmatP = Gamma_Exchange(PROTON , m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);

        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GmatN.t();
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GmatP.t();

        arma::mat &DmatN = Delta(NEUTRON, m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);
        arma::mat &DmatP = Delta(PROTON , m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);

        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) = DmatN.t();
        field(PROTON , PAIRING).submat(bras.filter, kets.filter) = DmatP.t();
      }//if
    }//j
  }//i

  // do we have symmetric matrices ?
  Tools::checkSymmetry(field(NEUTRON, PAIRING), "Pairing_Neutron");
  Tools::checkSymmetry(field(PROTON , PAIRING), "Pairing_Proton ");

  Tools::checkSymmetry(field(NEUTRON, DIRECT), "Mean_Central_Exchange_Neutron");
  Tools::checkSymmetry(field(PROTON , DIRECT), "Mean_Central_Exchange_Proton ");

  Tools::checkSymmetry(field(NEUTRON, EXCHANGE), "Mean_Central_Direct_Neutron");
  Tools::checkSymmetry(field(PROTON , EXCHANGE), "Mean_Central_Direct_Proton ");

  // Store the calculating length
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Ir for the direct and exchange parts.
 *  Its last element is an arma::vector and not an arma::cube because the n's dependence on m makes it ineffective.
 */

void FieldCentral::calcIr(void)
{
  DBG_ENTER;

  // dependencies
  calcJr();

  arma::wall_clock timer;
  timer.tic();

  if (!IrDirect.empty() && !IrExchange.empty()) DBG_LEAVE;

  INT M = basis.mMax;
  arma::mat MatDirect;
  arma::mat MatExchange;
  arma::mat MatExchange_sym;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates indices for the direct part.
  INT dimDirect  = 0;
  IVEC indicesDirect(M + 1);

  for (INT m = 0 ; m < M ; m++)
  {
    dimDirect += basis.nMax(m) * basis.nMax(m);
    indicesDirect(m + 1) = dimDirect;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////Calculates indices for the exchange part.
  INT dimExchange_0 = 0;
  IVEC indicesExchange_0(2 * M);

  for (INT m = -M + 1 ; m < M ; m++)
  {
    dimExchange_0 += basis.nMax(abs(m)) * basis.nMax(abs(m));
    indicesExchange_0(m + M) = dimExchange_0;
  }

  INT dimExchange_1 = 0;
  IVEC indicesExchange_1(2 * M - 1);

  for (INT m = -M + 2 ; m < M ; m++)
  {
    dimExchange_1 += basis.nMax(abs(m)) * basis.nMax(abs(m - 1));
    indicesExchange_1(m + M - 1) = dimExchange_1;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////Calculates IrDirect.
  for (int m = 0 ; m < M ; m++)
  {
    for (int n_alpha = 0 ; n_alpha < basis.nMax(abs(m)) ; n_alpha++)
    {
      for (int n_gamma = 0 ; n_gamma < basis.nMax(abs(m)) ; n_gamma++)
      {
        arma::vec &IrDvec = IrDirect(m, n_alpha, n_gamma);
        IrDvec = arma::zeros(dimDirect);
        arma::vec &jvec = Jr(m, n_alpha, m, n_gamma);

        for (int mp = 0 ; mp < M ; mp++)
        {
          MatDirect = arma::zeros(basis.nMax(mp), basis.nMax(mp));

          for (int n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
          {
            for (int n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
            {
              MatDirect(n_delta, n_beta) = ACCU(basis.talmanr(mp, n_beta, mp, n_delta) % jvec);
            }//n_delta
          }//n_beta

          IrDvec.subvec(indicesDirect(mp), indicesDirect(mp + 1) - 1) = arma::vectorise(MatDirect);
        }//mp
      }//n_gamma
    }//n_alpha
  }//m

  /////////////////////////////////////////////////////////////////////////////////////////Calculates IrExchange for m_alpha = m_gamma.
  for (int m = 0 ; m < M  ; m++)
  {
    for (int n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (int n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //One uses the symmetry gamma <-> alpha.
      {
        arma::vec &IrEvec =  IrExchange(0, m, n_alpha, n_gamma);
        arma::vec &IrEvec_sym =  IrExchange(0, m, n_gamma, n_alpha);
        IrEvec = arma::zeros(dimExchange_0);
        IrEvec_sym = arma::zeros(dimExchange_0);

        for (int mp = -M + 1 ; mp < M; mp++)
        {
          MatExchange = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp)));
          MatExchange_sym = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp)));

          for (int n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
          {
            arma::vec &jvec = Jr(m, n_alpha, mp, n_delta);

            for (int n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
            {
              MatExchange(n_delta, n_beta) = ACCU(basis.talmanr(mp, n_beta, m, n_gamma) % jvec);
              MatExchange_sym(n_beta, n_delta) = MatExchange(n_delta, n_beta);
            }// n_belta
          }// n_delta

          IrEvec.subvec(indicesExchange_0(mp - 1 + M), indicesExchange_0(mp + M) - 1) = arma::vectorise(MatExchange);
          IrEvec_sym.subvec(indicesExchange_0(mp - 1 + M), indicesExchange_0(mp + M) - 1) = arma::vectorise(MatExchange_sym);
        }//mp
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////Calculates IrExchange for m_alpha = m_gamma + 1.
  for (int m = 1 ; m < M  ; m++)
  {
    for (int n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (int n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        arma::vec &IrEvec =  IrExchange(1, m, n_alpha, n_gamma);
        IrEvec = arma::zeros(dimExchange_1);

        for (int mp = -M + 2 ; mp < M; mp++)
        {
          MatExchange = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp - 1)));

          for (int n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
          {
            arma::vec &jvec = Jr(m, n_alpha, mp, n_delta);

            for (int n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
            {
              MatExchange(n_delta, n_beta) = ACCU(basis.talmanr(mp - 1, n_beta, m - 1, n_gamma) % jvec);
            }//n_belta
          }//n_delta

          IrEvec.subvec(indicesExchange_1(mp - 2 + M), indicesExchange_1(mp + M - 1) - 1) = arma::vectorise(MatExchange);
        }//mp
      }//n_gamma
    }//n_alpha
  }//m

  ramTable["preIr"] =
    IrDirect.size()
  + IrExchange.size();

  timeTable["preIr"] = timer.toc();

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Iz for the direct and exchange parts.
 */

void FieldCentral::calcIz(void)
{
  DBG_ENTER;

  // dependencies
  calcJz();

  arma::wall_clock timer;
  timer.tic();

  if (!IzDirect.empty() && !IzExchange.empty()) DBG_LEAVE;

  INT dMax = basis.dMax;
  INT DMax = (dMax - 1) * 3 + 1;;
  INT n_zGMax = basis.n_zGlobalMax;
  INT Dd = 0;
  INT De = 0;
  INT De_sym = 0;
  INT Dd_sym = 0;

  /////////////////////////////////////////////////////////////Calculates IzDirect.
  for (int d_alpha = 0; d_alpha < dMax; d_alpha++)
  {
    for (int d_gamma = 0; d_gamma < dMax; d_gamma++)
    {
      for (int n_zalpha = 0; n_zalpha < n_zGMax; n_zalpha++)
      {
        for (int n_zgamma = 0; n_zgamma < n_zalpha + 1; n_zgamma++) //One uses the symmetry gamma <-> alpha.
        {
          IzDirect(d_alpha, d_gamma, n_zalpha, n_zgamma) = arma::zeros(n_zGMax, n_zGMax, DMax);
          IzDirect(d_gamma, d_alpha, n_zgamma, n_zalpha) = arma::zeros(n_zGMax, n_zGMax, DMax);
          arma::cube &IzCubeD_up = IzDirect(d_alpha, d_gamma, n_zalpha, n_zgamma);
          arma::cube &IzCubeD_low = IzDirect(d_gamma, d_alpha, n_zgamma, n_zalpha);
          arma::vec talvec = basis.talmanz(n_zalpha, d_alpha, n_zgamma, d_gamma).subvec(0, 2 * n_zGMax - 1);

          for (int d_delta = 0; d_delta < dMax; d_delta++)
          {
            INT D_gade = d_gamma + 2 * d_delta;

            for (int n_zdelta = 0; n_zdelta < n_zGMax; n_zdelta++)
            {
              for (INT d_beta = 0; d_beta < dMax; d_beta++)
              {
                Dd = 2 * d_beta + d_delta;
                De = 2 * d_beta + d_gamma;
                Dd_sym = 2 * d_delta + d_beta;
                INT D_albe = d_alpha + 2 * d_beta;

                for (int n_zbeta = 0; n_zbeta < n_zdelta + 1; n_zbeta++)
                {
                  IzCubeD_up(n_zdelta, n_zbeta, Dd) = ACCU(talvec % Jz(D_albe, D_gade, n_zbeta, n_zdelta));
                  IzCubeD_up(n_zbeta, n_zdelta, Dd_sym) = IzCubeD_up(n_zdelta, n_zbeta, Dd);
                  IzCubeD_low(n_zdelta, n_zbeta, Dd) = IzCubeD_up(n_zdelta, n_zbeta, Dd);
                  IzCubeD_low(n_zbeta, n_zdelta, Dd_sym) = IzCubeD_up(n_zdelta, n_zbeta, Dd);
                }//n_zbeta
              }//d_beta
            }//n_zdelta
          }//d_delta
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  /////////////////////////////////////////////////////////////Calculates IzExchange.
  for (int d_alpha = 0; d_alpha < dMax; d_alpha++)
  {
    for (int d_delta = 0; d_delta < dMax; d_delta++)
    {
      for (int n_zalpha = 0; n_zalpha < n_zGMax; n_zalpha++)
      {
        for (int n_zdelta = 0; n_zdelta < n_zalpha + 1; n_zdelta++) //One uses the symmetry delta <-> alpha.
        {
          IzExchange(d_alpha, d_delta, n_zalpha, n_zdelta) = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax, DMax);
          IzExchange(d_delta, d_alpha, n_zdelta, n_zalpha) = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax, DMax);
          arma::cube &IzCubeE = IzExchange(d_alpha, d_delta, n_zalpha, n_zdelta);
          arma::cube &IzCubeE_sym = IzExchange(d_delta, d_alpha, n_zdelta, n_zalpha);

          for (int d_gamma = 0; d_gamma < dMax; d_gamma++)
          {
            for (int d_beta = 0; d_beta < dMax; d_beta++)
            {
              Dd     = 2 * d_beta + d_delta;
              De     = 2 * d_beta + d_gamma;
              De_sym = 2 * d_gamma + d_beta;

              for (int n_zgamma = 0; n_zgamma < n_zGMax; n_zgamma++)
              {
                for (int n_zbeta = 0; n_zbeta < n_zGMax; n_zbeta++)
                {
                  IzCubeE(n_zgamma, n_zbeta, De) = IzDirect(d_alpha, d_gamma, n_zalpha, n_zgamma)(n_zdelta, n_zbeta, Dd);
                  IzCubeE_sym(n_zbeta, n_zgamma, De_sym) = IzCubeE(n_zgamma, n_zbeta, De);
                }//n_zbeta
              }//n_zgamma
            }//d_beta
          }//d_gamma
        }//n_zdelta
      }//n_zalpha
    }//d_delta
  }//d_alpha

  ramTable["preIz"] =
    IzDirect.size()
  + IzExchange.size();
  timeTable["preIz"] = timer.toc();

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Jr.
 */

void FieldCentral::calcJr(void)
{
  DBG_ENTER;

  // dependencies
  basis.calcTalmanr();
  basis.calcMoshinskyr();

  arma::wall_clock timer;
  timer.tic();

  if (!Jr.empty()) DBG_LEAVE;

  INT M = basis.mMax;
  INT maxn = basis.Nmaxr;
  INT nl = 0.0;
  double fac = 0.0;

  /////////////////////////////////////////////////////////////Calculates the Gaussian integral factor and merges it with a Moshinsky.
  Multi<arma::vec> mosGauss;

  fac =  pow(parameters["p"] / basis.b_r, 2);

  for (int ma = M * -2 + 2; ma < 2 * M - 1; ma++)
  {
    for (int nb = 0; nb < maxn; nb++)
    {
      mosGauss(ma, nb) = arma::zeros(maxn);
      arma::vec &mos = mosGauss(ma, nb);

      for (int na = 0; na < maxn; na++)
      {
        nl = na + nb + abs(ma);
        mos(na) = basis.moshinskyr(ma, nb)(na) * fac * (pow(2.0, nl) / pow(fac + 2.0, nl + 1));
      }//na
    }//nb
  }//ma

  /////////////////////////////////////////////////////////////Calculates Jr.
  for (int m_alpha = 0; m_alpha < M; m_alpha++)
  {
    for (int n_alpha = 0; n_alpha < basis.nMax(m_alpha); n_alpha++)
    {
      for (int m_beta = -M + 1; m_beta < M; m_beta++)
      {
        for (int n_beta = 0; n_beta < basis.nMax(abs(m_beta)); n_beta++)
        {
          Jr(m_alpha, n_alpha, m_beta, n_beta) = arma::zeros(maxn);
          Jr(-m_alpha, n_alpha, -m_beta, n_beta) = arma::zeros(maxn);

          arma::vec &talvec  = basis.talmanr(m_alpha, n_alpha, m_beta, n_beta);
          arma::mat &jvec  = Jr(m_alpha, n_alpha, m_beta, n_beta);
          arma::mat &jvec_min = Jr(-m_alpha, n_alpha, -m_beta, n_beta);

          for (int nb = 0 ; nb < maxn ; nb++)
          {
            jvec(nb) = ACCU(talvec % mosGauss(m_beta - m_alpha, nb));
            jvec_min(nb) = jvec(nb);
          }//nb
        }//n_beta
      }//m_beta
    }//n_alpha
  }//m_alpha

  ramTable["preJr"] = Jr.size();
  timeTable["preJr"] = timer.toc();

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Jz.
 */

void FieldCentral::calcJz(void)
{
  DBG_ENTER;

  // dependencies
  basis.calcTalmanz();
  basis.calcMoshinskyz();

  arma::wall_clock timer;
  timer.tic();

  if (!Jz.empty()) DBG_LEAVE;

  INT dMax = basis.dMax;
  INT n_zGlobalMax = basis.n_zGlobalMax;

  /////////////////////////////////////////////////////////////Gaussian integral factor's calculation.
  Multi<arma::vec> gaussianIntegralz;
  double fac0 = sqrt(sqrt(PI) * basis.b_z);

  double fac1 = parameters["p"] / basis.b_z;

  for (int d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
  {
    double zd_alpha = (0.5 - d_alpha) * basis.d_0;

    for (int d_beta = 0; d_beta < basis.dMax; d_beta++)
    {
      double zd_beta = (0.5 - d_beta) * basis.d_0;

      for (int d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
      {
        double zd_gamma = (0.5 - d_gamma) * basis.d_0;
        double d_alga = (zd_alpha + zd_gamma) / 2.0;

        for (int d_delta = 0; d_delta < basis.dMax; d_delta++)
        {
          double zd_delta = (0.5 - d_delta) * basis.d_0;
          double d_bede = (zd_beta + zd_delta) / 2.0;
          gaussianIntegralz(d_alpha, d_beta, d_gamma, d_delta) = arma::zeros(basis.n_zGlobalMax * 4);
          arma::vec &gauvec = gaussianIntegralz(d_alpha, d_beta, d_gamma, d_delta);

          for (int n_z = 0; n_z < basis.n_zGlobalMax * 4 - 3; n_z++)
          {
            double fac2 = sqrt(pow(2.0, n_z) / pow(2.0 + pow(fac1, 2), (n_z + 1)));
            double sigma = (d_bede - d_alga) / sqrt(2.0 + pow(fac1, 2));
            gauvec(n_z) = fac0 * fac1 * fac2 * exp(-0.5 * pow(sigma / basis.b_z, 2)) * basis.zPartScalar(sigma, n_z);
          }//n_z
        }//d_delta
      }//d_gamma
    }//d_beta
  }//d_alpha

  //////////////////////////////////////////////////////////////////////Merges a Moshinsky with the Gaussian integral factor.
  Multi<arma::vec> mosvec;

  for (int d_alpha = 0; d_alpha < dMax; d_alpha++)
  {
    for (int d_beta = 0; d_beta < dMax; d_beta++)
    {
      INT D_albe = d_alpha + 2 * d_beta;

      for (int d_gamma = 0; d_gamma < dMax; d_gamma++)
      {
        for (int d_delta = 0; d_delta < dMax; d_delta++)
        {
          INT D_gade = d_gamma + 2 * d_delta;
          arma::vec &gausVec = gaussianIntegralz(d_alpha, d_beta, d_gamma, d_delta);

          for (int n_za = 0; n_za < 2 * n_zGlobalMax; n_za ++)
          {
            mosvec(D_albe, D_gade, n_za) = arma::zeros(2 * n_zGlobalMax);
            arma::vec &mos = mosvec(D_albe, D_gade, n_za);

            for (int n_zb = 0; n_zb < 2 * n_zGlobalMax; n_zb ++)
            {
              mos(n_zb) = basis.moshinskyz(n_za, n_zb) * gausVec(n_za + n_zb);
            }//n_zb
          }//n_za
        }//d_delta
      }//d_gamma
    }//d_beta
  }//d_alpha

  //////////////////////////////////////////////////////////////////Calculates Jz.
  for (INT d_beta = 0; d_beta < dMax; d_beta++)
  {
    for (INT d_delta = 0; d_delta < dMax; d_delta++)
    {
      for (INT n_zbeta = 0; n_zbeta < n_zGlobalMax; n_zbeta++)
      {
        for (INT n_zdelta = 0; n_zdelta < n_zbeta + 1; n_zdelta++) //One uses the symmetry beta <-> delta.
        {
          arma::vec talvec = basis.talmanz(n_zbeta, d_beta, n_zdelta, d_delta).subvec(0, 2 * n_zGlobalMax - 1);

          for (int d_alpha = 0; d_alpha < dMax; d_alpha++)
          {
            INT D_albe = d_alpha + 2 * d_beta;
            INT D_albe_sym = d_alpha + 2 * d_delta;

            for (int d_gamma = 0; d_gamma < dMax; d_gamma++)
            {
              INT D_gade = d_gamma + 2 * d_delta;
              INT D_gade_sym = d_gamma + 2 * d_beta;
              Jz(D_albe, D_gade, n_zbeta, n_zdelta) = arma::zeros(basis.n_zGlobalMax * 2);
              Jz(D_albe_sym, D_gade_sym, n_zdelta, n_zbeta) = arma::zeros(basis.n_zGlobalMax * 2);
              arma::vec &jzvec_up = Jz(D_albe, D_gade, n_zbeta, n_zdelta);
              arma::vec &jzvec_low = Jz(D_albe_sym, D_gade_sym, n_zdelta, n_zbeta);

              for (int n_za = 0; n_za < n_zGlobalMax * 2; n_za++)
              {
                jzvec_up(n_za) = ACCU(talvec % mosvec(D_albe, D_gade, n_za));
                jzvec_low(n_za) = jzvec_up(n_za);
              }//n_za
            }//d_gamma
          }//d_alpha
        }//n_zdelta
      }//n_zbeta
    }//d_delta
  }//d_beta

  ramTable["preJz"] = Jz.size();
  timeTable["preJz"] = timer.toc();

  DBG_LEAVE;
}


