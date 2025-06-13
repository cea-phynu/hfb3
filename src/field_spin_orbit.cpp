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

#include "field_spin_orbit.h"
#include "tools.h"
#include "basis.h"
#include "interaction.h"
#include "discrete.h"
#include "qnumbers.h"

/** \file
 *  \brief Methods of the FieldSpinOrbit class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldSpinOrbit::FieldSpinOrbit(Field::Parameters fp, State *_state) : Field(fp, _state)
{
  DBG_ENTER;

  name = "spin-orbit";
  shortName = "SpiOrb";

  nGLA = 40; //Initially 20.
  nGHE = 100; //Initially 50.

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcField(void)
{
  DBG_ENTER;

  if (!field(NEUTRON, DIRECT).empty() && !field(PROTON, DIRECT).empty()) DBG_LEAVE;

  //dependencies
  calcXz();
  calcYz();
  calcYr();
  calcXr();

  warningStr = "";

  calculatingLength = -1.0;
  double startTime = Tools::clock();

  //Usefull quantities
  double wls = parameters["wso"];
  INT n_zMax = basis.n_zGlobalMax;
  INT Nmax = INT(basis.talmanr(0, 0, 0, 0).n_elem);
  INT dMax = basis.dMax;
  INT DMax = (dMax - 1) * 3 + 1;
  INT mMax = basis.mMax;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Kappa and Rho Multis with explicit quantum numbers.
  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});
  Multi<arma::mat> Rho; //Rho  (Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)
  //Multi<arma::mat,7> Kappa; //Kappa(Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)

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
        Rho(NEUTRON, m, mp, n, np, d + 2 * dp, s + 2 * sp) = rho(NEUTRON).submat(bras.filter, kets.filter);
        Rho(PROTON, m, mp, n, np, d + 2 * dp, s + 2 * sp) = rho(PROTON ).submat(bras.filter, kets.filter);
        //Kappa(NEUTRON,m,mp,n,np,d+2*dp,s+2*sp) = kappa(NEUTRON).submat(bras.filter, kets.filter);
        //Kappa(PROTON ,m,mp,n,np,d+2*dp,s+2*sp) = kappap.submat(bras.filter, kets.filter);
      }//if
    }//j
  }//i

  Multi<arma::cube> Rmat; //(index(3),Isospin,m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)
  //Multi<arma::cube,5> Kmat;//(index,Isospin,m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)

  for (INT iso: {NEUTRON, PROTON})
  {
    //Defines isospin and anti-isospin.
    INT non_iso = 0;

    if (iso == 0)
    {
      non_iso = 1;
    }

    //if mp = 0.
    for (INT n_beta = 0 ; n_beta < basis.nMax(0) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
      {
        arma::cube &Rc2 = Rmat(2, iso, 0, n_delta, n_beta);
        Rc2 = arma::zeros(basis.n_zMax(0, n_delta), basis.n_zMax(0, n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            Rc2.slice(d_delta + 2 * d_beta) = wls * (Rho(non_iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0) + 2 * Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0));
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta

    for (INT n_beta = 0 ; n_beta < basis.nMax(1) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
      {
        arma::cube &Rc0 = Rmat(0, iso, 0, n_delta, n_beta);
        Rc0 = arma::zeros(basis.n_zMax(0, n_delta), basis.n_zMax(1, n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            Rc0.slice(d_delta + 2 * d_beta) = 2 * sqrt(2) * wls * (Rho(non_iso, 0, 1, n_delta, n_beta, d_delta + 2 * d_beta, 2) + 2 * Rho(iso, 0, 1, n_delta, n_beta, d_delta + 2 * d_beta, 2));
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta

    for (INT mp = 1 ; mp < mMax - 1 ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::cube &Rc1 = Rmat(1, iso, mp, n_delta, n_beta);
          arma::cube &Rc2 = Rmat(2, iso, mp, n_delta, n_beta);
          Rc1 = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);
          Rc2 = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);

          for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
          {
            for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
            {
              Rc1.slice(d_delta + 2 * d_beta) = wls * (Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) - Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3)
                                                + 2 * (Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) - Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3)));
              Rc2.slice(d_delta + 2 * d_beta) = wls * (Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3)
                                                + 2 * (Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3)));
            }//d_delta
          }//d_beta
        }//n_delta
      }//n_beta

      for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::cube &Rc0 = Rmat(0, iso, mp, n_delta, n_beta);
          Rc0 = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp + 1, n_beta), DMax);

          for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
          {
            for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
            {
              Rc0.slice(d_delta + 2 * d_beta) = 2 * sqrt(2) * wls * (Rho(non_iso, mp, mp + 1, n_delta, n_beta, d_delta + 2 * d_beta, 2) + 2 * Rho(iso, mp, mp + 1, n_delta, n_beta, d_delta + 2 * d_beta, 2));
            }//d_delta
          }//d_beta
        }//n_delta
      }//n_beta
    }//mp

    //if mp = mMax -1.
    for (INT n_beta = 0 ; n_beta < basis.nMax(mMax - 1) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mMax - 1) ; n_delta++)
      {
        arma::cube &Rc1 = Rmat(1, iso, mMax - 1, n_delta, n_beta);
        arma::cube &Rc2 = Rmat(2, iso, mMax - 1, n_delta, n_beta);
        Rc1 = arma::zeros(basis.n_zMax(mMax - 1, n_delta), basis.n_zMax(mMax - 1, n_beta), DMax);
        Rc2 = arma::zeros(basis.n_zMax(mMax - 1, n_delta), basis.n_zMax(mMax - 1, n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            Rc1.slice(d_delta + 2 * d_beta) = wls * (Rho(non_iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 0) - Rho(non_iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 3)
                                              + 2 * (Rho(iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 0) - Rho(iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 3)));
            Rc2.slice(d_delta + 2 * d_beta) = wls * (Rho(non_iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(non_iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 3)
                                              + 2 * (Rho(iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(iso, mMax - 1, mMax - 1, n_delta, n_beta, d_delta + 2 * d_beta, 3)));
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with Xz.
  //One first considers n_beta <-> n_delta non-symmetric terms.
  // Multi<arma::mat> Rmatz; //(index(4),iso,n_za,mp,d_alpha+2*d_gamma)(n_delta,n_beta)

  FMulti<arma::mat> Rmatz({{0, 4}, TRANGE, {0, 2*basis.n_zGlobalMax+2}, MRANGE, DDRANGE});

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < d_alpha + 1 ; d_gamma++) //symmetry d_alpha <-> d_gamma.
      {
        for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
        {
          arma::cube &xz_0 = Xz(0, n_za, d_alpha + 2 * d_gamma);
          arma::cube &xz_1 = Xz(1, n_za, d_alpha + 2 * d_gamma);
          //if mp = 0.
          arma::mat &rz_00 = Rmatz(0, iso, n_za, 0, d_alpha + 2 * d_gamma);
          arma::mat &rz_10 = Rmatz(1, iso, n_za, 0, d_alpha + 2 * d_gamma);
          rz_00 = arma::zeros(basis.nMax(0), basis.nMax(1));
          rz_10 = arma::zeros(basis.nMax(0), basis.nMax(1));

          for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < basis.nMax(1) ; n_beta++)
            {
              rz_00(n_delta, n_beta) = arma::accu(xz_0.subcube(0, 0, 0, basis.n_zMax(0, n_delta) - 1, basis.n_zMax(1, n_beta) - 1, DMax - 1) % Rmat(0, iso, 0, n_delta, n_beta));
              rz_10(n_delta, n_beta) = arma::accu(xz_1.subcube(0, 0, 0, basis.n_zMax(0, n_delta) - 1, basis.n_zMax(1, n_beta) - 1, DMax - 1) % Rmat(0, iso, 0, n_delta, n_beta));
            }//n_beta
          }//n_delta

          Rmatz(0, iso, n_za, 0, d_gamma + 2 * d_alpha) = rz_00;
          Rmatz(1, iso, n_za, 0, d_gamma + 2 * d_alpha) = rz_10;

          for (INT mp = 1 ; mp < mMax - 1 ; mp++)
          {
            arma::mat &rz_0 = Rmatz(0, iso, n_za, mp, d_alpha + 2 * d_gamma);
            arma::mat &rz_1 = Rmatz(1, iso, n_za, mp, d_alpha + 2 * d_gamma);
            rz_0 = arma::zeros(basis.nMax(mp), basis.nMax(mp + 1));
            rz_1 = arma::zeros(basis.nMax(mp), basis.nMax(mp + 1));

            for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
            {
              for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
              {
                rz_0(n_delta, n_beta) = arma::accu(xz_0.subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp + 1, n_beta) - 1, DMax - 1) % Rmat(0, iso, mp, n_delta, n_beta));
                rz_1(n_delta, n_beta) = arma::accu(xz_1.subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp + 1, n_beta) - 1, DMax - 1) % Rmat(0, iso, mp, n_delta, n_beta));
              }//n_beta
            }//n_delta

            Rmatz(0, iso, n_za, mp, d_gamma + 2 * d_alpha) = rz_0;
            Rmatz(1, iso, n_za, mp, d_gamma + 2 * d_alpha) = rz_1;
          }//mp
        }//n_za
      }//d_gamma
    }//d_alpha
  }//iso


  //One considers n_beta <-> n_delta symmetric terms.
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < d_alpha + 1 ; d_gamma++) //symmetry d_alpha <-> d_gamma.
      {
        for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
        {
          arma::cube &xz_2 = Xz(2, n_za, d_alpha + 2 * d_gamma);
          //if mp = 0.
          arma::mat &rz_30 = Rmatz(3, iso, n_za, 0, d_alpha + 2 * d_gamma);
          rz_30 = arma::zeros(basis.nMax(0), basis.nMax(0));

          for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++) //symmetry n_beta <-> n_delta.
            {
              rz_30(n_delta, n_beta) = arma::accu(xz_2.subcube(0, 0, 0, basis.n_zMax(0, n_delta) - 1, basis.n_zMax(0, n_beta) - 1, DMax - 1) % Rmat(2, iso, 0, n_delta, n_beta));
              rz_30(n_beta, n_delta) = rz_30(n_delta, n_beta);
            }//n_beta
          }//n_delta

          Rmatz(3, iso, n_za, 0, d_gamma + 2 * d_alpha) = rz_30;

          for (INT mp = 1 ; mp < mMax - 1 ; mp++)
          {
            arma::mat &rz_2 = Rmatz(2, iso, n_za, mp, d_alpha + 2 * d_gamma);
            arma::mat &rz_3 = Rmatz(3, iso, n_za, mp, d_alpha + 2 * d_gamma);
            rz_2 = arma::zeros(basis.nMax(mp), basis.nMax(mp));
            rz_3 = arma::zeros(basis.nMax(mp), basis.nMax(mp));

            for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
            {
              for (INT n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++) //symmetry n_beta <-> n_delta.
              {
                rz_2(n_delta, n_beta) = arma::accu(xz_2.subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_beta) - 1, DMax - 1) % Rmat(1, iso, mp, n_delta, n_beta));
                rz_3(n_delta, n_beta) = arma::accu(xz_2.subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_beta) - 1, DMax - 1) % Rmat(2, iso, mp, n_delta, n_beta));
                rz_2(n_beta, n_delta) = rz_2(n_delta, n_beta);
                rz_3(n_beta, n_delta) = rz_3(n_delta, n_beta);
              }//n_beta
            }//n_delta

            Rmatz(2, iso, n_za, mp, d_gamma + 2 * d_alpha) = rz_2;
            Rmatz(3, iso, n_za, mp, d_gamma + 2 * d_alpha) = rz_3;
          }//mp

          //if mp = mMax-1;
          arma::mat &rz_2M = Rmatz(2, iso, n_za, mMax - 1, d_alpha + 2 * d_gamma);
          arma::mat &rz_3M = Rmatz(3, iso, n_za, mMax - 1, d_alpha + 2 * d_gamma);
          rz_2M = arma::zeros(basis.nMax(mMax - 1), basis.nMax(mMax - 1));
          rz_3M = arma::zeros(basis.nMax(mMax - 1), basis.nMax(mMax - 1));

          for (INT n_delta = 0 ; n_delta < basis.nMax(mMax - 1) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++) //symmetry n_beta <-> n_delta.
            {
              rz_2M(n_delta, n_beta) = arma::accu(xz_2.subcube(0, 0, 0, basis.n_zMax(mMax - 1, n_delta) - 1, basis.n_zMax(mMax - 1, n_beta) - 1, DMax - 1) % Rmat(1, iso, mMax - 1, n_delta, n_beta));
              rz_3M(n_delta, n_beta) = arma::accu(xz_2.subcube(0, 0, 0, basis.n_zMax(mMax - 1, n_delta) - 1, basis.n_zMax(mMax - 1, n_beta) - 1, DMax - 1) % Rmat(2, iso, mMax - 1, n_delta, n_beta));
              rz_2M(n_beta, n_delta) = rz_2M(n_delta, n_beta);
              rz_3M(n_beta, n_delta) = rz_3M(n_delta, n_beta);
            }//n_beta
          }//n_delta

          Rmatz(2, iso, n_za, mMax - 1, d_gamma + 2 * d_alpha) = rz_2M;
          Rmatz(3, iso, n_za, mMax - 1, d_gamma + 2 * d_alpha) = rz_3M;
        }//n_za
      }//d_gamma
    }//d_alpha
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with Xr.
  // Multi<arma::vec> Rmatr; //(index(4),iso,n_za,d_alpha+2*d_gamma)(n_a)

  FMulti<arma::vec> Rmatr({{0, 4}, TRANGE, {0, 2*basis.n_zGlobalMax+2}, DDRANGE});

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT n_za = 0 ; n_za < 2 * n_zMax + 2; n_za++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < d_alpha + 1 ; d_gamma++) //symmetry d_alpha <-> d_gamma.
        {
          arma::vec &rr_0 = Rmatr(0, iso, n_za, d_alpha + 2 * d_gamma);
          arma::vec &rr_1 = Rmatr(1, iso, n_za, d_alpha + 2 * d_gamma);
          arma::vec &rr_2 = Rmatr(2, iso, n_za, d_alpha + 2 * d_gamma);
          arma::vec &rr_3 = Rmatr(3, iso, n_za, d_alpha + 2 * d_gamma);
          rr_0 = arma::zeros(Nmax);
          rr_1 = arma::zeros(Nmax);
          rr_2 = arma::zeros(Nmax);
          rr_3 = arma::zeros(Nmax);

          for (INT n_a = 0 ; n_a < Nmax ; n_a++)
          {
            double sum_0 = 0;
            double sum_1 = 0;
            double sum_2 = 0;
            double sum_3 = 0;
            //for mp = 0.
            sum_0 += arma::accu(Xr(0, n_a, 0) % Rmatz(0, iso, n_za, 0, d_alpha + 2 * d_gamma));
            sum_1 += arma::accu(Xr(1, n_a, 0) % Rmatz(1, iso, n_za, 0, d_alpha + 2 * d_gamma));
            sum_3 += arma::accu(Xr(3, n_a, 0) % Rmatz(3, iso, n_za, 0, d_alpha + 2 * d_gamma));
            //for mp = mMax - 1.
            sum_2 += arma::accu(Xr(2, n_a, mMax - 1) % Rmatz(2, iso, n_za, mMax - 1, d_alpha + 2 * d_gamma));
            sum_3 += arma::accu(Xr(3, n_a, mMax - 1) % Rmatz(3, iso, n_za, mMax - 1, d_alpha + 2 * d_gamma));

            for (INT mp = 1 ; mp < mMax - 1 ; mp++)
            {
              sum_0 += arma::accu(Xr(0, n_a, mp) % Rmatz(0, iso, n_za, mp, d_alpha + 2 * d_gamma));
              sum_1 += arma::accu(Xr(1, n_a, mp) % Rmatz(1, iso, n_za, mp, d_alpha + 2 * d_gamma));
              sum_2 += arma::accu(Xr(2, n_a, mp) % Rmatz(2, iso, n_za, mp, d_alpha + 2 * d_gamma));
              sum_3 += arma::accu(Xr(3, n_a, mp) % Rmatz(3, iso, n_za, mp, d_alpha + 2 * d_gamma));
            }//mp

            rr_0(n_a) = sum_0;
            rr_1(n_a) = sum_1;
            rr_2(n_a) = sum_2;
            rr_3(n_a) = sum_3;
          }//n_a

          Rmatr(0, iso, n_za, d_gamma + 2 * d_alpha) = rr_0;
          Rmatr(1, iso, n_za, d_gamma + 2 * d_alpha) = rr_1;
          Rmatr(2, iso, n_za, d_gamma + 2 * d_alpha) = rr_2;
          Rmatr(3, iso, n_za, d_gamma + 2 * d_alpha) = rr_3;
        }//d_gamma
      }//d_alpha
    }//n_za
  }//iso


  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with Yr.
  // Multi<arma::vec> Gmatr; //(index(5),iso,m,n_alpha,n_gamma,d_alpha+2*d_gamma)(n_za)

  FMulti<arma::vec> Gmatr({{0, 5}, TRANGE, MRANGE, NRANGE, NRANGE, DDRANGE});

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        //Gmatr(0) and Gmatr(1).
        for (INT m = 0 ; m < mMax ; m++)
        {
          for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
          {
            for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //symmetry n_alpha <-> n_gamma.
            {
              arma::vec &gr_0 = Gmatr(0, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              arma::vec &gr_1 = Gmatr(1, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              gr_0 = arma::zeros(2 * n_zMax + 2);
              gr_1 = arma::zeros(2 * n_zMax + 2);

              arma::vec &yr_0 = Yr(0, m, n_alpha, n_gamma);

              for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
              {
                arma::vec sumRmat_01 = Rmatr(0, iso, n_za, d_alpha + 2 * d_gamma) - Rmatr(1, iso, n_za, d_alpha + 2 * d_gamma);
                gr_0(n_za) = arma::accu(yr_0 % sumRmat_01);
                gr_1(n_za) = arma::accu(yr_0 % Rmatr(2, iso, n_za, d_alpha + 2 * d_gamma));
              }//n_za

              Gmatr(0, iso, m, n_gamma, n_alpha, d_alpha + 2 * d_gamma) = gr_0;
              Gmatr(1, iso, m, n_gamma, n_alpha, d_alpha + 2 * d_gamma) = gr_1;
            }//n_gamma
          }//n_alpha
        }//m

        //Gmatr(2).
        for (INT m = 1 ; m < mMax ; m++)
        {
          for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
          {
            for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //symmetry n_alpha <-> n_gamma.
            {
              arma::vec &gr_2 = Gmatr(2, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              gr_2 = arma::zeros(2 * n_zMax + 2);

              arma::vec &yr_1 = Yr(1, m, n_alpha, n_gamma);

              for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
              {
                gr_2(n_za) = arma::accu(yr_1 % Rmatr(3, iso, n_za, d_alpha + 2 * d_gamma));
              }//n_za

              Gmatr(2, iso, m, n_gamma, n_alpha, d_alpha + 2 * d_gamma) = gr_2;
            }//n_gamma
          }//n_alpha
        }//m

        //Gmatr(3) and Gmatr(4).
        for (INT m = 1 ; m < mMax ; m++)
        {
          for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
          {
            for (INT n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
            {
              arma::vec &gr_3 = Gmatr(3, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              arma::vec &gr_4 = Gmatr(4, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              gr_3 = arma::zeros(2 * n_zMax + 2);
              gr_4 = arma::zeros(2 * n_zMax + 2);

              arma::vec &yr_2 = Yr(2, m, n_alpha, n_gamma);
              arma::vec &yr_3 = Yr(3, m, n_alpha, n_gamma);

              for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
              {
                gr_3(n_za) = arma::accu(yr_2 % Rmatr(3, iso, n_za, d_alpha + 2 * d_gamma));
                gr_4(n_za) = arma::accu(yr_3 % Rmatr(3, iso, n_za, d_alpha + 2 * d_gamma));
              }//n_za
            }//n_gamma
          }//n_alpha
        }//m
      }//d_gamma
    }//d_alpha
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with Yz.
  // Multi<arma::cube> Gamma; //(iso,m,n_alpha,n_gamma,s_alpha+2*s_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma)

  FMulti<arma::cube> Gamma({TRANGE, MPRANGE, NRANGE, NRANGE, SSRANGE});


  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //symmetry n_alpha <-> n_gamma.
        {
          arma::cube &g_pp = Gamma(iso, m, n_alpha, n_gamma, 0);
          g_pp = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);

          for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
            {
              arma::vec  &gr_0 = Gmatr(0, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              arma::vec  &gr_1 = Gmatr(1, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);

              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  arma::vec &yz_0 = Yz(0, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
                  g_pp(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += arma::accu(yz_0 % gr_0);
                  g_pp(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += arma::accu(yz_0 % gr_1);
                }//n_zgamma
              }//n_zalpha
            }//d_alpha
          }//d_gamma

          Gamma(iso, m, n_alpha, n_gamma, 3) = g_pp;
        }//n_gamma
      }//d_alpha
    }//m

    for (INT m = 1 ; m < mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //symmetry n_alpha <-> n_gamma.
        {
          arma::cube &g_pp = Gamma(iso, m, n_alpha, n_gamma, 0);
          arma::cube &g_mm = Gamma(iso, m, n_alpha, n_gamma, 3);

          for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
            {
              arma::vec  &gr_2 = Gmatr(2, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);

              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  arma::vec &yz_0 = Yz(0, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
                  double value = arma::accu(yz_0 % gr_2);
                  g_pp(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) +=  value;
                  g_mm(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += -value;
                }//n_zgamma
              }//n_zalpha
            }//d_alpha
          }//d_gamma
        }//n_gamma
      }//d_alpha
    }//m

    for (INT m = 1 ; m < mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
        {
          arma::cube &g_pm = Gamma(iso, m, n_alpha, n_gamma, 1);
          g_pm = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);

          for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
            {
              arma::vec  &gr_3 = Gmatr(3, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);
              arma::vec  &gr_4 = Gmatr(4, iso, m, n_alpha, n_gamma, d_alpha + 2 * d_gamma);

              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m - 1, n_gamma) ; n_zgamma++)
                {
                  arma::vec &yz_1 = Yz(1, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
                  arma::vec &yz_2 = Yz(2, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
                  g_pm(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) +=   arma::accu(yz_1 % gr_3);
                  g_pm(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) +=   -arma::accu(yz_2 % gr_4);
                }//n_zgamma
              }//n_zalpha
            }//d_alpha
          }//d_gamma

          g_pm = g_pm * sqrt(2);
        }//n_gamma
      }//d_alpha
    }//m
  }//iso

  //Fills with symmetric terms.
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha ; n_gamma++) //symmetry n_alpha <-> n_gamma.
        {
          arma::cube &g_pp = Gamma(iso, m, n_alpha, n_gamma, 0);
          arma::cube &g_pp_sym = Gamma(iso, m, n_gamma, n_alpha, 0);
          g_pp_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

          arma::cube &g_mm = Gamma(iso, m, n_alpha, n_gamma, 3);
          arma::cube &g_mm_sym = Gamma(iso, m, n_gamma, n_alpha, 3);
          g_mm_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

          for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
            {
              g_pp_sym.slice(d_gamma + 2 * d_alpha) = (g_pp.slice(d_alpha + 2 * d_gamma)).t();
              g_mm_sym.slice(d_gamma + 2 * d_alpha) = (g_mm.slice(d_alpha + 2 * d_gamma)).t();
            }//d_alpha
          }//d_gamma
        }//n_gamma
      }//d_alpha
    }//m
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Gamma and Delta.
  field(NEUTRON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  for (UINT i = 0; i < nbBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT n_alpha =   bras(1, 0);
    INT d_alpha =   bras(3, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT n_gamma = kets(1, 0);
      INT d_gamma = kets(3, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0)//Mean Exchange part ++.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean Exchange part --.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma + 1 && m_alpha > 0 && s_alpha == 1 && s_gamma == 0) //Mean part -+.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma - 1 && m_alpha < mMax - 1 && s_alpha == 0 && s_gamma == 1) //Mean part +-.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = (Gamma(NEUTRON, m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha)).t();
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = (Gamma(PROTON, m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha)).t();
      }//if
    }//j
  }//i

  //Stores the calculating length
  calculatingLength = Tools::clock() - startTime;
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcXz(void)
{
  DBG_ENTER;

  if (!Xz.empty()) DBG_LEAVE;

  //dependencies
  basis.calcTalmanz();
  basis.calcMoshinskyz();

  //Usefull quantities
  INT n_zMax = basis.n_zGlobalMax;
  double d_0 = basis.d_0;
  INT dMax = basis.dMax;
  INT DMax = (basis.dMax - 1) * 3 + 1;
  double b_z = basis.b_z;
  double fac_xz01 = 1 / (2 * std::pow(basis.b_z, 1.5) * std::pow(PI, 0.25));
  double fac_xz2  = 1 / sqrt(2 * basis.b_z * sqrt(PI));

  //////////////////////////////////////////////////////////////////////////////////////////////////////Adds factors to the Moshinsky-z.
  Multi<arma::mat> addMosz;

  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    double dv_alpha = (d_alpha - 0.5) * d_0;

    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      double dv_gamma = (d_gamma - 0.5) * d_0;
      double k_ag = 0.5 * (dv_alpha + dv_gamma);

      for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
      {
        double dv_delta = (d_delta - 0.5) * d_0;

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          double dv_beta = (d_beta - 0.5) * d_0;
          double k_bd = 0.5 * (dv_beta + dv_delta);
          double d_agbd = (k_ag - k_bd) / sqrt(2);
          double fac_exp = exp(-0.5 * std::pow(d_agbd / b_z, 2));
          arma::mat &mosMatz = addMosz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta);
          mosMatz = basis.moshinskyz * fac_exp;

          for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
          {
            for (INT n_zb = 0 ; n_zb < 2 * n_zMax + 2 ; n_zb++)
            {
              mosMatz(n_za, n_zb) = mosMatz(n_za, n_zb) * basis.zPartScalar(d_agbd, n_za + n_zb, -1);
            }//n_zb
          }//n_za
        }//d_beta
      }//d_delta
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-z for Xz(0).
  Multi<arma::vec> newTz_0;

  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
      {
        //if n_zbeta = 0.
        newTz_0(n_zdelta, 0, d_delta + 2 * d_beta) = -fac_xz01 * basis.talmanz(n_zdelta, d_delta, 1, d_beta);

        for (INT n_zbeta = 1 ; n_zbeta < n_zMax ; n_zbeta++)
        {
          newTz_0(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = fac_xz01 * (sqrt(n_zbeta) * basis.talmanz(n_zdelta, d_delta, n_zbeta - 1, d_beta) - sqrt(n_zbeta + 1) * basis.talmanz(n_zdelta, d_delta, n_zbeta + 1, d_beta));
        }//n_zbeta
      }//n_zdelta
    }//d_beta
  }//d_delta

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-z for Xz(1).
  Multi<arma::vec> newTz_1;

  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
      {
        //if n_zdelta = 0.
        newTz_1(0, n_zbeta, d_delta + 2 * d_beta) = -fac_xz01 * basis.talmanz(1, d_delta, n_zbeta, d_beta);

        for (INT n_zdelta = 1 ; n_zdelta < n_zMax ; n_zdelta++)
        {
          newTz_1(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = fac_xz01 * (sqrt(n_zdelta) * basis.talmanz(n_zdelta - 1, d_delta, n_zbeta, d_beta) - sqrt(n_zdelta + 1) * basis.talmanz(n_zdelta + 1, d_delta, n_zbeta, d_beta));
        }//n_zdelta
      }//n_zbeta
    }//d_beta
  }//d_delta

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates no-reshaped form for Xz(0), Xz(1) and Xz(2).
  Multi<arma::vec> preXz_0;
  Multi<arma::vec> preXz_1;
  Multi<arma::vec> preXz_2;

  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
      {
        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          arma::mat &mosMatz = addMosz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta);

          for (INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
            {
              preXz_0(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta) = mosMatz * newTz_0(n_zdelta, n_zbeta, d_delta + 2 * d_beta);
              preXz_1(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta) = mosMatz * newTz_1(n_zdelta, n_zbeta, d_delta + 2 * d_beta);
              preXz_2(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta) = fac_xz2 * mosMatz * basis.talmanz(n_zdelta, d_delta, n_zbeta, d_beta);
            }//n_zbeta
          }//n_zdelta
        }//d_beta
      }//d_delta
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Reshapes all elements.
  for (INT n_za = 0 ; n_za < 2 * n_zMax + 2 ; n_za++)
  {
    for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
      {
        arma::cube &xzCube_0 = Xz(0, n_za, d_alpha + 2 * d_gamma);
        arma::cube &xzCube_1 = Xz(1, n_za, d_alpha + 2 * d_gamma);
        arma::cube &xzCube_2 = Xz(2, n_za, d_alpha + 2 * d_gamma);
        xzCube_0 = arma::zeros(n_zMax, n_zMax, DMax);
        xzCube_1 = arma::zeros(n_zMax, n_zMax, DMax);
        xzCube_2 = arma::zeros(n_zMax, n_zMax, DMax);

        for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
        {
          for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
          {
            for (INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
            {
              for (INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
              {
                xzCube_0(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = preXz_0(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta)(n_za);
                xzCube_1(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = preXz_1(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta)(n_za);
                xzCube_2(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = preXz_2(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma, d_delta + 2 * d_beta)(n_za);
              }//n_zbeta
            }//n_zdelta
          }//d_beta
        }//d_delta
      }//d_gamma
    }//d_alpha
  }//n_za

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcYz(void)
{
  DBG_ENTER;

  if (!Yz.empty()) DBG_LEAVE;

  //dependencies
  basis.calcTalmanz();
  basis.calcMoshinskyz();

  //Usefull quantities
  INT n_zMax = basis.n_zGlobalMax;
  INT dMax = basis.dMax;
  double fac_yz12 = 1 / (basis.b_z * sqrt(2));

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Yz(0).
  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
        {
          Yz(0, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = basis.talmanz(n_zalpha, d_alpha, n_zgamma, d_gamma);
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Yz(1).
  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT n_zgamma = 0 ; n_zgamma < n_zMax ; n_zgamma++)
      {
        //if n_zalpha = 0.
        Yz(1, 0, n_zgamma, d_alpha + 2 * d_gamma) = - fac_yz12 * basis.talmanz(1, d_alpha, n_zgamma, d_gamma);

        for (INT n_zalpha = 1 ; n_zalpha < n_zMax ; n_zalpha++)
        {
          Yz(1, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_yz12 * (sqrt(n_zalpha) * basis.talmanz(n_zalpha - 1, d_alpha, n_zgamma, d_gamma) - sqrt(n_zalpha + 1) * basis.talmanz(n_zalpha + 1, d_alpha, n_zgamma, d_gamma));
        }//n_zalpha
      }//n_zgamma
    }//d_gamma
  }//d_alpha


  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Yz(2).
  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT n_zalpha = 0 ; n_zalpha < n_zMax ; n_zalpha++)
      {
        //if n_zgamma = 0.
        Yz(2, n_zalpha, 0, d_alpha + 2 * d_gamma) = - fac_yz12 * basis.talmanz(n_zalpha, d_alpha, 1, d_gamma);

        for (INT n_zgamma = 1 ; n_zgamma < n_zMax ; n_zgamma++)
        {
          Yz(2, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_yz12 * (sqrt(n_zgamma) * basis.talmanz(n_zalpha, d_alpha, n_zgamma - 1, d_gamma) - sqrt(n_zgamma + 1) * basis.talmanz(n_zalpha, d_alpha, n_zgamma + 1, d_gamma));
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcXr(void)
{
  DBG_ENTER;

  if (!Xr.empty()) DBG_LEAVE;

  //dependencies
  basis.calcTalmanr();
  basis.calcMoshinskyr();

  //Usefull quantities
  INT Nmax = INT(basis.talmanr(0, 0, 0, 0).n_elem);
  INT mMax = basis.mMax;
  double b_r = basis.b_r;
  double fac_xr01 = 1 / (std::pow(b_r, 3) * PI * std::pow(2, 1.5));
  double fac_xr2  = 1 / (4 * std::pow(b_r, 4) * PI);

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Xr(0).
  //mp = m_delta = m_beta-1
  Multi<arma::vec> newTr_0;

  for (INT mp = 0 ; mp < mMax - 1 ; mp++)
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
    {
      //if n_delta = 0.
      newTr_0(mp, 0, n_beta) = fac_xr01 * sqrt(mp + 1) * basis.talmanr(mp + 1, 0, mp + 1, n_beta);

      for (INT n_delta = 1 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        newTr_0(mp, n_delta, n_beta) = fac_xr01 * ( sqrt(mp + n_delta + 1) * basis.talmanr(mp + 1, n_delta, mp + 1, n_beta) + sqrt(n_delta) * basis.talmanr(mp + 1, n_delta - 1, mp + 1, n_beta));
      }//n_delta
    }//n_beta
  }//mp

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Xr(1).
  //mp = m_delta = m_beta-1
  Multi<arma::vec> newTr_1;

  for (INT mp = 0 ; mp < mMax - 1 ; mp++)
  {
    for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
      {
        newTr_1(mp, n_delta, n_beta) = -fac_xr01 * (sqrt(mp + n_beta + 1) * basis.talmanr(mp, n_delta, mp, n_beta) + sqrt(n_beta + 1) * basis.talmanr(mp, n_delta, mp, n_beta + 1));
      }//n_beta
    }//n_delta
  }//mp

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Xr(2).
  //mp = m_delta = m_beta
  Multi<arma::vec> newTr_2;

  for (INT mp = 1 ; mp < mMax ; mp++)
  {
    //if n_delta = 0 and n_beta>0.
    for (INT n_beta = 1 ; n_beta < basis.nMax(mp) ; n_beta++)
    {
      newTr_2(mp, 0, n_beta) = fac_xr2 * (-sqrt((mp + n_beta) * mp) * basis.talmanr(mp - 1, 0, mp - 1, n_beta) - sqrt(mp + n_beta) * basis.talmanr(mp - 1, 1, mp - 1, n_beta)
                                          - sqrt((n_beta + 1) * mp) * basis.talmanr(mp - 1, 0, mp - 1, n_beta + 1) - sqrt(n_beta + 1) * basis.talmanr(mp - 1, 1, mp - 1, n_beta + 1)
                                          + sqrt((mp + n_beta + 1) * (mp + 1)) * basis.talmanr(mp + 1, 0, mp + 1, n_beta)
                                          + sqrt(n_beta * (mp + 1)) * basis.talmanr(mp + 1, 0, mp + 1, n_beta - 1));
    }//n_beta

    //if n_delta = 0 and n_beta = 0.
    newTr_2(mp, 0, 0) = fac_xr2 * (-sqrt(mp * mp) * basis.talmanr(mp - 1, 0, mp - 1, 0) - sqrt(mp) * basis.talmanr(mp - 1, 1, mp - 1, 0)
                                   - sqrt(mp) * basis.talmanr(mp - 1, 0, mp - 1, 1) - basis.talmanr(mp - 1, 1, mp - 1, 1)
                                   + sqrt((mp + 1) * (mp + 1)) * basis.talmanr(mp + 1, 0, mp + 1, 0));

    for (INT n_delta = 1 ; n_delta < basis.nMax(mp) ; n_delta++)
    {
      //if n_beta = 0.
      newTr_2(mp, n_delta, 0) = fac_xr2 * (-sqrt(mp * (mp + n_delta)) * basis.talmanr(mp - 1, n_delta, mp - 1, 0) - sqrt(mp * (n_delta + 1)) * basis.talmanr(mp - 1, n_delta + 1, mp - 1, 0)
                                           - sqrt(mp + n_delta) * basis.talmanr(mp - 1, n_delta, mp - 1, 1) - sqrt(n_delta + 1) * basis.talmanr(mp - 1, n_delta + 1, mp - 1, 1)
                                           + sqrt((mp + 1) * (mp + n_delta + 1)) * basis.talmanr(mp + 1, n_delta, mp + 1, 0) + sqrt((mp + 1) * n_delta) * basis.talmanr(mp + 1, n_delta - 1, mp + 1, 0));

      for (INT n_beta = 1 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        newTr_2(mp, n_delta, n_beta) = fac_xr2 * (-sqrt((mp + n_beta) * (mp + n_delta)) * basis.talmanr(mp - 1, n_delta, mp - 1, n_beta) - sqrt((mp + n_beta) * (n_delta + 1)) * basis.talmanr(mp - 1, n_delta + 1, mp - 1, n_beta)
                                       - sqrt((n_beta + 1) * (mp + n_delta)) * basis.talmanr(mp - 1, n_delta, mp - 1, n_beta + 1) - sqrt((n_beta + 1) * (n_delta + 1)) * basis.talmanr(mp - 1, n_delta + 1, mp - 1, n_beta + 1)
                                       + sqrt((mp + n_beta + 1) * (mp + n_delta + 1)) * basis.talmanr(mp + 1, n_delta, mp + 1, n_beta) + sqrt((mp + n_beta + 1) * n_delta) * basis.talmanr(mp + 1, n_delta - 1, mp + 1, n_beta)
                                       + sqrt(n_beta * (mp + n_delta + 1)) * basis.talmanr(mp + 1, n_delta, mp + 1, n_beta - 1) + sqrt(n_beta * n_delta) * basis.talmanr(mp + 1, n_delta - 1, mp + 1, n_beta - 1));
      }//n_beta
    }//n_delta
  }//mp

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates no-reshaped form for Xr(0), Xr(1) and Xr(2).
  Multi<arma::vec> preXr_0;
  Multi<arma::vec> preXr_1;
  Multi<arma::vec> preXr_2;

  for (INT mp = 0 ; mp < mMax - 1 ; mp++)
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        arma::vec &prexr_0 = preXr_0(mp, n_delta, n_beta);
        arma::vec &prexr_1 = preXr_1(mp, n_delta, n_beta);
        prexr_0 = arma::zeros(Nmax);
        prexr_1 = arma::zeros(Nmax);
        arma::vec &newtr_0 = newTr_0(mp, n_delta, n_beta);
        arma::vec &newtr_1 = newTr_1(mp, n_delta, n_beta);

        for (INT n_a = 0 ; n_a < Nmax ; n_a++)
        {
          prexr_0(n_a) = arma::accu(newtr_0 % basis.moshinskyr(0, n_a));
          prexr_1(n_a) = arma::accu(newtr_1 % basis.moshinskyr(0, n_a));
        }//n_a
      }//n_delta
    }//n_beta
  }//mp

  for (INT mp = 1 ; mp < mMax ; mp++)
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        arma::vec &prexr_2 = preXr_2(mp, n_delta, n_beta);
        prexr_2 = arma::zeros(Nmax);
        arma::vec &newtr_2 = newTr_2(mp, n_delta, n_beta);

        for (INT n_a = 0 ; n_a < Nmax ; n_a++)
        {
          prexr_2(n_a) = arma::accu(newtr_2 % basis.moshinskyr(0, n_a));
        }//n_a
      }//n_delta
    }//n_beta
  }//mp

  //////////////////////////////////////////////////////////////////////////////////////////////////////Reshapes all elements.
  for (INT n_a = 0 ; n_a < Nmax ; n_a++)
  {
    Xr(0, n_a, mMax - 1) = arma::zeros(0, 0);
    Xr(1, n_a, mMax - 1) = arma::zeros(0, 0);

    for (INT mp = 0 ; mp < mMax - 1  ; mp++)
    {
      arma::mat &xr_0 = Xr(0, n_a, mp);
      arma::mat &xr_1 = Xr(1, n_a, mp);
      xr_0 = arma::zeros(basis.nMax(mp), basis.nMax(mp + 1));
      xr_1 = arma::zeros(basis.nMax(mp), basis.nMax(mp + 1));

      for (INT n_beta = 0 ; n_beta < basis.nMax(mp + 1) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          xr_0(n_delta, n_beta) = preXr_0(mp, n_delta, n_beta)(n_a);
          xr_1(n_delta, n_beta) = preXr_1(mp, n_delta, n_beta)(n_a);
        }//n_delta
      }//n_beta
    }//mp
  }//n_a

  for (INT n_a = 0 ; n_a < Nmax ; n_a++)
  {
    Xr(2, n_a, 0) = arma::zeros(0, 0);

    for (INT mp = 1 ; mp < mMax ; mp++)
    {
      arma::mat &xr_2 = Xr(2, n_a, mp);
      xr_2 = arma::zeros(basis.nMax(mp), basis.nMax(mp));

      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          xr_2(n_delta, n_beta) = preXr_2(mp, n_delta, n_beta)(n_a);
        }//n_delta
      }//n_beta
    }//mp
  }//n_a

  for (INT n_a = 0 ; n_a < Nmax ; n_a++)
  {
    for (INT mp = 0 ; mp < mMax ; mp++)
    {
      arma::mat &xr_3 = Xr(3, n_a, mp);
      xr_3 = arma::zeros(basis.nMax(mp), basis.nMax(mp));

      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          xr_3(n_delta, n_beta) = basis.talmanr(mp, n_delta, mp, n_beta)(n_a);
        }//n_delta
      }//n_beta
    }//mp
  }//n_a

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcYr(void)
{
  DBG_ENTER;

  if (!Yr.empty()) DBG_LEAVE;

  //dependencies
  basis.calcTalmanr();
  basis.calcMoshinskyr();

  //Usefull quantities
  INT Nmax = INT(basis.talmanr(0, 0, 0, 0).n_elem);
  INT mMax = basis.mMax;
  double b_r = basis.b_r;
  double fac_yr23 = 1 / (std::pow(b_r, 3) * PI * std::pow(2, 1.5));
  double fac_yr1  = 1 / (4 * std::pow(b_r, 4) * PI);

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Yr(1).
  //m = m_alpha = m_gamma
  Multi<arma::vec> newTr_1;

  for (INT m = 1 ; m < mMax ; m++)
  {
    //if n_alpha = 0 and n_gamma>0.
    for (INT n_gamma = 1 ; n_gamma < basis.nMax(m) ; n_gamma++)
    {
      newTr_1(m, 0, n_gamma) = fac_yr1 * (-sqrt((m + n_gamma) * m) * basis.talmanr(m - 1, 0, m - 1, n_gamma) - sqrt(m * (n_gamma + 1)) * basis.talmanr(m - 1, 0, m - 1, n_gamma + 1)
                                          - sqrt(m + n_gamma) * basis.talmanr(m - 1, 1, m - 1, n_gamma) - sqrt(n_gamma + 1) * basis.talmanr(m - 1, 1, m - 1, n_gamma + 1)
                                          + sqrt((m + n_gamma + 1) * (m + 1)) * basis.talmanr(m + 1, 0, m + 1, n_gamma)
                                          + sqrt(n_gamma * (m + 1)) * basis.talmanr(m + 1, 0, m + 1, n_gamma - 1));
    }//n_gamma

    //if n_alpha = 0 and n_gamma = 0.
    newTr_1(m, 0, 0) = fac_yr1 * (-sqrt(m * m) * basis.talmanr(m - 1, 0, m - 1, 0) - sqrt(m) * basis.talmanr(m - 1, 0, m - 1, 1)
                                  - sqrt(m) * basis.talmanr(m - 1, 1, m - 1, 0) - basis.talmanr(m - 1, 1, m - 1, 1)
                                  + sqrt((m + 1) * (m + 1)) * basis.talmanr(m + 1, 0, m + 1, 0));

    for (INT n_alpha = 1 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      //if n_gamma = 0.
      newTr_1(m, n_alpha, 0) = fac_yr1 * (-sqrt(m * (m + n_alpha)) * basis.talmanr(m - 1, n_alpha, m - 1, 0) - sqrt(m + n_alpha) * basis.talmanr(m - 1, n_alpha, m - 1, 1)
                                          - sqrt((n_alpha + 1) * m) * basis.talmanr(m - 1, n_alpha + 1, m - 1, 0) - sqrt(n_alpha + 1) * basis.talmanr(m - 1, n_alpha + 1, m - 1, 1)
                                          + sqrt((m + 1) * (m + n_alpha + 1)) * basis.talmanr(m + 1, n_alpha, m + 1, 0) + sqrt((m + 1) * n_alpha) * basis.talmanr(m + 1, n_alpha - 1, m + 1, 0));

      for (INT n_gamma = 1 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        newTr_1(m, n_alpha, n_gamma) = fac_yr1 * (-sqrt((m + n_alpha) * (m + n_gamma)) * basis.talmanr(m - 1, n_alpha, m - 1, n_gamma) - sqrt((m + n_alpha) * (n_gamma + 1)) * basis.talmanr(m - 1, n_alpha, m - 1, n_gamma + 1)
                                       - sqrt((n_alpha + 1) * (m + n_gamma)) * basis.talmanr(m - 1, n_alpha + 1, m - 1, n_gamma) - sqrt((n_alpha + 1) * (n_gamma + 1)) * basis.talmanr(m - 1, n_alpha + 1, m - 1, n_gamma + 1)
                                       + sqrt((m + n_alpha + 1) * (m + n_gamma + 1)) * basis.talmanr(m + 1, n_alpha, m + 1, n_gamma) + sqrt((m + n_alpha + 1) * n_gamma) * basis.talmanr(m + 1, n_alpha, m + 1, n_gamma - 1)
                                       + sqrt(n_alpha * (m + n_gamma + 1)) * basis.talmanr(m + 1, n_alpha - 1, m + 1, n_gamma) + sqrt(n_alpha * n_gamma) * basis.talmanr(m + 1, n_alpha - 1, m + 1, n_gamma - 1));
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Yr(2).
  //m = m_alpha = m_gamma+1
  Multi<arma::vec> newTr_2;

  for (INT m = 1 ; m < mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      //if n_gamma = 0.
      newTr_2(m, n_alpha, 0) = fac_yr23 * sqrt(m) * basis.talmanr(m, n_alpha, m, 0);

      for (INT n_gamma = 1 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        newTr_2(m, n_alpha, n_gamma) = fac_yr23 * (sqrt(m + n_gamma) * basis.talmanr(m, n_alpha, m, n_gamma) + sqrt(n_gamma) * basis.talmanr(m, n_alpha, m, n_gamma - 1));
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates new Talman-r for Yr(3).
  //m = m_alpha = m_gamma+1
  Multi<arma::vec> newTr_3;

  for (INT m = 1 ; m < mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        newTr_3(m, n_alpha, n_gamma) = -fac_yr23 * (sqrt(m + n_alpha) * basis.talmanr(m - 1, n_alpha, m - 1, n_gamma) + sqrt(n_alpha + 1) * basis.talmanr(m - 1, n_alpha + 1, m - 1, n_gamma));
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Yr(0), Yr(1), Yr(2) and Yr(3).
  for (INT m = 1 ; m < mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        arma::vec &yr_2 = Yr(2, m, n_alpha, n_gamma);
        arma::vec &yr_3 = Yr(3, m, n_alpha, n_gamma);
        yr_2 = arma::zeros(Nmax);
        yr_3 = arma::zeros(Nmax);
        arma::vec &newtr_2 = newTr_2(m, n_alpha, n_gamma);
        arma::vec &newtr_3 = newTr_3(m, n_alpha, n_gamma);

        for (INT n_a = 0 ; n_a < Nmax ; n_a++)
        {
          yr_2(n_a) = arma::accu(newtr_2 % basis.moshinskyr(0, n_a));
          yr_3(n_a) = arma::accu(newtr_3 % basis.moshinskyr(0, n_a));
        }//n_a
      }//n_gamma
    }//n_alpha
  }//m

  for (INT m = 1 ; m < mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        arma::vec &yr_1 = Yr(1, m, n_alpha, n_gamma);
        yr_1 = arma::zeros(Nmax);
        arma::vec &newtr_1 = newTr_1(m, n_alpha, n_gamma);

        for (INT n_a = 0 ; n_a < Nmax ; n_a++)
        {
          yr_1(n_a) = arma::accu(newtr_1 % basis.moshinskyr(0, n_a));
        }//n_a
      }//n_gamma
    }//n_alpha
  }//m

  for (INT m = 0 ; m < mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        Yr(0, m, n_alpha, n_gamma) = basis.talmanr(m, n_alpha, m, n_gamma);
      }//n_gamma
    }//n_alpha
  }//m

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldSpinOrbit::calcLegacy(void)
{
  DBG_ENTER;

  if (!fieldLegacy(NEUTRON).empty() && !fieldLegacy(PROTON).empty()) DBG_LEAVE;

  double wso = parameters["wso"];

  warningStr = "";

  calculatingLength = -1.0;
  double startTime = Tools::clock();

  // Please replace me by the new method !
  Discrete discrete2(&basis);
  discrete2.calcDensit(rho(NEUTRON), rho(PROTON ), nGLA, nGHE);
  arma::mat sodens1n = discrete2.BSODensit1(0);
  arma::mat sodens1p = discrete2.BSODensit1(1);

#ifdef NEW_DENSIT
  // new method, equivalent and FASTER
  Mesh specialMesh  = Mesh::gaussLaguerreHermite(nGLA, nGHE);
  double facx = basis.b_r * sqrt(PI);
  double facz = sqrt(basis.b_z);
  double fact = 2.0 / (basis.b_r * basis.b_r * basis.b_z);
  specialMesh.ax.p  = basis.b_r * arma::sqrt(specialMesh.ax.p);
  specialMesh.az.p *= basis.b_z;
  Discrete discrete3(&basis, specialMesh);

  discrete3.isNegative = false;
  arma::mat densn  = discrete3.getLocalXZ(rho(NEUTRON), true) * facx * facx * facz * facz * fact;
  arma::mat densp  = discrete3.getLocalXZ(rho(PROTON ), true) * facx * facx * facz * facz * fact;

  if (discrete3.isNegative) warningStr = "ND";

  // black magic to avoid using Discrete::calcDensit() (it is too slow)
  arma::mat mask = rho(NEUTRON) * 0.0;

  UINT nbBlocks = basis.HOqn.calcBlocks({0, 4});

  for (UINT msb = 0; msb < nbBlocks; msb++)
  {
    Qnumbers &qnb = basis.HOqn.blocks[msb];

    INT mb = qnb(0, 0);
    INT sb = qnb(4, 0);

    if (mb == 0) continue;

    for (UINT msd = 0; msd < nbBlocks; msd++)
    {
      Qnumbers &qnd = basis.HOqn.blocks[msd];

      INT md = qnd(0, 0);
      INT sd = qnd(4, 0);

      if (md == 0) continue;

      if (mb != md) continue;

      if (sb != sd) continue;

      double fac = mb * (sb == 0 ? 1.0 : -1.0) * 0.5 * facx * facx * facz * facz * fact;
      mask.submat(basis.HOqn.blocks[msb].filter, basis.HOqn.blocks[msd].filter) += fac;
    }
  }

  Multi<arma::mat> specialrho;
  specialrho(NEUTRON) = mask % rho(NEUTRON);
  specialrho(PROTON ) = mask % rho(PROTON );

  arma::mat sodensn = discrete3.getLocalXZ(specialrho(NEUTRON), false);
  arma::mat sodensp = discrete3.getLocalXZ(specialrho(PROTON ), false);
  // end of black magic
#else
  // old method
  arma::mat densn    = discrete2.BDensit(0);
  arma::mat densp    = discrete2.BDensit(1);
  arma::mat sodensn  = discrete2.BSODensit(0);
  arma::mat sodensp  = discrete2.BSODensit(1);
#endif


  Discrete discrete(&basis, Mesh::gaussLaguerreHermite(nGLA, nGHE));
  INT Mmax = basis.mMax;
  arma::vec eta = discrete.mesh.ax.p; // radial coordinate
  arma::vec zeta = discrete.mesh.az.p; // axial coordinate
  arma::mat zetashift = arma::zeros(discrete.mesh.az.nb, basis.dMax);

  for (INT dd = 0; dd < basis.dMax; dd++)
  {
    zetashift.col(dd) = zeta - (0.5 - dd) * basis.d_0 / basis.b_z;
  }

  for (INT iiso = 0; iiso < 2; iiso ++) fieldLegacy(iiso) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  arma::mat IntWeight = discrete.mesh.ax.we * discrete.mesh.az.we.t();
  double facint2 = 1. / (basis.b_r * basis.b_z);
  double facSO = 6.0 * 1.04719;
  //------------- Calculation of the sa=sc contribution ------------------------
  INT mc;
  INT Nmax;
  arma::vec Ra;
  arma::vec L0a;
  arma::vec Rc;
  arma::vec L0c;
  arma::vec Rac;
  arma::vec dRaRc;
  arma::vec RadRc;
  arma::vec dRac;
  INT NZmaxa;
  INT NZmaxc;
  arma::vec zda;
  arma::vec zdc;
  UINT ia;
  arma::vec Za;
  UINT ic;
  arma::vec Zc;
  arma::vec Zac;
  arma::mat dRZ;
  arma::mat RZ;
  arma::mat dSO_nup;
  arma::mat dSO_ndown;
  arma::mat dSO_pup;
  arma::mat dSO_pdown;
  arma::mat integ_nup;
  arma::mat integ_ndown;
  arma::mat integ_pup;
  arma::mat integ_pdown;
  double val00ndown;
  double val00nup;
  double val00pdown;
  double val00pup;

  for (INT ma = 0; ma < Mmax; ma++ )
  {
    mc = ma;
    Nmax = basis.nMax(ma);

    for (INT na = 0; na < Nmax; na++ )
    {
      Ra = basis.rPartNorm(eta, ma, na);
      L0a = basis.rPartL0(eta, ma, na);

      for (INT nc = 0; nc < Nmax; nc++ )
      {
        Rc = basis.rPartNorm(eta, mc, nc);
        L0c = basis.rPartL0(eta, mc, nc);
        Rac = Ra % Rc;
        dRaRc = L0a % Rc / arma::sqrt(eta);
        RadRc = Ra % L0c / arma::sqrt(eta);
        dRac = dRaRc + RadRc;
        NZmaxa = basis.n_zMax(ma, na);
        NZmaxc = basis.n_zMax(mc, nc);

        for (INT da = 0; da < basis.dMax; da++ )
        {
          zda = zetashift.col(da);

          for (INT dc = 0; dc < basis.dMax; dc++ )
          {
            zdc = zetashift.col(dc);

            for (INT n_za = 0; n_za < NZmaxa; n_za++ )
            {
              ia = basis.HOqn.find({ma, na, n_za, da, 0}); // spin up
              Za = basis.zPartNorm(zda, n_za);
              ic = basis.HOqn.find({mc, nc, 0, dc, 0}); // spin up

              for (INT n_zc = 0; n_zc < NZmaxc; n_zc ++ )
              {
                Zc = basis.zPartNorm(zdc, n_zc);
                Zac = Za % Zc;
                dRZ = dRac * Zac.t();
                RZ = Rac * Zac.t();
                val00nup   = arma::accu((dRZ % (sodensn - ma / 2. * densn) / pow(basis.b_r, 2) - RZ % sodens1n) % IntWeight);
                val00ndown = arma::accu((dRZ % (sodensn + ma / 2. * densn) / pow(basis.b_r, 2) - RZ % sodens1n) % IntWeight);
                val00pup   = arma::accu((dRZ % (sodensp - ma / 2. * densp) / pow(basis.b_r, 2) - RZ % sodens1p) % IntWeight);
                val00pdown = arma::accu((dRZ % (sodensp + ma / 2. * densp) / pow(basis.b_r, 2) - RZ % sodens1p) % IntWeight);
                fieldLegacy(NEUTRON)(ia, ic) = 2. * val00nup + val00pup;
                fieldLegacy(PROTON)(ia, ic)  = 2. * val00pup + val00nup;

                if (ma > 0)
                {
                  fieldLegacy(NEUTRON)(ia + 1, ic + 1) = 2. * val00ndown + val00pdown;
                  fieldLegacy(PROTON)(ia + 1, ic + 1)  = 2. * val00pdown + val00ndown;
                }

                ic += basis.dMax * (ma ? 2 : 1);
              } // n_zc
            } // n_za
          } // dc
        } // da
      } // nc
    } // na
  } // ma

  fieldLegacy(NEUTRON) *= wso / facSO;
  fieldLegacy(PROTON ) *= wso / facSO;
  //-------------- Calculation of SODensit compo -+ ----------------------------
  INT Nmaxa;
  INT Nmaxc;
  arma::vec Lavecma;
  arma::vec Lavecmc;
  arma::vec A_r1;
  arma::vec A_r2;
  arma::vec S_r1;
  arma::vec S_r2;
  arma::vec dZa;
  arma::vec dZc;
  arma::mat Aac;
  arma::mat Sac;
  arma::mat Intmpn;
  arma::mat Intmpp;
  double valmpn;
  double valmpp;

  for (INT ma = 1; ma < Mmax; ma++ )
  {
    mc = ma - 1;
    Nmaxa = basis.nMax(ma);
    Nmaxc = basis.nMax(mc);

    for (INT na = 0; na < Nmaxa; na++ )
    {
      Ra = basis.rPartNorm(eta, ma, na);
      L0a = basis.rPartL0(eta, ma, na);
      Lavecma = basis.rPartLavecm(eta, ma, na);

      for (INT nc = 0; nc < Nmaxc; nc++ )
      {
        Rc = basis.rPartNorm(eta, mc, nc);
        L0c = basis.rPartL0(eta, mc, nc);
        Lavecmc = basis.rPartLavecm(eta, mc, nc);
        A_r1 = L0a % Rc;
        A_r2 = Ra % L0c;
        S_r1 = Lavecma % Rc;
        S_r2 = Ra % Lavecmc;
        NZmaxa = basis.n_zMax(ma, na);
        NZmaxc = basis.n_zMax(mc, nc);

        for (INT da = 0; da < basis.dMax; da++ )
        {
          zda = zetashift.col(da);

          for (INT dc = 0; dc < basis.dMax; dc++ )
          {
            zdc = zetashift.col(dc);

            for (INT n_za = 0; n_za < NZmaxa; n_za++ )
            {
              ia = basis.HOqn.find({ma, na, n_za, da, 1}); // spin down
              Za = basis.zPartNorm(zda, n_za);
              dZa = basis.zPartNormd(zda, n_za);
              ic = basis.HOqn.find({mc, nc, 0, dc, 0}); // spin up

              for (INT n_zc = 0; n_zc < NZmaxc; n_zc ++ )
              {
                Zc = basis.zPartNorm(zdc, n_zc);
                dZc = basis.zPartNormd(zdc, n_zc);
                Aac = A_r1 * (Za % dZc).t() - A_r2 * (dZa % Zc).t();
                Sac = S_r1 * (Za % dZc).t() + S_r2 * (dZa % Zc).t();
                Intmpn = (Sac + Aac) % densn;
                Intmpp = (Sac + Aac) % densp;
                valmpn = arma::accu(Intmpn % IntWeight);
                valmpp = arma::accu(Intmpp % IntWeight);
                fieldLegacy(NEUTRON)(ia, ic) = facint2 * wso / 2. * (2.*valmpn + valmpp) / facSO;
                fieldLegacy(PROTON )(ia, ic) = facint2 * wso / 2. * (2.*valmpp + valmpn) / facSO;
                ic += basis.dMax * (mc ? 2 : 1);
              } // n_zc
            } // n_za
          } // dc
        } // da
      } // nc
    } // na
  } // ma

  //---------------Calculation of SODensit compo +- ----------------------------
  arma::mat Intpmn;
  arma::mat Intpmp;
  double valpmn;
  double valpmp;

  for (INT ma = 0; ma < Mmax - 1; ma++ )
  {
    mc = ma + 1;
    Nmaxa = basis.nMax(ma);
    Nmaxc = basis.nMax(mc);

    for (INT na = 0; na < Nmaxa; na++ )
    {
      Ra = basis.rPartNorm(eta, ma, na);
      L0a = basis.rPartL0(eta, ma, na);
      Lavecma = basis.rPartLavecm(eta, ma, na);

      for (INT nc = 0; nc < Nmaxc; nc++ )
      {
        Rc = basis.rPartNorm(eta, mc, nc);
        L0c = basis.rPartL0(eta, mc, nc);
        Lavecmc = basis.rPartLavecm(eta, mc, nc);
        A_r1 = L0a % Rc;
        A_r2 = Ra % L0c;
        S_r1 = Lavecma % Rc;
        S_r2 = Ra % Lavecmc;
        NZmaxa = basis.n_zMax(ma, na);
        NZmaxc = basis.n_zMax(mc, nc);

        for (INT da = 0; da < basis.dMax; da++ )
        {
          zda = zetashift.col(da);

          for (INT dc = 0; dc < basis.dMax; dc++ )
          {
            zdc = zetashift.col(dc);

            for (INT n_za = 0; n_za < NZmaxa; n_za++ )
            {
              ia = basis.HOqn.find({ma, na, n_za, da, 0}); // spin up
              Za = basis.zPartNorm(zda, n_za);
              dZa = basis.zPartNormd(zda, n_za);

              for (INT n_zc = 0; n_zc < NZmaxc; n_zc ++ )
              {
                ic = basis.HOqn.find({mc, nc, n_zc, dc, 1}); // spin down
                Zc = basis.zPartNorm(zdc, n_zc);
                dZc = basis.zPartNormd(zdc, n_zc);
                Aac = A_r1 * (Za % dZc).t() - A_r2 * (dZa % Zc).t();
                Sac = S_r1 * (Za % dZc).t() + S_r2 * (dZa % Zc).t();
                Intpmn = (Sac - Aac) % densn;
                Intpmp = (Sac - Aac) % densp;
                valpmn = arma::accu(Intpmn % IntWeight);
                valpmp = arma::accu(Intpmp % IntWeight);
                fieldLegacy(NEUTRON)(ia, ic) = facint2 * wso / 2. * (2.*valpmn + valpmp) / facSO;
                fieldLegacy(PROTON )(ia, ic) = facint2 * wso / 2. * (2.*valpmp + valpmn) / facSO;
              } // n_zc
            } // n_za
          } // dc
        } // da
      } // nc
    } // na
  } // ma

  // Store the calculating length
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

