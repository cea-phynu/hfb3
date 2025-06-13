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

#include "field_density_fr.h"
#include "tools.h"
#include "interaction.h"
#include "discrete.h"

/** \file
 *  \brief Methods of the FieldDensityFR class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldDensityFR::FieldDensityFR(Field::Parameters fp, State *_state) :
  Field(fp, _state)
{
  DBG_ENTER;

  name = "density FR";
  shortName = "DensFR";

  //Select the number of nodes.
  if (basis.dMax == 1)
  {
    nGLA = 16;
    nGHE = 42;
  }
  else
  {
    nGLA = 16;
    nGHE = 38;
  }

  //Builds Discrete instances.
  Discrete _discrete (&_state->basis);
  discrete = _discrete;
  discrete.mesh = Mesh::gaussLaguerreHermite(nGLA, nGHE);
  discrete0 = discrete;
  discrete0.mesh.ax.p = basis.b_r * arma::sqrt(discrete0.mesh.ax.p);
  discrete0.mesh.az.p *= basis.b_z;

  // Multi allocations
  // INT M = basis.mMax;
  // INT dMax = basis.dMax;
  // INT N_p = basis.nGlobalMax;
  // INT N_z = basis.n_zGlobalMax;
  // INT DMax = (dMax - 1) * 3 + 1;

  rVals = FMulti<arma::vec>({MRANGE, NRANGE});

  rValsExchange = FMulti<arma::vec>({MRANGE, {0, (INT) nGLA}});
  rValsExchange2 = FMulti<arma::vec>({MRANGE, NRANGE});

  IzDirect = FMulti<arma::cube>({{0, (INT) nGHE}});
  IzExchange = FMulti<arma::mat>({TRANGE, NZRANGE, {0, (INT) nGHE}});

  RmatDirect        = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  RmatExchange_0    = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});
  Kmat_0            = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});
  RmatExchange_1_p1 = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});
  RmatExchange_1_p2 = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});
  Kmat_1_p1         = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});
  Kmat_1_p2         = FMulti<arma::mat >({TRANGE, MRANGE, NRANGE, NRANGE});

  // Multi<arma::cube> Rzze_0; //Rzze_0(iso, d_gamma,n_zgamma,mp)(n_beta,n_delta,index_z).
  Rzze_0    = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});
  Kzze_0    = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});
  Rzze_1_p1 = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});
  Kzze_1_p1 = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});
  Rzze_1_p2 = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});
  Kzze_1_p2 = FMulti<arma::cube>({TRANGE, DRANGE, NZRANGE, MRANGE});

  // Multi<arma::mat> Rfinal_0; //Rfinal_0(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  Rfinal_0    = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Kfinal_0    = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Rfinal_3    = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Kfinal_3    = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Rfinal_1_p1 = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Kfinal_1_p1 = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Rfinal_1_p2 = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});
  Kfinal_1_p2 = FMulti<arma::mat>({TRANGE, DRANGE, MRANGE, NRANGE, NZRANGE});

  //Rz(iso,n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)
  Rz      = FMulti<arma::vec>({TRANGE, NZRANGE, NZRANGE, DDRANGE});
  Rz2     = FMulti<arma::vec>({TRANGE, NZRANGE, NZRANGE, DDRANGE});

  // Rz_0(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma)
  Rz_0    = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Kz_0    = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Rz_3    = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Kz_3    = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Rz_1_p1 = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Kz_1_p1 = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Rz_1_p2 = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});
  Kz_1_p2 = FMulti<arma::vec>({TRANGE, MRANGE, NRANGE, NZRANGE, NZRANGE, DDRANGE});

  //Gamma(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  Gamma    = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  GammaE_0 = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  Delta_0  = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  GammaE_3 = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  Delta_3  = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  GammaE_1 = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});
  Delta_1  = FMulti<arma::cube>({TRANGE, MRANGE, NRANGE, NRANGE});

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcField(void)
{
  DBG_ENTER;

  if (!field(NEUTRON, DIRECT   ).empty()) DBG_LEAVE;

  //Dependencies.
  calcIz();
  calcIr();
  calcWaveZ();
  calcWaveMatZ();
  calcRVals();
  calcPreR();
  calcPreZ();

  //Useful tools.
  double w3 = parameters["w"];
  double b3 = parameters["b"];
  double h3 = parameters["h"];
  double m3 = parameters["m"];
  double a  = parameters["a"];

  //Usefull quantities
  // INT M = basis.mMax;
  INT dMax = basis.dMax;
  // INT N_p = basis.nGlobalMax;
  // INT N_z = basis.n_zGlobalMax;
  INT DMax = (dMax - 1) * 3 + 1;

  calculatingLength = -1.0;
  double startTime = Tools::clock();

  ///////////////////////////////////////////////////////////Builds Rho Multi with explicit quantum numbers.
  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});
  Multi<arma::mat> Rho;  //Rho(Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)
  Multi<arma::mat> Kappa;  //Kappa(Isospin,m_0,m_1,n_0,n_1,d_0+2*d_1,s_0+2*s_1)(n_z0,n_z1)

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
        Rho(PROTON , m, mp, n, np, d + 2 * dp, s + 2 * sp) = rho(PROTON ).submat(bras.filter, kets.filter);
        Kappa(NEUTRON, m, mp, n, np, d + 2 * dp, s + 2 * sp) = kappa(NEUTRON).submat(bras.filter, kets.filter);
        Kappa(PROTON , m, mp, n, np, d + 2 * dp, s + 2 * sp) = kappa(PROTON ).submat(bras.filter, kets.filter);
      }//if
    }//j
  }//i

  ////////////////////////////////////////////////////////////////////////////////////////////////////////One builds RmatDirect.
  // Multi<arma::cube> RmatDirect;
  // Multi<arma::mat> RmatExchange_0; //RmatExchange_0(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  // Multi<arma::mat> Kmat_0; //KmatExchange_0(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  // Multi<arma::mat> RmatExchange_1_p1; //RmatExchange_1_p1(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  // Multi<arma::mat> RmatExchange_1_p2; //RmatExchange_1_p2(iso,m_delta,n_delta)(VECTORISE(n_zdelta,d_delta),VECTORISE(n_zbeta,d_beta)).
  // Multi<arma::mat> Kmat_1_p1; //KmatExchange_1_p1(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  // Multi<arma::mat> Kmat_1_p2; //KmatExchange_1_p2(iso,m_delta,n_delta)(VECTORISE(n_zdelta,d_delta),VECTORISE(n_zbeta,d_beta)).


  ///////////////////////////////////RmatDirect.
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
        arma::cube &rm = RmatDirect(iso, 0, n_delta, n_beta);
        rm = arma::zeros(basis.n_zMax(0,n_delta),basis.n_zMax(0,n_beta),DMax);
        for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
          {
            rm.slice(d_delta + 2 * d_beta) = (2*w3 + b3) * ( Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(non_iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0))
            - (2*h3 + m3) * Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
    for (INT mp = 1 ; mp < basis.mMax; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::cube &rm = RmatDirect(iso, mp, n_delta, n_beta);
          rm = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);
          for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
          {
            for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
            {
              rm.slice(d_delta + 2 * d_beta) = (2*w3 + b3) * ( Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0)
              + Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3) + Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3))
              - (2*h3 + m3) * (Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3));
            }//d_delta
          }//d_beta
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  ///////////////////////////////////RmatExchange_0 and Kmat_0.
  for (INT iso: {NEUTRON, PROTON})
  {
    //Defines isospin and anti-isospin.
    INT non_iso = 0;
    if (iso == 0)
    {
      non_iso = 1;
    }
    for (INT n_beta = 0 ; n_beta < basis.nMax(0) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
      {
        arma::mat &rme_0 = RmatExchange_0(iso, 0, n_delta, n_beta);
        arma::mat &km_0 = Kmat_0(iso, 0, n_delta, n_beta);
        rme_0 = arma::zeros(basis.dMax * basis.n_zMax(0,n_beta),basis.dMax * basis.n_zMax(0,n_delta));
        km_0 = arma::zeros(basis.dMax * basis.n_zMax(0,n_beta),basis.dMax * basis.n_zMax(0,n_delta));
        INT count = 0; 
        for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
        {
          for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(0,n_beta) ; n_zbeta++)
          {
            arma::mat r = arma::zeros(basis.n_zMax(0,n_delta),basis.dMax);         
            arma::mat k = arma::zeros(basis.n_zMax(0,n_delta),basis.dMax);         
            for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
            {
              r.col(d_delta) = - (w3 - h3 + 2*(b3 - m3)) * Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta)
              + (h3 + 2*m3) * Rho(non_iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
              k.col(d_delta) = (w3 - h3 - b3 + m3) * Kappa(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
            }//d_delta
            rme_0.row(count) = arma::vectorise(r).t();
            km_0.row(count)= arma::vectorise(k).t();
            count += 1;
          }//n_zbeta
        }//d_beta
      }//n_delta
    }//n_beta
    for (INT mp = 1 ; mp < basis.mMax; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &rme_0_plus = RmatExchange_0(iso, mp, n_delta, n_beta);
          arma::mat &rme_0_minus = RmatExchange_0(iso, -mp, n_delta, n_beta);
          arma::mat &km_0_plus = Kmat_0(iso, mp, n_delta, n_beta);
          arma::mat &km_0_minus = Kmat_0(iso, -mp, n_delta, n_beta);
          rme_0_plus = arma::zeros(basis.dMax*basis.n_zMax(mp,n_beta),basis.dMax*basis.n_zMax(mp,n_delta));            
          rme_0_minus = arma::zeros(basis.dMax*basis.n_zMax(mp,n_beta),basis.dMax*basis.n_zMax(mp,n_delta));            
          km_0_plus = arma::zeros(basis.dMax*basis.n_zMax(mp,n_beta),basis.dMax*basis.n_zMax(mp,n_delta)); 
          km_0_minus = arma::zeros(basis.dMax*basis.n_zMax(mp,n_beta),basis.dMax*basis.n_zMax(mp,n_delta)); 
          INT count = 0;
          for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(mp,n_beta) ; n_zbeta++)
            {
              arma::mat r_plus = arma::zeros(basis.n_zMax(mp,n_delta),basis.dMax);
              arma::mat r_minus = arma::zeros(basis.n_zMax(mp,n_delta),basis.dMax);
              arma::mat k_plus = arma::zeros(basis.n_zMax(mp,n_delta),basis.dMax);
              arma::mat k_minus = arma::zeros(basis.n_zMax(mp,n_delta),basis.dMax);
              for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
              {
                r_plus.col(d_delta) = - (w3 - h3 + b3 - m3) * Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta) + (h3 + m3) * Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta)
                - (b3 - m3) * Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta) + m3 * Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta);
                r_minus.col(d_delta) = - (w3 - h3 + b3 - m3) * Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta) + (h3 + m3) * Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta)
                - (b3 - m3) * Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta) + m3 * Rho(non_iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
                k_plus.col(d_delta) = (w3 - h3) * Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta) - (b3 - m3) * Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta);
                k_minus.col(d_delta) = (w3 - h3) * Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta) - (b3 - m3) * Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
              }//d_delta
              rme_0_plus.row(count) = arma::vectorise(r_plus).t();
              rme_0_minus.row(count) = arma::vectorise(r_minus).t();
              km_0_plus.row(count)= arma::vectorise(k_plus).t();
              km_0_minus.row(count)= arma::vectorise(k_minus).t();
              count += 1;
            }//n_zbeta
          }//d_beta
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  ///////////////////////////////////RmatExchange_1 and Kmat_1.
  for (INT iso: {NEUTRON, PROTON})
  {
    //Defines isospin and anti-isospin.
    INT non_iso = 0;
    if (iso == 0)
    {
      non_iso = 1;
    }
    for (INT mp = -basis.mMax + 2 ; mp < 1 ; mp++) //mp <= 0.
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
        {
          arma::mat &rme_1 = RmatExchange_1_p1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1_p1(iso, mp, n_delta, n_beta);
          rme_1 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          km_1 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          INT count = 0;
          for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(abs(mp-1),n_beta) ; n_zbeta++)
            {
              arma::mat r = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              arma::mat k = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
              {
                r.col(d_delta) =  (w3 - h3) * Rho(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).col(n_zbeta) - h3 * Rho(non_iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).col(n_zbeta);
                k.col(d_delta) = - (w3 - h3 + b3 - m3) * Kappa(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).col(n_zbeta);
              }//d_delta
              rme_1.row(count) = arma::vectorise(r).t();
              km_1.row(count) = arma::vectorise(k).t();
              count += 1;
            }//n_zbeta
          }//d_beta
          arma::mat &rme_2 = RmatExchange_1_p2(iso, mp, n_delta, n_beta);
          arma::mat &km_2 = Kmat_1_p2(iso, mp, n_delta, n_beta);
          rme_2 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp), n_delta),basis.dMax*basis.n_zMax(abs(mp-1), n_beta));
          km_2 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp), n_delta),basis.dMax*basis.n_zMax(abs(mp-1), n_beta));
          count = 0;
          for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
          {
            for (INT n_zdelta = 0 ; n_zdelta < basis.n_zMax(abs(mp),n_delta) ; n_zdelta++)
            {
              arma::mat r = arma::zeros(basis.n_zMax(abs(mp-1),n_beta),basis.dMax);
              arma::mat k = arma::zeros(basis.n_zMax(abs(mp-1),n_beta),basis.dMax);
              for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
              {
                r.col(d_beta) =  (w3 - h3) * Rho(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).row(n_zdelta).t() - h3 * Rho(non_iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).row(n_zdelta).t();
                k.col(d_beta) = - (w3 - h3 + b3 - m3) * Kappa(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 2).row(n_zdelta).t();
              }//d_beta
              rme_2.row(count) = arma::vectorise(r).t();
              km_2.row(count) = arma::vectorise(k).t();
              count += 1;
            }//n_zdelta
          }//d_delta
        }//n_delta
      }//n_beta
    }//mp <= 0
    for (INT mp = 1 ; mp < basis.mMax ; mp++) //mp > 0.
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp - 1) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &rme_1 = RmatExchange_1_p1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1_p1(iso, mp, n_delta, n_beta);
          rme_1 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          km_1 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          INT count = 0;
          for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(mp-1,n_beta) ; n_zbeta++)
            {
              arma::mat r = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              arma::mat k = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
              {
                r.col(d_delta) =  -(w3 - h3) * Rho(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).col(n_zbeta) + h3 * Rho(non_iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).col(n_zbeta);
                k.col(d_delta) = (w3 - h3 + b3 - m3) * Kappa(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).col(n_zbeta);
              }//d_delta
              rme_1.row(count) = arma::vectorise(r).t();
              km_1.row(count) = arma::vectorise(k).t();
              count += 1;
            }//n_zbeta
          }//d_beta
          arma::mat &rme_2 = RmatExchange_1_p2(iso, mp, n_delta, n_beta);
          arma::mat &km_2 = Kmat_1_p2(iso, mp, n_delta, n_beta);
          rme_2 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp), n_delta),basis.dMax*basis.n_zMax(abs(mp-1), n_beta));
          km_2 = arma::zeros(basis.dMax*basis.n_zMax(abs(mp), n_delta),basis.dMax*basis.n_zMax(abs(mp-1), n_beta));
          count = 0;
          for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
          {
            for (INT n_zdelta = 0 ; n_zdelta < basis.n_zMax(mp,n_delta) ; n_zdelta++)
            {
              arma::mat r = arma::zeros(basis.n_zMax(abs(mp-1),n_beta),basis.dMax);
              arma::mat k = arma::zeros(basis.n_zMax(abs(mp-1),n_beta),basis.dMax);
              for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
              {
                r.col(d_beta) =  -(w3 - h3) * Rho(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).row(n_zdelta).t() + h3 * Rho(non_iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).row(n_zdelta).t();
                k.col(d_beta) = (w3 - h3 + b3 - m3) * Kappa(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).row(n_zdelta).t();
              }//d_beta
              rme_2.row(count) = arma::vectorise(r).t();
              km_2.row(count) = arma::vectorise(k).t();
              count += 1;
            }//n_zdelta
          }//d_delta
        }//n_delta
      }//n_beta
    }//mp > 0
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////// Builds Rphis and Ris.
  // Multi<arma::mat> RphiDirect; //RphiDirect(iso)(r,z).
  // Multi<arma::mat> RiDirect;   //RiDirect(iso)(r,z).

  FMulti<arma::mat> RphiDirect({TRANGE});
  FMulti<arma::mat> RiDirect({TRANGE});

  arma::mat localRho = arma::pow(discrete0.getLocalXZ(rho(NEUTRON), true)  + discrete0.getLocalXZ(rho(PROTON), true), a);

  //////////////////////////////////RphiDirect.
  for(INT iso: {NEUTRON, PROTON})
  {
    RphiDirect(iso) = arma::zeros(nGLA,nGHE);
    for (INT mp = 0 ; mp < basis.mMax ; mp++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        for (INT n_beta = 0 ; n_beta < n_delta ; n_beta++)
        {
          arma::vec zVec = arma::zeros(nGHE);
          arma::vec rVec = rVals(mp, n_delta) % rVals(mp, n_beta);
          arma::cube &rm = RmatDirect(iso, mp, n_delta, n_beta);

          for (INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            zVec(index_z) = arma::accu(waveZ(index_z).subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_beta) - 1, DMax - 1) % rm);
          }//index_z
          RphiDirect(iso) += 2 * rVec * zVec.t();
        }//n_beta
        arma::vec zVec = arma::zeros(nGHE);
        arma::vec rVec = rVals(mp, n_delta) % rVals(mp, n_delta);
        arma::cube &rm = RmatDirect(iso, mp, n_delta, n_delta);
        for (INT index_z = 0 ; index_z < nGHE ; index_z++)
        {
          zVec(index_z) = arma::accu(waveZ(index_z).subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_delta) - 1, DMax - 1) % rm);
        }//index_z
        RphiDirect(iso) += rVec * zVec.t();
      }//n_delta
    }//mp
    RphiDirect(iso) = RphiDirect(iso) % localRho;
  }//iso

  //////////////////////////////////RiDirect.
  // Multi<arma::mat> Rimatz;

  FMulti<arma::mat> Rimatz({TRANGE, {0, (INT) nGHE}, MRANGE});

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = 0 ; mp < basis.mMax ; mp++)
    {
      for (INT index_z = 0 ; index_z < nGHE ; index_z++)
      {
        Rimatz(iso,index_z,mp) = arma::zeros(basis.nMax(mp),basis.nMax(mp));
      }//index_z
    }//mp
  }//iso
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = 0 ; mp < basis.mMax ; mp++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        for (INT n_beta = 0 ; n_beta < n_delta+1 ; n_beta++)
        {
          arma::cube &rm = RmatDirect(iso, mp, n_delta, n_beta);
          for (INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            Rimatz(iso,index_z,mp)(n_delta,n_beta) = arma::accu(rm % IzDirect(index_z).subcube(0, 0, 0, basis.n_zMax(mp, n_delta) - 1, basis.n_zMax(mp, n_beta) - 1, DMax - 1));
            Rimatz(iso,index_z,mp)(n_beta,n_delta) = Rimatz(iso,index_z,mp)(n_delta,n_beta);
          }//index_z
        }//n_beta
      }//n_delta
    }//mp
  }//iso
  for(INT iso: {NEUTRON, PROTON})
  {
    RiDirect(iso) = arma::zeros(nGLA,nGHE);
    for (INT index_r = 0 ; index_r < nGLA ; index_r++)
    {
      for (INT index_z = 0 ; index_z < nGHE ; index_z++)
      {
        double sum_m = 0.0;
        for (INT mp = 0 ; mp < basis.mMax ; mp++)
        {
          sum_m += arma::accu(Rimatz(iso,index_z,mp)%IrDirect(index_r,mp));
        }//mp
        RiDirect(iso)(index_r,index_z) = sum_m;
      }//index_z
    }//index_r
    RiDirect(iso) = RiDirect(iso) % localRho;
  }//iso

  //////////////////////////////////Builds Rze_0.
  // Multi<arma::cube> Rze_0; //Rze_0(iso, mp, n_beta)(n_delta,n_zbeta,index_z)
  // Multi<arma::cube> Kze_0; //Kze_0(iso, mp, n_beta)(n_delta,n_zbeta,index_z)

  FMulti<arma::cube> Rze_0({TRANGE, MRANGE, NRANGE});
  FMulti<arma::cube> Kze_0({TRANGE, MRANGE, NRANGE});

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
      {

        arma::cube &rze_0 = Rze_0(iso, mp, n_beta);
        arma::cube &kze_0 = Kze_0(iso, mp, n_beta);
        rze_0 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp),n_beta),nGHE);
        kze_0 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp),n_beta),nGHE);

        for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
        {
          arma::mat &rme_0 = RmatExchange_0(iso, mp, n_delta, n_beta);
          arma::mat &km_0 = Kmat_0(iso, mp, n_delta, n_beta);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rze_0.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp),n_beta)-1,index_z) = rme_0 * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kze_0.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp),n_beta)-1,index_z) = km_0  * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  //////////////////////////////////Builds Rze_1.
  // Multi<arma::cube> Rze_1_p1; //Rze_1_p1(iso, mp, n_beta)(n_delta,n_zbeta,index_z)
  // Multi<arma::cube> Kze_1_p1; //Kze_1_p1(iso, mp, n_beta)(n_delta,n_zbeta,index_z)
  // Multi<arma::cube> Rze_1_p2; //Rze_1_p1(iso, mp, n_delta)(n_beta,n_zdelta,index_z)
  // Multi<arma::cube> Kze_1_p2; //Kze_1_p1(iso, mp, n_delta)(n_beta,n_zdelta,index_z)

  FMulti<arma::cube> Rze_1_p1({TRANGE, MRANGE, NRANGE});
  FMulti<arma::cube> Kze_1_p1({TRANGE, MRANGE, NRANGE});
  FMulti<arma::cube> Rze_1_p2({TRANGE, MRANGE, NRANGE});
  FMulti<arma::cube> Kze_1_p2({TRANGE, MRANGE, NRANGE});

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
      {
        arma::cube &rze_1 = Rze_1_p1(iso, mp, n_beta);
        arma::cube &kze_1 = Kze_1_p1(iso, mp, n_beta);
        rze_1 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp-1),n_beta),nGHE);
        kze_1 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp-1),n_beta),nGHE);       
        for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
        {
          arma::mat &rme_1 = RmatExchange_1_p1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1_p1(iso, mp, n_delta, n_beta);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rze_1.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp-1),n_beta)-1,index_z) = rme_1 * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kze_1.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp-1),n_beta)-1,index_z) = km_1  * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_delta
      }//n_beta
      for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
      {
        arma::cube &rze_1 = Rze_1_p2(iso, mp, n_delta);
        arma::cube &kze_1 = Kze_1_p2(iso, mp, n_delta);
        rze_1 = arma::zeros(basis.nMax(abs(mp-1)),basis.dMax * basis.n_zMax(abs(mp),n_delta),nGHE);
        kze_1 = arma::zeros(basis.nMax(abs(mp-1)),basis.dMax * basis.n_zMax(abs(mp),n_delta),nGHE);
        for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
        {
          arma::mat &rme_1 = RmatExchange_1_p2(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1_p2(iso, mp, n_delta, n_beta);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rze_1.subcube(n_beta,0,index_z,n_beta,basis.dMax * basis.n_zMax(abs(mp),n_delta)-1,index_z) = rme_1 * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1));
            kze_1.subcube(n_beta,0,index_z,n_beta,basis.dMax * basis.n_zMax(abs(mp),n_delta)-1,index_z) = km_1  * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1));
          }//index_z
        }//n_beta
      }//n_delta
    }//mp
  }//iso

  ////////////////////////////////////////////////////////////////////////////////////////////Builds Rzze and Kzze exchange 0.
  // Multi<arma::cube> Rzze_0; //Rzze_0(iso, d_gamma,n_zgamma,mp)(n_beta,n_delta,index_z).
  // Multi<arma::cube> Kzze_0; //Kzze_0(iso, d_gamma,n_zgamma,mp)(n_beta,n_deltaindex_z).


  for(INT iso: {NEUTRON, PROTON})
  { 
    for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
    { 
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_0 = Rzze_0(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_0 = Kzze_0(iso, d_gamma,n_zgamma,mp);
          rzze_0 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp)),nGHE);
          kzze_0 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp)),nGHE);

          for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
          {
            arma::cube &rze_0 = Rze_0(iso, mp, n_beta);
            arma::cube &kze_0 = Kze_0(iso, mp, n_beta);
            for(INT index_z = 0 ; index_z < nGHE ; index_z++)
            {
              rzze_0.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = rze_0.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_beta)-1,basis.dMax-1));
              kzze_0.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = kze_0.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_beta)-1,basis.dMax-1));
            }//index_z
          }//n_beta
        }//mp
      }//d_gamma
    }//n_zgamma
  }//iso

  ///////////////////////////////////////////////////////////////////////////////////////////Builds Rzze and Kzze exchange 1.

  // Multi<arma::cube> Rzze_1_p1; //Rzze_1_p1(iso, d_gamma,n_zgamma,mp)(n_beta,n_delta,index_z).
  // Multi<arma::cube> Kzze_1_p1; //Kzze_1_p1(iso, d_gamma,n_zgamma,mp)(n_beta,n_delta,index_z).
  // Multi<arma::cube> Rzze_1_p2; //Rzze_1_p2(iso, d_gamma,n_zgamma,mp)(n_delta,n_beta,index_z).
  // Multi<arma::cube> Kzze_1_p2; //Kzze_1_p2(iso, d_gamma,n_zgamma,mp)(n_delta,n_beta,index_z).

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
    { 
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_1_p1 = Rzze_1_p1(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_1_p1 = Kzze_1_p1(iso, d_gamma,n_zgamma,mp);
          arma::cube &rzze_1_p2 = Rzze_1_p2(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_1_p2 = Kzze_1_p2(iso, d_gamma,n_zgamma,mp);

          rzze_1_p1 = arma::zeros(basis.nMax(abs(mp-1)),basis.nMax(abs(mp)),nGHE);
          kzze_1_p1 = arma::zeros(basis.nMax(abs(mp-1)),basis.nMax(abs(mp)),nGHE);
          rzze_1_p2 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp-1)),nGHE);
          kzze_1_p2 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp-1)),nGHE);
          for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
          {
            arma::cube &rze_1 = Rze_1_p1(iso, mp, n_beta);
            arma::cube &kze_1 = Kze_1_p1(iso, mp, n_beta);
            for(INT index_z = 0 ; index_z < nGHE ; index_z++)
            {
              rzze_1_p1.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = rze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1)); //ICI ETRE PLUS PRECIS sur IzE.
              kzze_1_p1.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = kze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1)); //ICI ETRE PLUS PRECIS SUR IzE.
            }//index_z
          }//n_beta
          for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
          {
            arma::cube &rze_1 = Rze_1_p2(iso, mp, n_delta); 
            arma::cube &kze_1 = Kze_1_p2(iso, mp, n_delta); 
            for(INT index_z = 0 ; index_z < nGHE ; index_z++) 
            {
              rzze_1_p2.subcube(n_delta,0,index_z,n_delta,basis.nMax(abs(mp-1))-1,index_z) = rze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1)); //ICI ETRE PLUS PRECIS sur IzE.
              kzze_1_p2.subcube(n_delta,0,index_z,n_delta,basis.nMax(abs(mp-1))-1,index_z) = kze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1)); //ICI ETRE PLUS PRECIS SUR IzE.
            }//index_z
          }//n_delta
        }//mp
      }//d_gamma
    }//n_zgamma
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Multi<arma::cube> Rzzre_0; //Rzzre_0(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_delta),index_r).
  // Multi<arma::cube> Kzzre_0; //Kzzre_0(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_delta),index_r).

  FMulti<arma::cube> Rzzre_0({TRANGE, DRANGE, NZRANGE});
  FMulti<arma::cube> Kzzre_0({TRANGE, DRANGE, NZRANGE});

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
    { 
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        arma::cube &r_0 = Rzzre_0(iso, d_gamma,n_zgamma);
        arma::cube &k_0 = Kzzre_0(iso, d_gamma,n_zgamma);
        r_0 = arma::zeros(nGHE,index_mn_0,nGLA);
        k_0 = arma::zeros(nGHE,index_mn_0,nGLA);
        for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_0 = Rzze_0(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_0 = Kzze_0(iso, d_gamma,n_zgamma,mp);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              r_0.subcube(index_z,index_M_0(mp - (- basis.mMax + 1 )) ,index_r,index_z, index_M_0(mp+1 - (- basis.mMax + 1) )-1 ,index_r) = rzze_0.slice(index_z) * rValsExchange(abs(mp),index_r);
              k_0.subcube(index_z,index_M_0(mp - (- basis.mMax + 1 )) ,index_r,index_z, index_M_0(mp+1 - (- basis.mMax + 1) )-1 ,index_r) = kzze_0.slice(index_z) * rValsExchange(abs(mp),index_r);
            }//index_r
          }//index_z
        }//mp
      }//d_gamma
    }//n_zgamma
  }//iso

  /////////////////////////////////////////////////////////////////////////////////////////////////
  // Multi<arma::cube> Rzzre_1_p1; //Rzzre_1_p1(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_delta),index_r).
  // Multi<arma::cube> Kzzre_1_p1; //Kzzre_1_p1(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_delta),index_r).
  // Multi<arma::cube> Rzzre_1_p2; //Rzzre_1_p1(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_beta),index_r).
  // Multi<arma::cube> Kzzre_1_p2; //Kzzre_1_p1(iso, d_gamma, n_zgamma)(index_z,VECTOR(mp,n_beta),index_r).

  FMulti<arma::cube> Rzzre_1_p1({TRANGE, DRANGE, NZRANGE});
  FMulti<arma::cube> Kzzre_1_p1({TRANGE, DRANGE, NZRANGE});
  FMulti<arma::cube> Rzzre_1_p2({TRANGE, DRANGE, NZRANGE});
  FMulti<arma::cube> Kzzre_1_p2({TRANGE, DRANGE, NZRANGE});

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
    { 
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        arma::cube &r_1_p1 = Rzzre_1_p1(iso, d_gamma,n_zgamma);
        arma::cube &k_1_p1 = Kzzre_1_p1(iso, d_gamma,n_zgamma);
        arma::cube &r_1_p2 = Rzzre_1_p2(iso, d_gamma,n_zgamma);
        arma::cube &k_1_p2 = Kzzre_1_p2(iso, d_gamma,n_zgamma);

        r_1_p1 = arma::zeros(nGHE,index_mn_1_p1,nGLA);
        k_1_p1 = arma::zeros(nGHE,index_mn_1_p1,nGLA);
        r_1_p2 = arma::zeros(nGHE,index_mn_1_p2,nGLA);
        k_1_p2 = arma::zeros(nGHE,index_mn_1_p2,nGLA);
        for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_1_p1 = Rzze_1_p1(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_1_p1 = Kzze_1_p1(iso, d_gamma,n_zgamma,mp);
          arma::cube &rzze_1_p2 = Rzze_1_p2(iso, d_gamma,n_zgamma,mp);
          arma::cube &kzze_1_p2 = Kzze_1_p2(iso, d_gamma,n_zgamma,mp);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              r_1_p1.subcube(index_z,index_M_1_p1(mp - (- basis.mMax + 2 )) ,index_r,index_z, index_M_1_p1(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = rzze_1_p1.slice(index_z) * rValsExchange(abs(mp),index_r);
              k_1_p1.subcube(index_z,index_M_1_p1(mp - (- basis.mMax + 2 )) ,index_r,index_z, index_M_1_p1(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = kzze_1_p1.slice(index_z) * rValsExchange(abs(mp),index_r);

              r_1_p2.subcube(index_z,index_M_1_p2(mp - (- basis.mMax + 2 )) ,index_r,index_z, index_M_1_p2(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = rzze_1_p2.slice(index_z) * rValsExchange(abs(mp-1),index_r);
              k_1_p2.subcube(index_z,index_M_1_p2(mp - (- basis.mMax + 2 )) ,index_r,index_z, index_M_1_p2(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = kzze_1_p2.slice(index_z) * rValsExchange(abs(mp-1),index_r);
            }//index_r
          }//index_z
        }//mp
      }//d_gamma
    }//n_zgamma
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////
  // Multi<arma::mat> Rfinal_0; //Rfinal_0(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Kfinal_0; //Kfinal_0(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Rfinal_3; //Rfinal_3(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Kfinal_3; //Kfinal_3(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_0 = Rfinal_0(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_0 = Kfinal_0(iso,d_gamma, m,n_gamma,n_zgamma);
            rf_0 = arma::zeros(nGLA,nGHE);
            kf_0 = arma::zeros(nGLA,nGHE);
            arma::cube &r_0 = Rzzre_0(iso, d_gamma,n_zgamma);
            arma::cube &k_0 = Kzzre_0(iso, d_gamma,n_zgamma);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              rf_0.row(index_r) = (r_0.slice(index_r) * IrExchange(index_r,m,n_gamma)).t();
              kf_0.row(index_r) = (k_0.slice(index_r) * IrExchange(index_r,m,n_gamma)).t();
            }//index_r
            rf_0 = rf_0 % localRho;
            kf_0 = kf_0 % localRho;
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_3 = Rfinal_3(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_3 = Kfinal_3(iso,d_gamma, m,n_gamma,n_zgamma);
            rf_3 = arma::zeros(nGLA,nGHE);
            kf_3 = arma::zeros(nGLA,nGHE);
            arma::cube &r_0 = Rzzre_0(iso, d_gamma,n_zgamma);
            arma::cube &k_0 = Kzzre_0(iso, d_gamma,n_zgamma);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              rf_3.row(index_r) = (r_0.slice(index_r) * IrExchange_3(index_r,m,n_gamma)).t();
              kf_3.row(index_r) = (k_0.slice(index_r) * IrExchange_3(index_r,m,n_gamma)).t();
            }//index_z
            rf_3 = rf_3 % localRho;
            kf_3 = kf_3 % localRho;
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////Fourth part ?
  // Multi<arma::mat> Rfinal_1_p1; //Rfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Kfinal_1_p1; //Kfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Rfinal_1_p2; //Rfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  // Multi<arma::mat> Kfinal_1_p2; //Kfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r,index_z).
  
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m-1) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m-1,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_1 = Rfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_1 = Kfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma);
            rf_1 = arma::zeros(nGLA,nGHE);
            kf_1 = arma::zeros(nGLA,nGHE);
            arma::cube &r_1_p1 = Rzzre_1_p1(iso, d_gamma,n_zgamma);
            arma::cube &k_1_p1 = Kzzre_1_p1(iso, d_gamma,n_zgamma);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              rf_1.row(index_r) = (r_1_p1.slice(index_r) * IrExchange_1_p1(index_r,m,n_gamma)).t();
              kf_1.row(index_r) = (k_1_p1.slice(index_r) * IrExchange_1_p1(index_r,m,n_gamma)).t();
            }//index_z
            rf_1 = rf_1 % localRho;
            kf_1 = kf_1 % localRho;
          }//d_gamma
        }//n_zgamma
      }//n_gamma
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m,n_alpha) ; n_zalpha++)
        { 
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            arma::mat &rf_1 = Rfinal_1_p2(iso,d_alpha, m,n_alpha,n_zalpha);
            arma::mat &kf_1 = Kfinal_1_p2(iso,d_alpha, m,n_alpha,n_zalpha);
            rf_1 = arma::zeros(nGLA,nGHE);
            kf_1 = arma::zeros(nGLA,nGHE);
            arma::cube &r_1_p2 = Rzzre_1_p2(iso, d_alpha,n_zalpha);
            arma::cube &k_1_p2 = Kzzre_1_p2(iso, d_alpha,n_zalpha);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              rf_1.row(index_r) = (r_1_p2.slice(index_r) * IrExchange_1_p2(index_r,m,n_alpha)).t();
              kf_1.row(index_r) = (k_1_p2.slice(index_r) * IrExchange_1_p2(index_r,m,n_alpha)).t();
            }//index_z
            rf_1 = rf_1 % localRho;
            kf_1 = kf_1 % localRho;
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
  }//iso

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Numerical Integration (z).
  // Multi<arma::vec> Rz;  //Rz(iso,n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)
  // Multi<arma::vec> Rz2; //Rz2(iso,n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)


  for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
  {
    for (INT n_zgamma = 0 ; n_zgamma < n_zalpha+1 ; n_zgamma++)
    {
      for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
      {
        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          Rz(NEUTRON,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = RiDirect(NEUTRON) * preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz(PROTON ,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = RiDirect(PROTON ) * preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz(NEUTRON,n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz(NEUTRON,n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz(PROTON ,n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz(PROTON ,n_zalpha,n_zgamma,d_alpha+2*d_gamma);

          Rz2(NEUTRON,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = RphiDirect(NEUTRON) * newIzDirect(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz2(PROTON ,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = RphiDirect(PROTON ) * newIzDirect(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz2(NEUTRON,n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz2(NEUTRON,n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz2(PROTON ,n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz2(PROTON ,n_zalpha,n_zgamma,d_alpha+2*d_gamma);
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  // Multi<arma::vec> Rz_0;//Rz_0(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Kz_0;//Kz_0(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Rz_3;//Rz_3(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Kz_3;//Kz_3(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)

  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_0 = Rfinal_0(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_0 = Kfinal_0(iso,d_gamma, m,n_gamma,n_zgamma);
            for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m,0) ; n_zalpha++)
            { 
              for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
              {
                Rz_0(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = rf_0 * zValsExchange(n_zalpha,d_alpha);
                Kz_0(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = kf_0 * zValsExchange(n_zalpha,d_alpha);
              }//d_alpha
            }//n_zalpha
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_3 = Rfinal_3(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_3 = Kfinal_3(iso,d_gamma, m,n_gamma,n_zgamma);
            for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m,0) ; n_zalpha++)
            { 
              for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
              {
                Rz_3(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = rf_3 * zValsExchange(n_zalpha,d_alpha);
                Kz_3(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = kf_3 * zValsExchange(n_zalpha,d_alpha);
              }//d_alpha
            }//n_zalpha
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
  }//iso

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Multi<arma::vec> Rz_1_p1;//Rz_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Kz_1_p1;//Kz_1_p1(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Rz_1_p2;//Rz_1_p2(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)
  // Multi<arma::vec> Kz_1_p2;//Kz_1_p2(iso,d_gamma, m,n_gamma,n_zgamma)(index_r)


  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m-1) ; n_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m-1,n_gamma) ; n_zgamma++)
        { 
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            arma::mat &rf_1 = Rfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma);
            arma::mat &kf_1 = Kfinal_1_p1(iso,d_gamma, m,n_gamma,n_zgamma);
            for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m,0) ; n_zalpha++)
            { 
              for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
              {
                Rz_1_p1(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = rf_1 * zValsExchange(n_zalpha,d_alpha);
                Kz_1_p1(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = kf_1 * zValsExchange(n_zalpha,d_alpha);
              }//d_alpha
            }//n_zalpha
          }//d_gamma
        }//n_zgamma
      }//n_gamma
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m,n_alpha) ; n_zalpha++)
        { 
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            arma::mat &rf_1 = Rfinal_1_p2(iso,d_alpha, m,n_alpha,n_zalpha);
            arma::mat &kf_1 = Kfinal_1_p2(iso,d_alpha, m,n_alpha,n_zalpha);
            for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m-1,0) ; n_zgamma++)
            { 
              for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
              {
                Rz_1_p2(iso,m,n_alpha,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = rf_1 * zValsExchange(n_zgamma,d_gamma);
                Kz_1_p2(iso,m,n_alpha,n_zalpha,n_zgamma,d_alpha+2*d_gamma) = kf_1 * zValsExchange(n_zgamma,d_gamma);
              }//d_alpha
            }//n_zalpha
          }//d_gamma
        }//n_zgamma
      }//n_gamma
    }//m
  }//iso

  ///////////////////////////////////////////////////////////////////////////////////////////////////////Numerical Integration (r).
  // Multi<arma::cube> Gamma; //Gamma(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < basis.mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha+1 ; n_gamma++)
        {
          arma::vec &prer = preR(m, n_alpha, n_gamma);
          arma::vec &newird = newIrDirect(m, n_alpha, n_gamma);
          arma::cube &gamma = Gamma(iso, m, n_alpha, n_gamma);
          arma::cube &gamma_sym = Gamma(iso, m, n_gamma, n_alpha);
          gamma = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
          gamma_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
            {
              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(prer % Rz(iso, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                  gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += arma::accu(newird % Rz2(iso, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                }//n_zgamma
              }//n_zalpha
              gamma_sym.slice(d_gamma+2*d_alpha) = gamma.slice(d_alpha+2*d_gamma).t();
            }//d_gamma
          }//d_alpha
        }//n_gamma
      }//n_alpha
    }//m
  }//iso

  // Multi<arma::cube> GammaE_0;//GammaE_0(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Delta_0;//Delta_0(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> GammaE_3;//GammaE_3(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Delta_3;//Delta_3(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < basis.mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
        {
          arma::cube &gammae_0 = GammaE_0(iso, m, n_alpha, n_gamma);
          arma::cube &delta_0 = Delta_0(iso, m, n_alpha, n_gamma);
          gammae_0 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m, n_gamma),DMax);
          delta_0 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m, n_gamma),DMax);
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
            {
              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  gammae_0(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Rz_0(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  delta_0(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Kz_0(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//n_gamma
      }//n_alpha
    }//m
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
        {
          arma::cube &gammae_3 = GammaE_3(iso, m, n_alpha, n_gamma);
          arma::cube &delta_3 = Delta_3(iso, m, n_alpha, n_gamma);
          gammae_3 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m, n_gamma),DMax);
          delta_3 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m, n_gamma),DMax);
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
            {
              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  gammae_3(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Rz_3(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  delta_3(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Kz_3(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//n_gamma
      }//n_alpha
    }//m
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Last Part.
  // Multi<arma::cube> GammaE_1;//GammaE_1(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Delta_1;//Delta_1(iso, m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).


  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 1 ; m < basis.mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m-1) ; n_gamma++)
        {
          arma::cube &gammae_1 = GammaE_1(iso, m, n_alpha, n_gamma);
          arma::cube &delta_1 = Delta_1(iso, m, n_alpha, n_gamma);
          gammae_1 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m-1, n_gamma),DMax);
          delta_1 = arma::zeros(basis.n_zMax(m, n_alpha),basis.n_zMax(m-1, n_gamma),DMax);
          for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
          {
            for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
            {
              for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m-1, n_gamma) ; n_zgamma++)
                {
                  gammae_1(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Rz_1_p1(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  gammae_1(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += arma::accu(rValsExchange2(m-1,n_gamma) % Rz_1_p2(iso,m,n_alpha,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  delta_1(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(rValsExchange2(m,n_alpha) % Kz_1_p1(iso,m,n_gamma,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                  delta_1(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) += arma::accu(rValsExchange2(m-1,n_gamma) % Kz_1_p2(iso,m,n_alpha,n_zalpha,n_zgamma,d_alpha+2*d_gamma));
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//n_gamma
      }//n_alpha
    }//m
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds the interaction..

  field(NEUTRON, DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(NEUTRON, EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(NEUTRON, PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

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
        field(NEUTRON, DIRECT  ).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , DIRECT  ).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GammaE_0(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GammaE_0(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(NEUTRON, PAIRING ).submat(bras.filter, kets.filter) = Delta_0(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , PAIRING ).submat(bras.filter, kets.filter) = Delta_0(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean Exchange part --.
      {
        field(NEUTRON, DIRECT  ).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , DIRECT  ).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GammaE_3(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GammaE_3(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(NEUTRON, PAIRING ).submat(bras.filter, kets.filter) = Delta_3(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , PAIRING ).submat(bras.filter, kets.filter) = Delta_3(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if
    }//j
  }//i

  //Transposition for exchange and kappa ++ and --.
  field(NEUTRON, EXCHANGE) = field(NEUTRON, EXCHANGE) + field(NEUTRON, EXCHANGE).t();
  field(PROTON , EXCHANGE) = field(PROTON , EXCHANGE) + field(PROTON , EXCHANGE).t();
  field(NEUTRON, PAIRING ) = field(NEUTRON, PAIRING ) + field(NEUTRON, PAIRING ).t();
  field(PROTON , PAIRING ) = field(PROTON , PAIRING ) + field(PROTON , PAIRING ).t();

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

      if (m_alpha == m_gamma + 1 && m_alpha > 0 && s_alpha == 1 && s_gamma == 0) //Mean exchange part -+ and Pairing part -+.
      {
        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GammaE_1(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GammaE_1(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(NEUTRON, PAIRING ).submat(bras.filter, kets.filter) = Delta_1(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON , PAIRING ).submat(bras.filter, kets.filter) = Delta_1(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma - 1 && m_alpha < basis.mMax - 1 && s_alpha == 0 && s_gamma == 1) //Mean exchange part +- and Pairing part +-.
      {
        field(NEUTRON, EXCHANGE).submat(bras.filter, kets.filter) = GammaE_1(NEUTRON, m_alpha+1, n_gamma, n_alpha).slice(d_gamma + 2 * d_alpha).t();
        field(PROTON , EXCHANGE).submat(bras.filter, kets.filter) = GammaE_1(PROTON, m_alpha+1, n_gamma, n_alpha).slice(d_gamma + 2 * d_alpha).t();
        field(NEUTRON, PAIRING ).submat(bras.filter, kets.filter) = Delta_1(NEUTRON, m_alpha+1, n_gamma, n_alpha).slice(d_gamma + 2 * d_alpha).t();
        field(PROTON , PAIRING ).submat(bras.filter, kets.filter) = Delta_1(PROTON, m_alpha+1, n_gamma, n_alpha).slice(d_gamma + 2 * d_alpha).t();
      }//if

    }//j
  }//i

  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcIz(void)
{
  DBG_ENTER;

  if (!IzDirect.empty()) DBG_LEAVE;

  //Dependencies.
  basis.calcTalmanz();

  //Useful tools.
  double p3   =  parameters["p"];
  double d2_factor = 1.0/(2.0 * std::pow(p3*sqrt(PI),3));
  double z_factor  = p3 * sqrt(PI)/ (sqrt(basis.b_z * sqrt(PI)));
  double Kz = (p3*p3 + basis.b_z * basis.b_z)/(basis.b_z * basis.b_z);
  INT DMax = (basis.dMax - 1) * 3 + 1;  

  //Builds the "gaussian integral".
  Multi<arma::vec> gaussianPhi; //gaussianPhi(z,d_delta+2*d_beta)(n_za).
  for (INT index_z = 0; index_z < nGHE ; index_z++)
  {
    double value_z = discrete0.mesh.az.p(index_z);
    for (INT d_delta = 0; d_delta < basis.dMax; d_delta++)
    {
      double zd_delta = -(0.5 - d_delta) * basis.d_0;
      for (INT d_beta = 0; d_beta < basis.dMax; d_beta++)
      {
        double zd_beta = -(0.5 - d_beta) * basis.d_0;
        double k_db = (zd_delta + zd_beta)/2.0;
        double exp_factor = exp(-0.5 * std::pow((value_z + k_db)/(basis.b_z * sqrt(Kz)),2));
        arma::vec &gphi = gaussianPhi(index_z,d_delta+2*d_beta);
        gphi = arma::zeros(2 * basis.n_zGlobalMax + 2);
        for (INT n_za = 0; n_za < 2 * basis.n_zGlobalMax + 2; n_za++)
        {
          gphi(n_za) = d2_factor * z_factor * exp_factor * basis.zPartScalar((value_z + k_db)/sqrt(Kz), n_za,-1) / std::pow(sqrt(Kz),n_za + 1);
        }//n_za
      }//d_beta
    }//d_delta
  }//index_z

  for (INT index_z = 0; index_z < nGHE ; index_z++)
  {
    arma::cube &izdirect = IzDirect(index_z);
    izdirect = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax,DMax);
    for (INT d_delta = 0; d_delta < basis.dMax; d_delta++)
    {
      for (INT d_beta = 0; d_beta < basis.dMax; d_beta++)
      {
        arma::vec &gphi = gaussianPhi(index_z,d_delta+2*d_beta);
        for (INT n_zdelta = 0; n_zdelta < basis.n_zGlobalMax; n_zdelta++)
        {
          for (INT n_zbeta = 0; n_zbeta < n_zdelta+1; n_zbeta++)
          {
            izdirect(n_zdelta,n_zbeta,d_delta+2*d_beta) = arma::accu(basis.talmanz(n_zdelta, d_delta, n_zbeta, d_beta) % gphi);
            izdirect(n_zbeta,n_zdelta,d_beta+2*d_delta) = izdirect(n_zdelta,n_zbeta,d_delta+2*d_beta);
          }//n_zbeta
        }//n_zdelta
      }//d_beta
    }//d_delta
  }//index_z

  //Reshape IzDirect.
  for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
  {
    for (INT n_zgamma = 0 ; n_zgamma < n_zalpha+1 ; n_zgamma++) 
    {
      for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
      {
        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          arma::vec &niz = newIzDirect(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          niz = arma::zeros(nGHE);
          for (INT index_z = 0; index_z < nGHE ; index_z++)
          {
            niz(index_z) = IzDirect(index_z)(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          }//index_z
          niz = discrete0.mesh.az.we % niz;
          newIzDirect(n_zgamma,n_zalpha,d_gamma+2*d_alpha) = niz;
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  //Reshape IzExchange.
  for (INT index_z = 0; index_z < nGHE ; index_z++)
  {
    for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
    {
      for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
      {
        arma::mat &ize = IzExchange(d_alpha,n_zalpha,index_z);
        ize = arma::zeros(basis.n_zGlobalMax,basis.dMax);
        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++) 
          {
            ize(n_zgamma,d_gamma) = IzDirect(index_z)(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          }//n_zgamma
        }//d_gamma
      }//d_alpha
    }//n_zalpha
  }//index_z

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcIr(void)
{
  DBG_ENTER;

  if (!IrDirect.empty()) DBG_LEAVE;

  //Dependencies.
  basis.calcTalmanr();

  //Useful tools.
  double p3   =  parameters["p"];
  double Kr = (p3*p3 + basis.b_r * basis.b_r)/(basis.b_r * basis.b_r);
  double r_factor = p3 * p3 * sqrt(PI) / basis.b_r;

  //Builds the "gaussian integral".
  Multi<arma::vec> gaussianPhi; //gaussianPhi(r)(n_a).
  for (INT index_r = 0; index_r < nGLA ; index_r++)
  {
    double value_r = discrete0.mesh.ax.p(index_r);
    double exp_factor = exp(-0.5 * std::pow(value_r/(basis.b_r * sqrt(Kr)),2));
    arma::vec &gphi = gaussianPhi(index_r);
    gphi = arma::zeros(basis.Nmaxr);    
    for(INT n_a = 0 ; n_a < basis.Nmaxr ; n_a++)
    {
      gphi(n_a) = r_factor * exp_factor * basis.rPartScalar(value_r/sqrt(Kr),0,n_a) /std::pow(sqrt(Kr),2*n_a+2);
    }//n_a
  }//index_r

  for(INT index_r = 0; index_r < nGLA ; index_r++)
  {
    arma::vec &gphi = gaussianPhi(index_r);
    for(INT mp = 0 ; mp < basis.mMax ; mp++)
    {
      arma::mat &irdirect = IrDirect(index_r,mp);
      irdirect = arma::zeros(basis.nMax(mp),basis.nMax(mp));
      for(INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        for(INT n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++)
        {
          irdirect(n_delta,n_beta) = arma::accu(basis.talmanr(mp, n_beta, mp, n_delta) % gphi);
          irdirect(n_beta,n_delta) = irdirect(n_delta,n_beta);
        }//n_beta
      }//n_delta
    }//mp
  }//index_r

  //Reshape Ir.
  double fac_r_integral = PI * basis.b_z * std::pow(basis.b_r,2);
  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha+1 ; n_gamma++)
      {
        arma::vec &nir = newIrDirect(m,n_alpha,n_gamma);
        nir = arma::zeros(nGLA);
        for (INT index_r = 0 ; index_r < nGLA ; index_r++)
        {
          nir(index_r) = IrDirect(index_r,m)(n_alpha,n_gamma);
        }//index_r
        nir = fac_r_integral * discrete0.mesh.ax.we % nir;
        newIrDirect(m,n_gamma,n_alpha) = nir;
      }//n_gamma
    }//n_alpha
  }//m

  ////////////////////////////////////////////////////////////////////////////////Exchange part.
  index_mn_0 = 0;
  index_M_0 = arma::zeros<IVEC>(2*basis.mMax);
  index_M_0(0) = 0;
  for(INT m = -basis.mMax + 1; m < basis.mMax ; m++)
  {
    for(INT n = 0; n < basis.nMax(abs(m)) ; n++)
    {
      index_mn_0 += 1;
    }//n
    index_M_0(m - (-basis.mMax + 1) + 1) = index_mn_0;
  }//m

  index_mn_1_p1 = 0;
  index_mn_1_p2 = 0;
  index_M_1_p1 = arma::zeros<IVEC>(2*basis.mMax-1);
  index_M_1_p2 = arma::zeros<IVEC>(2*basis.mMax-1);
  index_M_1_p1(0) = 0;
  index_M_1_p2(0) = 0;
  for(INT m = -basis.mMax + 2; m < basis.mMax ; m++)
  {
    for(INT n = 0; n < basis.nMax(abs(m)) ; n++)
    {
      index_mn_1_p1 += 1;
    }//n
    index_M_1_p2(m - (-basis.mMax + 2) + 1) = index_mn_1_p1;
    for(INT n = 0; n < basis.nMax(abs(m-1)) ; n++)
    {
      index_mn_1_p2 += 1;
    }//n
     index_M_1_p1(m - (-basis.mMax + 2) + 1) = index_mn_1_p2;
  }//m
  //Builds the Exchange "gaussian integral".
  Multi<arma::vec> gaussianPhiExchange; //gaussianPhiExchange(r,ma)(n_a).
  for (INT index_r = 0; index_r < nGLA ; index_r++)
  {
    double value_r = discrete0.mesh.ax.p(index_r);
    double exp_factor = exp(-0.5 * std::pow(value_r/(basis.b_r * sqrt(Kr)),2));
    for (INT ma = -2*basis.mMax + 2; ma < 2 * basis.mMax - 1; ma++)
    {
      arma::vec &gphi = gaussianPhiExchange(index_r,ma);
      gphi = arma::zeros(basis.Nmaxr); 
      for(INT n_a = 0 ; n_a < basis.Nmaxr ; n_a++)
      {
        gphi(n_a) = r_factor * exp_factor * basis.rPartScalar(value_r/sqrt(Kr),ma,n_a) /std::pow(sqrt(Kr),2*n_a+2 + abs(ma));
      }//n_a
    }//ma
  }//index_r
  for(INT index_r = 0; index_r < nGLA ; index_r++)
  {
    for(INT m = 0 ; m < basis.mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        arma::vec &irexchange = IrExchange(index_r,m,n_alpha);
        irexchange = arma::zeros(index_mn_0);
        INT count = 0;
        for(INT mp = -basis.mMax + 1; mp < basis.mMax ; mp++)
        {
          arma::vec &gphi = gaussianPhiExchange(index_r,m-mp);
          for(INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
          {
            irexchange(count) = arma::accu(basis.talmanr(mp, n_beta, m, n_alpha) % gphi);
            count += 1;
          }//n_beta
        }//n_alpha
      }//mp
    }//m
    for(INT m = 1 ; m < basis.mMax ; m++)
    {
      for(INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        arma::vec &irexchange = IrExchange_3(index_r,m,n_alpha);
        irexchange = arma::zeros(index_mn_0);
        INT count = 0;
        for(INT mp = basis.mMax -1; mp > -basis.mMax ; mp+=-1)
        {
          arma::vec &gphi = gaussianPhiExchange(index_r,m-mp);
          for(INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
          {
            irexchange(count) = arma::accu(basis.talmanr(mp, n_beta, m, n_alpha) % gphi);
            count += 1;
          }//n_beta
        }//n_alpha
      }//mp
    }//m
  }//index_r

  /////////////////////////////////////////////////////// -+ Part.
  for(INT index_r = 0; index_r < nGLA ; index_r++)
  {
    for(INT m = 1 ; m < basis.mMax ; m++)
    {
      for(INT n_gamma = 0 ; n_gamma < basis.nMax(m-1) ; n_gamma++)
      {
        arma::vec &irexchange = IrExchange_1_p1(index_r,m,n_gamma);
        irexchange = arma::zeros(index_mn_1_p1);
        INT count = 0;
        for(INT mp = -basis.mMax + 2; mp < basis.mMax ; mp++)
        {
          arma::vec &gphi = gaussianPhiExchange(index_r,m-mp);
          for(INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
          {
            irexchange(count) = arma::accu(basis.talmanr(mp-1, n_beta, m-1, n_gamma) % gphi);
            count += 1;
          }//n_beta
        }//n_gamma
      }//mp
      for(INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        arma::vec &irexchange = IrExchange_1_p2(index_r,m,n_alpha);
        irexchange = arma::zeros(index_mn_1_p2);
        INT count = 0;
        for(INT mp = -basis.mMax + 2; mp < basis.mMax ; mp++)
        {
          arma::vec &gphi = gaussianPhiExchange(index_r,m-mp);
          for(INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
          {
            irexchange(count) = arma::accu(basis.talmanr(mp, n_delta, m, n_alpha) % gphi);
            count += 1;
          }//n_delta
        }//n_gamma
      }//mp
    }//m
  }//index_r

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcWaveZ(void)
{
  DBG_ENTER;

  if (!waveZ.empty()) DBG_LEAVE;

  //Useful tools.
  INT DMax = (basis.dMax - 1) * 3 + 1;

  //z wave function values.
  Multi<arma::vec> zVals;
  for (INT n_z = 0; n_z < basis.n_zGlobalMax; n_z++)
  {
    for (INT d = 0; d < basis.dMax; d++)
    { 
      double dv = -(0.5 - d) * basis.d_0;
      zVals(n_z,d) = arma::zeros(nGHE);
      for(INT index_z = 0 ; index_z < nGHE ; index_z++)
      {
        zVals(n_z,d)(index_z) = basis.zPartScalar(discrete0.mesh.az.p(index_z) + dv, n_z);
      }//index_z
      zValsExchange(n_z,d) = zVals(n_z,d) % discrete0.mesh.az.we;
    }//d
  }//n_z
  
  for (INT index_z = 0 ; index_z < nGHE ; index_z++)
  {
    arma::cube &wavecube = waveZ(index_z);
    wavecube = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax, DMax);

    for (INT d = 0; d < basis.dMax ; d++)
    {
      for (INT dp = 0; dp < basis.dMax ; dp++)
      {
        for (INT n_z = 0; n_z < basis.n_zGlobalMax ; n_z++)
        {
          for (INT n_zp = 0; n_zp < n_z+1 ; n_zp++)
          {
            wavecube(n_z, n_zp, d + 2 * dp) = zVals(n_z, d)(index_z) * zVals(n_zp, dp)(index_z);
            wavecube(n_zp, n_z, dp + 2 * d) = wavecube(n_z, n_zp, d + 2 * dp);
          }//n_zp
        }//n_z
      }//d
    }//dp
  }//index_z

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcWaveMatZ(void)
{
  DBG_ENTER;

  if (!waveMatZ.empty()) DBG_LEAVE;

  for(INT index_z = 0 ; index_z < nGHE ; index_z++)
  {
    arma::mat &zvm = waveMatZ(index_z);
    zvm = arma::zeros(basis.n_zGlobalMax,basis.dMax);
    for (INT n_z = 0; n_z < basis.n_zGlobalMax; n_z++)
    {
      for (INT d = 0; d < basis.dMax; d++)
      {
        double dv = -(0.5 - d) * basis.d_0;
        zvm(n_z,d) = basis.zPartScalar(discrete0.mesh.az.p(index_z) + dv, n_z);
      }//d
    }//n_z
  }//index_z

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcRVals(void)
{
  DBG_ENTER;

  if (!rVals.empty()) DBG_LEAVE;

  double fac_r_integral = PI * basis.b_z * std::pow(basis.b_r,2);

  //r wave function values.
  for (INT m = 0; m < basis.mMax; m++)
  {
    for (INT n = 0; n < basis.nMax(m) ; n++)
    {
      rVals(m, n) = basis.rPart(discrete0.mesh.ax.p, m, n);
    }//n
  }//m

  for (INT m = 0; m < basis.mMax; m++)
  {
    for(INT index_r = 0 ; index_r < nGLA ; index_r++)
    {
      arma::vec &rve = rValsExchange(m, index_r);
      rve = arma::zeros(basis.nMax(m));
      for (INT n = 0; n < basis.nMax(m) ; n++)
      {
        rve(n) = rVals(m, n)(index_r);
      }//n
    }//index_r
  }//m

  for (INT m = 0; m < basis.mMax; m++)
  {
    for (INT n = 0; n < basis.nMax(m) ; n++)
    {
      rValsExchange2(m, n) = fac_r_integral * rVals(m,n) % discrete0.mesh.ax.we;
    }//n
  }//m

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcPreR(void)
{
  DBG_ENTER;

  if (!preR.empty()) DBG_LEAVE;

  double fac_r_integral = PI * basis.b_z * std::pow(basis.b_r,2);
  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha+1 ; n_gamma++)
      {
        preR(m, n_alpha, n_gamma) = fac_r_integral * discrete0.mesh.ax.we % basis.rPart(discrete0.mesh.ax.p,m,n_alpha) % basis.rPart(discrete0.mesh.ax.p,m,n_gamma);
        preR(m, n_gamma, n_alpha) = preR(m, n_alpha, n_gamma);
      }//n_gamma
    }//n_alpha
  }//m

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensityFR::calcPreZ(void)
{
  DBG_ENTER;

  if (!preZ.empty()) DBG_LEAVE;
  
  for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
  {
    for (INT n_zgamma = 0 ; n_zgamma < n_zalpha+1 ; n_zgamma++)
    {
      for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
      {
        double zd_alpha = -(0.5 - d_alpha) * basis.d_0;
        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          double zd_gamma = -(0.5 - d_gamma) * basis.d_0;
          arma::vec &prez = preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          prez = arma::zeros(nGHE);
          for (INT index_z = 0; index_z < nGHE ; index_z++)
          {
            double value_z = discrete0.mesh.az.p(index_z);
            prez(index_z) = basis.zPartScalar(value_z + zd_alpha, n_zalpha) * basis.zPartScalar(value_z + zd_gamma, n_zgamma);
          }//index_z
          prez = prez % discrete0.mesh.az.we;
          preZ(n_zgamma,n_zalpha,d_gamma+2*d_alpha) = prez;
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  DBG_LEAVE;
}

