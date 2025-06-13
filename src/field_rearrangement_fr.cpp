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

#include "field_rearrangement_fr.h"
#include "tools.h"
#include "basis.h"
#include "interaction.h"
#include "discrete.h"

/** \file
 *  \brief Methods of the FieldRearrangementFR class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldRearrangementFR::FieldRearrangementFR(Field::Parameters fp, State *_state) :
  Field(fp, _state)
{
  DBG_ENTER;

  name = "rearrangement FR";
  shortName = "RearFR";

  contributeToEnergy = false;

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

  Discrete _discrete (&_state->basis);
  discrete = _discrete;
  discrete.mesh = Mesh::gaussLaguerreHermite(nGLA, nGHE);

  discrete0 = discrete;
  discrete0.mesh.ax.p = basis.b_r * arma::sqrt(discrete0.mesh.ax.p);
  discrete0.mesh.az.p *= basis.b_z;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangementFR::calcField(void)
{
  DBG_ENTER;

  if (!field(NEUTRON, DIRECT  ).empty()) DBG_LEAVE;

  double startTime = Tools::clock();

  //Dependencies.
  calcIz();
  calcIr();
  calcWaveZ();
  calcWaveMatZ();
  calcIndexZ();
  calcRVals();
  calcPreR();
  calcPreZ();

  //Useful tools.
  double w3 =  parameters["w"];
  double b3 =  parameters["b"];
  double h3 =  parameters["h"];
  double m3 =  parameters["m"];
  double x  =  parameters["a"];
  INT DMax = (basis.dMax - 1) * 3 + 1;

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
  Multi<arma::cube> RmatDirect;

  Multi<arma::mat> Rmat_pp;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Rmat_mm;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Kmat_pp;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Kmat_mm;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> RmatExchange_0; //RmatExchange_0(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Kmat_0; //KmatExchange_0(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).

  Multi<arma::mat> Rmat_mp;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Kmat_mp;//(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> RmatExchange_1; //RmatExchange_1(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).
  Multi<arma::mat> Kmat_1; //KmatExchange_1(iso,m_delta,n_delta)(VECTORISE(n_zbeta,d_beta),VECTORISE(n_zdelta,d_delta)).

  ///////////////////////////////////RmatDirect.
  for (INT iso: {NEUTRON, PROTON})
  {
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
            rm.slice(d_delta + 2 * d_beta) = Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0);
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
              rm.slice(d_delta + 2 * d_beta) = Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0) + Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3);
            }//d_delta
          }//d_beta
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  ///////////////////////////////////Rmat_pp, Kmat_pp, Rmat_mm et Kmat_mm.
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(0) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
      {
        arma::mat &r_pp = Rmat_pp(iso, 0, n_delta, n_beta);
        arma::mat &k_pp = Kmat_pp(iso, 0, n_delta, n_beta);
        r_pp = arma::zeros(basis.dMax * basis.n_zMax(0,n_beta),basis.dMax * basis.n_zMax(0,n_delta));
        k_pp = arma::zeros(basis.dMax * basis.n_zMax(0,n_beta),basis.dMax * basis.n_zMax(0,n_delta));
        INT count = 0;
        for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
        {
          for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(0,n_beta) ; n_zbeta++)
          {
            arma::mat r = arma::zeros(basis.n_zMax(0,n_delta),basis.dMax);
            arma::mat k = arma::zeros(basis.n_zMax(0,n_delta),basis.dMax);
            for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
            {
              r.col(d_delta) = Rho(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
              k.col(d_delta) = Kappa(iso, 0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
            }//d_delta
            r_pp.row(count) = arma::vectorise(r).t();
            k_pp.row(count) = arma::vectorise(k).t();
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
          arma::mat &r_pp = Rmat_pp(iso, mp, n_delta, n_beta);
          arma::mat &k_pp = Kmat_pp(iso, mp, n_delta, n_beta);
          arma::mat &r_mm = Rmat_mm(iso, mp, n_delta, n_beta);
          arma::mat &k_mm = Kmat_mm(iso, mp, n_delta, n_beta);
          r_pp = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.dMax * basis.n_zMax(mp,n_delta));
          k_pp = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.dMax * basis.n_zMax(mp,n_delta));
          r_mm = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.dMax * basis.n_zMax(mp,n_delta));
          k_mm = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.dMax * basis.n_zMax(mp,n_delta));
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
                r_plus.col(d_delta) = Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
                r_minus.col(d_delta) = Rho(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta);
                k_plus.col(d_delta) = Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0).col(n_zbeta);
                k_minus.col(d_delta) = Kappa(iso, mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3).col(n_zbeta);
              }//d_delta
              r_pp.row(count) = arma::vectorise(r_plus).t();
              r_mm.row(count) = arma::vectorise(r_minus).t();
              k_pp.row(count) = arma::vectorise(k_plus).t();
              k_mm.row(count) = arma::vectorise(k_minus).t();
              count += 1;
            }//n_zbeta
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

  ///////////////////////////////////Rmat_mp and Kmat_mp.
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = 1 ; mp < basis.mMax ; mp++) //mp > 0.
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp - 1) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &rm_mp = Rmat_mp(iso, mp, n_delta, n_beta);
          arma::mat &km_mp = Kmat_mp(iso, mp, n_delta, n_beta);
          rm_mp = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          km_mp = arma::zeros(basis.dMax*basis.n_zMax(abs(mp-1), n_beta),basis.dMax*basis.n_zMax(abs(mp), n_delta));
          INT count = 0;
          for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(mp-1,n_beta) ; n_zbeta++)
            {
              arma::mat r = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              arma::mat k = arma::zeros(basis.n_zMax(abs(mp),n_delta),basis.dMax);
              for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
              {
                r.col(d_delta) =  Rho(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).col(n_zbeta);
                k.col(d_delta) =  Kappa(iso, abs(mp), abs(mp - 1), n_delta, n_beta, d_delta+2*d_beta, 1).col(n_zbeta);
              }//d_delta
              rm_mp.row(count) = arma::vectorise(r).t();
              km_mp.row(count) = arma::vectorise(k).t();
              count += 1;
            }//n_zbeta
          }//d_beta
        }//n_delta
      }//n_beta
    }//mp > 0
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
          arma::mat &rme_1 = RmatExchange_1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1(iso, mp, n_delta, n_beta);
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
        }//n_delta
      }//n_beta
    }//mp <= 0
    for (INT mp = 1 ; mp < basis.mMax ; mp++) //mp > 0.
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp - 1) ; n_beta++)
      {
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &rme_1 = RmatExchange_1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1(iso, mp, n_delta, n_beta);
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
        }//n_delta
      }//n_beta
    }//mp > 0
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////// Builds Rphis and Ris.
  Multi<arma::mat> RphiDirect; //RphiDirect(iso)(r,z).
  Multi<arma::mat> RiDirect;   //RiDirect(iso)(r,z).
  arma::mat localRho = arma::pow(discrete0.getLocalXZ(rho(NEUTRON), true)  + discrete0.getLocalXZ(rho(PROTON), true),x-1.0);

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
  }//iso

  //////////////////////////////////RiDirect.
  Multi<arma::mat> Rimatz;
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
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////// Rphi_pp and Rphi_mm.
  Multi<arma::cube> Rphiz_pp; // Rphiz_pp(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  Multi<arma::cube> Kphiz_pp; // Kphiz_pp(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  Multi<arma::cube> Rphiz_mm; // Rphiz_mm(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  Multi<arma::cube> Kphiz_mm; // Kphiz_mm(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = 0 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        arma::cube &rz_pp = Rphiz_pp(iso,mp,n_beta);
        arma::cube &kz_pp = Kphiz_pp(iso,mp,n_beta);
        rz_pp = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.nMax(mp),nGHE);
        kz_pp = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.nMax(mp),nGHE);
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &r_pp = Rmat_pp(iso, mp, n_delta, n_beta);
          arma::mat &k_pp = Kmat_pp(iso, mp, n_delta, n_beta);
          for (INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rz_pp.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp,n_beta)-1,n_delta,index_z) = r_pp * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kz_pp.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp,n_beta)-1,n_delta,index_z) = k_pp * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_beta
      }//n_delta
    }//mp
    for (INT mp = 1 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
      {
        arma::cube &rz_mm = Rphiz_mm(iso,mp,n_beta);
        arma::cube &kz_mm = Kphiz_mm(iso,mp,n_beta);
        rz_mm = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.nMax(mp),nGHE);
        kz_mm = arma::zeros(basis.dMax * basis.n_zMax(mp,n_beta),basis.nMax(mp),nGHE);
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &r_mm = Rmat_mm(iso, mp, n_delta, n_beta);
          arma::mat &k_mm = Kmat_mm(iso, mp, n_delta, n_beta);
          for (INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rz_mm.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp,n_beta)-1,n_delta,index_z) = r_mm * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kz_mm.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp,n_beta)-1,n_delta,index_z) = k_mm * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_beta
      }//n_delta
    }//mp
  }//iso

  Multi<arma::vec> Rphi_pp; // Rphi_pp(iso,index_r,index_z)(index_gamma).
  Multi<arma::vec> Rphi_mm; // Rphi_mm(iso,index_r,index_z)(index_gamma).
  Multi<arma::vec> Kphi_pp; // Rphi_pp(iso,index_r,index_z)(index_gamma).
  Multi<arma::vec> Kphi_mm; // Rphi_mm(iso,index_r,index_z)(index_gamma).
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      for (INT index_r = 0 ; index_r < nGLA ; index_r++)
      {
        arma::vec &rphi_pp = Rphi_pp(iso,index_r,index_z);
        arma::vec &rphi_mm = Rphi_mm(iso,index_r,index_z);
        arma::vec &kphi_pp = Kphi_pp(iso,index_r,index_z);
        arma::vec &kphi_mm = Kphi_mm(iso,index_r,index_z);
        rphi_pp = arma::zeros(fullIndex_pp);
        rphi_mm = arma::zeros(fullIndex_mm);
        kphi_pp = arma::zeros(fullIndex_pp);
        kphi_mm = arma::zeros(fullIndex_mm);
        INT count = 0;
        for (INT mp = 0 ; mp < basis.mMax ; mp++)
        {
          for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
          {
            arma::cube &rz_pp = Rphiz_pp(iso,mp,n_beta);
            arma::cube &kz_pp = Kphiz_pp(iso,mp,n_beta);
            rphi_pp.subvec(index_Z_pp(count),index_Z_pp(count+1)-1) = rz_pp.slice(index_z) * rValsExchange(mp,index_r);
            kphi_pp.subvec(index_Z_pp(count),index_Z_pp(count+1)-1) = kz_pp.slice(index_z) * rValsExchange(mp,index_r);
            count += 1;
          }//n_beta.
        }//mp.
        count = 0;
        for (INT mp = 1 ; mp < basis.mMax ; mp++)
        {
          for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
          {
            arma::cube &rz_mm = Rphiz_mm(iso,mp,n_beta);
            arma::cube &kz_mm = Kphiz_mm(iso,mp,n_beta);
            rphi_mm.subvec(index_Z_mm(count),index_Z_mm(count+1)-1) = rz_mm.slice(index_z) * rValsExchange(mp,index_r);
            kphi_mm.subvec(index_Z_mm(count),index_Z_mm(count+1)-1) = kz_mm.slice(index_z) * rValsExchange(mp,index_r);
            count += 1;
          }//n_beta
        }//mp
      }//index_r
    }//index_z
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////// Rphi_mp.
  Multi<arma::cube> Rphiz_mp; // Rphiz_mp(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  Multi<arma::cube> Kphiz_mp; // Kphiz_mp(iso,mp,n_beta)(VECTORISE(n_zbeta,d_beta),n_delta,index_z).
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = 1 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(mp-1) ; n_beta++)
      {
        arma::cube &rz_mp = Rphiz_mp(iso,mp,n_beta);
        arma::cube &kz_mp = Kphiz_mp(iso,mp,n_beta);
        rz_mp = arma::zeros(basis.dMax * basis.n_zMax(mp-1,n_beta),basis.nMax(mp),nGHE);
        kz_mp = arma::zeros(basis.dMax * basis.n_zMax(mp-1,n_beta),basis.nMax(mp),nGHE);
        for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
        {
          arma::mat &r_mp = Rmat_mp(iso, mp, n_delta, n_beta);
          arma::mat &k_mp = Kmat_mp(iso, mp, n_delta, n_beta);
          for (INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rz_mp.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp-1,n_beta)-1,n_delta,index_z) = r_mp * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kz_mp.subcube(0,n_delta,index_z,basis.dMax * basis.n_zMax(mp-1,n_beta)-1,n_delta,index_z) = k_mp * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_beta
      }//n_delta
    }//mp
  }//iso
  Multi<arma::vec> Rphi_mp; // Rphi_mp(iso,index_r,index_z)(index_gamma).
  Multi<arma::vec> Kphi_mp; // Kphi_mp(iso,index_r,index_z)(index_gamma).
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      for (INT index_r = 0 ; index_r < nGLA ; index_r++)
      {
        arma::vec &rphi_mp = Rphi_mp(iso,index_r,index_z);
        arma::vec &kphi_mp = Kphi_mp(iso,index_r,index_z);
        rphi_mp = arma::zeros(fullIndex_mp);
        kphi_mp = arma::zeros(fullIndex_mp);
        INT count = 0;
        for (INT mp = 1 ; mp < basis.mMax ; mp++)
        {
          for (INT n_beta = 0 ; n_beta < basis.nMax(mp-1) ; n_beta++)
          {
            arma::cube &rz_mp = Rphiz_mp(iso,mp,n_beta);
            arma::cube &kz_mp = Kphiz_mp(iso,mp,n_beta);
            rphi_mp.subvec(index_Z_mp(count),index_Z_mp(count+1)-1) = rz_mp.slice(index_z) * rValsExchange(mp,index_r);
            kphi_mp.subvec(index_Z_mp(count),index_Z_mp(count+1)-1) = kz_mp.slice(index_z) * rValsExchange(mp,index_r);
            count += 1;
          }//n_beta.
        }//mp.
      }//index_r
    }//index_z
  }//iso


  //////////////////////////////////////////////////////////////////////////////////////////////////////// Ri_pp and Ri_mm and Ri_mp.
  Multi<arma::cube> Rze_0; //Rze_0(iso, mp, n_mu)(n_eta,VECTORISE(n_zmu,d_mu),index_z)
  Multi<arma::cube> Kze_0; //Kze_0(iso, mp, n_mu)(n_eta,VECTORISE(n_zmu,d_mu),index_z)
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
    {
      for (INT n_mu = 0 ; n_mu < basis.nMax(abs(mp)) ; n_mu++)
      {
        arma::cube &rze_0 = Rze_0(iso, mp, n_mu);
        arma::cube &kze_0 = Kze_0(iso, mp, n_mu);
        rze_0 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp),n_mu),nGHE);
        kze_0 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp),n_mu),nGHE);

        for (INT n_eta = 0 ; n_eta < basis.nMax(abs(mp)) ; n_eta++)
        {
          arma::mat &rme_0 = RmatExchange_0(iso, mp, n_mu, n_eta);
          arma::mat &km_0 = Kmat_0(iso, mp, n_mu, n_eta);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rze_0.subcube(n_eta,0,index_z,n_eta,basis.dMax * basis.n_zMax(abs(mp),n_mu)-1,index_z) = (arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_eta)-1,basis.dMax-1)).t() * rme_0).t();
            kze_0.subcube(n_eta,0,index_z,n_eta,basis.dMax * basis.n_zMax(abs(mp),n_mu)-1,index_z) = (arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_eta)-1,basis.dMax-1)).t() * km_0).t();
          }//index_z
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  //////////////////////////////////Builds Rze_1.
  Multi<arma::cube> Rze_1; //Rze_1(iso, mp, n_beta)(n_delta,VECTORISE(n_zbeta,d_beta),index_z)
  Multi<arma::cube> Kze_1; //Kze_1(iso, mp, n_beta)(n_delta,VECTORISE(n_zbeta,d_beta),index_z)
  for(INT iso: {NEUTRON, PROTON})
  {
    for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
    {
      for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
      {
        arma::cube &rze_1 = Rze_1(iso, mp, n_beta);
        arma::cube &kze_1 = Kze_1(iso, mp, n_beta);
        rze_1 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp-1),n_beta),nGHE);
        kze_1 = arma::zeros(basis.nMax(abs(mp)),basis.dMax * basis.n_zMax(abs(mp-1),n_beta),nGHE);
        for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
        {
          arma::mat &rme_1 = RmatExchange_1(iso, mp, n_delta, n_beta);
          arma::mat &km_1 = Kmat_1(iso, mp, n_delta, n_beta);
          for(INT index_z = 0 ; index_z < nGHE ; index_z++)
          {
            rze_1.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp-1),n_beta)-1,index_z) = rme_1 * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
            kze_1.subcube(n_delta,0,index_z,n_delta,basis.dMax * basis.n_zMax(abs(mp-1),n_beta)-1,index_z) = km_1  * arma::vectorise(waveMatZ(index_z).submat(0,0,basis.n_zMax(abs(mp),n_delta)-1,basis.dMax-1));
          }//index_z
        }//n_delta
      }//n_beta
    }//mp
  }//iso

  ////////////////////////////////////////////////////////////////////////////////////////////Builds Rzze and Kzze exchange 0.
  Multi<arma::cube> Rzze_0; //Rzze_0(iso, VECTORISE(n_zgamma,d_gamma),mp)(n_mu,n_eta,index_z).
  Multi<arma::cube> Kzze_0; //Kzze_0(iso, VECTORISE(n_zgamma,d_gamma),mp)(n_mu,n_eta,index_z).
  for(INT iso: {NEUTRON, PROTON})
  {
    INT count_z = 0;
    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
      {
        for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_0 = Rzze_0(iso,count_z,mp);
          arma::cube &kzze_0 = Kzze_0(iso,count_z,mp);
          rzze_0 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp)),nGHE);
          kzze_0 = arma::zeros(basis.nMax(abs(mp)),basis.nMax(abs(mp)),nGHE);

          for (INT n_mu = 0 ; n_mu < basis.nMax(abs(mp)) ; n_mu++)
          {
            arma::cube &rze_0 = Rze_0(iso, mp, n_mu);
            arma::cube &kze_0 = Kze_0(iso, mp, n_mu);
            for(INT index_z = 0 ; index_z < nGHE ; index_z++)
            {
              rzze_0.subcube(n_mu,0,index_z,n_mu,basis.nMax(abs(mp))-1,index_z) = rze_0.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_mu)-1,basis.dMax-1));
              kzze_0.subcube(n_mu,0,index_z,n_mu,basis.nMax(abs(mp))-1,index_z) = kze_0.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp),n_mu)-1,basis.dMax-1));
            }//index_z
          }//n_eta
        }//mp
        count_z += 1;
      }//n_zgamma
    }//d_gamma
  }//iso

  Multi<arma::cube> Rzze_1; //Rzze_1(iso, VECTORISE(n_zgamma,d_gamma),mp)(n_beta,n_delta,index_z).
  Multi<arma::cube> Kzze_1; //Kzze_1(iso, VECTORISE(n_zgamma,d_gamma),mp)(n_beta,n_delta,index_z).
  for(INT iso: {NEUTRON, PROTON})
  {
    INT count_z = 0;
    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
      {
        for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
        {
          arma::cube &rzze_1 = Rzze_1(iso, count_z,mp);
          arma::cube &kzze_1 = Kzze_1(iso, count_z,mp);
          rzze_1 = arma::zeros(basis.nMax(abs(mp-1)),basis.nMax(abs(mp)),nGHE);
          kzze_1 = arma::zeros(basis.nMax(abs(mp-1)),basis.nMax(abs(mp)),nGHE);
          for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp-1)) ; n_beta++)
          {
            arma::cube &rze_1 = Rze_1(iso, mp, n_beta);
            arma::cube &kze_1 = Kze_1(iso, mp, n_beta);
            for(INT index_z = 0 ; index_z < nGHE ; index_z++)
            {
              rzze_1.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = rze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1));
              kzze_1.subcube(n_beta,0,index_z,n_beta,basis.nMax(abs(mp))-1,index_z) = kze_1.slice(index_z) * arma::vectorise(IzExchange(d_gamma,n_zgamma,index_z).submat(0,0,basis.n_zMax(abs(mp-1),n_beta)-1,basis.dMax-1));
            }//index_z
          }//n_beta
        }//mp
        count_z += 1;
      }//n_zgamma
    }//d_gamma
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////
  Multi<arma::mat> Rzzre_0; //Rzzre_0(iso,index_z,VECTORISE(n_zgamma,d_gamma))(VECTORISE(mp,n_eta),index_r).
  Multi<arma::mat> Kzzre_0; //Kzzre_0(iso,index_z,VECTORISE(n_zgamma,d_gamma))(VECTORISE(mp,n_eta),index_r).
  for(INT iso: {NEUTRON, PROTON})
  {
    for(INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      INT count_z = 0;
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
        {
          arma::mat &r_0 = Rzzre_0(iso,index_z,count_z);
          arma::mat &k_0 = Kzzre_0(iso,index_z,count_z);
          r_0 = arma::zeros(index_mn_0,nGLA);
          k_0 = arma::zeros(index_mn_0,nGLA);
          for (INT mp = - basis.mMax + 1 ; mp < basis.mMax ; mp++)
          {
            arma::cube &rzze_0 = Rzze_0(iso,count_z,mp);
            arma::cube &kzze_0 = Kzze_0(iso,count_z,mp);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              r_0.submat(index_M_0(mp - (- basis.mMax + 1 )) ,index_r,index_M_0(mp+1 - (- basis.mMax + 1) )-1 ,index_r) = rzze_0.slice(index_z) * rValsExchange(abs(mp),index_r);
              k_0.submat(index_M_0(mp - (- basis.mMax + 1 )) ,index_r,index_M_0(mp+1 - (- basis.mMax + 1) )-1 ,index_r) = kzze_0.slice(index_z) * rValsExchange(abs(mp),index_r);
            }//index_r
          }//mp
          count_z += 1;
        }//n_zgamma
      }//d_gamma
    }//index_z
  }//iso

  /////////////////////////////////////////////////////////////////////////////////////////////////
  Multi<arma::mat> Rzzre_1; //Rzzre_1(iso, index_z,VECTORISE(n_zgamma,d_gamma))(VECTORISE(mp,n_delta),index_r).
  Multi<arma::mat> Kzzre_1; //Kzzre_1(iso, index_z,VECTORISE(n_zgamma,d_gamma))(VECTORISE(mp,n_delta),index_r).
  for(INT iso: {NEUTRON, PROTON})
  {
    for(INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      INT count_z = 0;
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
        {
          arma::mat &r_1 = Rzzre_1(iso, index_z,count_z);
          arma::mat &k_1 = Kzzre_1(iso, index_z,count_z);

          r_1 = arma::zeros(index_mn_1_p1,nGLA);
          k_1 = arma::zeros(index_mn_1_p1,nGLA);
          for (INT mp = - basis.mMax + 2 ; mp < basis.mMax ; mp++)
          {
            arma::cube &rzze_1 = Rzze_1(iso, count_z,mp);
            arma::cube &kzze_1 = Kzze_1(iso, count_z,mp);
            for(INT index_r = 0 ; index_r < nGLA ; index_r++)
            {
              r_1.submat(index_M_1_p1(mp - (- basis.mMax + 2 )) ,index_r,index_M_1_p1(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = rzze_1.slice(index_z) * rValsExchange(abs(mp),index_r);
              k_1.submat(index_M_1_p1(mp - (- basis.mMax + 2 )) ,index_r,index_M_1_p1(mp+1 - (- basis.mMax + 2) )-1 ,index_r) = kzze_1.slice(index_z) * rValsExchange(abs(mp),index_r);
            }//index_r
          }//mp
          count_z += 1;
        }//n_zgamma
      }//d_gamma
    }//index_z
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////
  Multi<arma::mat> Rifinal_0; //Rifinal_0(iso,index_z)(index_gamma,index_r).
  Multi<arma::mat> Kifinal_0; //Kifinal_0(iso,index_z)(index_gamma,index_r).
  Multi<arma::mat> Rifinal_3; //Rifinal_3(iso,index_z)(index_gamma,index_r).
  Multi<arma::mat> Kifinal_3; //Kifinal_3(iso,index_z)(index_gamma,index_r).

  for(INT iso: {NEUTRON, PROTON})
  {
    for(INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      arma::mat &ri_0 = Rifinal_0(iso,index_z);
      arma::mat &ki_0 = Kifinal_0(iso,index_z);
      ri_0 = arma::zeros(fullIndex_pp,nGLA);
      ki_0 = arma::zeros(fullIndex_pp,nGLA);
      INT count = 0;
      for (INT m = 0 ; m < basis.mMax ; m++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
        {
          INT count_z = 0;
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            count_z = d_gamma * basis.n_zGlobalMax;
            for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
            {
              arma::mat &r_0 = Rzzre_0(iso, index_z,count_z);
              arma::mat &k_0 = Kzzre_0(iso, index_z,count_z);
              for(INT index_r = 0 ; index_r < nGLA ; index_r++)
              {
                ri_0(count,index_r) = arma::accu(r_0.col(index_r) % IrExchange(index_r,m,n_gamma));
                ki_0(count,index_r) = arma::accu(k_0.col(index_r) % IrExchange(index_r,m,n_gamma));
              }//index_r
              count += 1;
              count_z += 1;
            }//n_zgamma
          }//d_gamma
        }//n_gamma
      }//m
      arma::mat &ri_3 = Rifinal_3(iso,index_z);
      arma::mat &ki_3 = Kifinal_3(iso,index_z);
      ri_3 = arma::zeros(fullIndex_mm,nGLA);
      ki_3 = arma::zeros(fullIndex_mm,nGLA);
      count = 0;
      for (INT m = 1 ; m < basis.mMax ; m++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
        {
          INT count_z = 0;
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            count_z = d_gamma * basis.n_zGlobalMax;
            for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m,n_gamma) ; n_zgamma++)
            {
              arma::mat &r_0 = Rzzre_0(iso, index_z,count_z);
              arma::mat &k_0 = Kzzre_0(iso, index_z,count_z);
              for(INT index_r = 0 ; index_r < nGLA ; index_r++)
              {
                ri_3(count,index_r) = arma::accu(r_0.col(index_r) % IrExchange_3(index_r,m,n_gamma));
                ki_3(count,index_r) = arma::accu(k_0.col(index_r) % IrExchange_3(index_r,m,n_gamma));
              }//index_r
              count += 1;
              count_z += 1;
            }//n_zgamma
          }//d_gamma
        }//n_gamma
      }//m
    }//index_z
  }//iso

  Multi<arma::mat> Rifinal_1; //Rifinal_1(iso,index_z)(index_gamma,index_r).
  Multi<arma::mat> Kifinal_1; //Kifinal_1(iso,index_z)(index_gamma,index_r).
  for(INT iso: {NEUTRON, PROTON})
  {
    for(INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      arma::mat &ri_1 = Rifinal_1(iso,index_z);
      arma::mat &ki_1 = Kifinal_1(iso,index_z);
      ri_1 = arma::zeros(fullIndex_mp,nGLA);
      ki_1 = arma::zeros(fullIndex_mp,nGLA);
      INT count = 0;
      for (INT m = 1 ; m < basis.mMax ; m++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m-1) ; n_gamma++)
        {
          INT count_z = 0;
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            count_z = d_gamma * basis.n_zGlobalMax;
            for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m-1,n_gamma) ; n_zgamma++)
            {
              arma::mat &r_1 = Rzzre_1(iso, index_z,count_z);
              arma::mat &k_1 = Kzzre_1(iso, index_z,count_z);
              for(INT index_r = 0 ; index_r < nGLA ; index_r++)
              {
                ri_1(count,index_r) = arma::accu(r_1.col(index_r) % IrExchange_1_p1(index_r,m,n_gamma));
                ki_1(count,index_r) = arma::accu(k_1.col(index_r) % IrExchange_1_p1(index_r,m,n_gamma));
              }//index_r
              count += 1;
              count_z += 1;
            }//n_zgamma
          }//d_gamma
        }//n_gamma
      }//m
    }//index_z
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////// Builds Rfinal.
  arma::mat Rfinal; //Rfinal(r,z).
  Rfinal = - (2*h3 + m3) * localRho % RphiDirect(NEUTRON) % RiDirect(NEUTRON) - (2*h3 + m3) * localRho % RphiDirect(PROTON) % RiDirect(PROTON)
  + (2 * w3 + b3) * localRho % (RphiDirect(NEUTRON) + RphiDirect(PROTON)) % (RiDirect(NEUTRON) + RiDirect(PROTON));
  Rfinal *= x;


  arma::mat Rfinal_Exchange = arma::zeros(nGLA,nGHE);
  arma::mat Kfinal = arma::zeros(nGLA,nGHE);
  for(INT iso: {NEUTRON, PROTON})
  {
    for(INT index_z = 0 ; index_z < nGHE ; index_z++)
    {
      arma::mat &ri_0 = Rifinal_0(iso,index_z);
      arma::mat &ri_1 = Rifinal_1(iso,index_z);
      arma::mat &ri_3 = Rifinal_3(iso,index_z);
      arma::mat &ki_0 = Kifinal_0(iso,index_z);
      arma::mat &ki_1 = Kifinal_1(iso,index_z);
      arma::mat &ki_3 = Kifinal_3(iso,index_z);
      for(INT index_r = 0 ; index_r < nGLA ; index_r++)
      {
        Rfinal_Exchange(index_r,index_z) += arma::accu(Rphi_pp(iso,index_r,index_z) % ri_0.col(index_r));
        Rfinal_Exchange(index_r,index_z) += arma::accu(Rphi_mm(iso,index_r,index_z) % ri_3.col(index_r));
        Rfinal_Exchange(index_r,index_z) += 2 * arma::accu(Rphi_mp(iso,index_r,index_z) % ri_1.col(index_r));

        Kfinal(index_r,index_z) += arma::accu(Kphi_pp(iso,index_r,index_z) % ki_0.col(index_r));
        Kfinal(index_r,index_z) += arma::accu(Kphi_mm(iso,index_r,index_z) % ki_3.col(index_r));
        Kfinal(index_r,index_z) += 2 * arma::accu(Kphi_mp(iso,index_r,index_z) % ki_1.col(index_r));
      }//index_r
    }//index_z
  }//iso

  Rfinal_Exchange = x * Rfinal_Exchange % localRho;
  Kfinal = x * Kfinal % localRho;
  ////////////////////////////////////////////////////////////////////////////////////////////////////////Numerical integration(z).
  Multi<arma::vec> Rz_Exchange;  //Rz_Exchange(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)
  Multi<arma::vec> Kz;  //Kz(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)
  Multi<arma::vec> Rz;  //Rz(n_zalpha,n_zgamma,d_alpha+2*d_gamma)(index_r)
  for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
  {
    for (INT n_zgamma = 0 ; n_zgamma < n_zalpha+1 ; n_zgamma++)
    {
      for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
      {
        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          Rz(n_zalpha,n_zgamma,d_alpha+2*d_gamma) = Rfinal * preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz(n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz_Exchange(n_zalpha,n_zgamma,d_alpha+2*d_gamma) = Rfinal_Exchange * preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Rz_Exchange(n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Rz_Exchange(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Kz(n_zalpha,n_zgamma,d_alpha+2*d_gamma) =  Kfinal * preZ(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
          Kz(n_zgamma,n_zalpha,d_gamma+2*d_alpha) = Kz(n_zalpha,n_zgamma,d_alpha+2*d_gamma);
        }//d_gamma
      }//d_alpha
    }//n_zgamma
  }//n_zalpha

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Numerical Integration (r).
  Multi<arma::cube> Delta; //Delta(m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  Multi<arma::cube> Gamma_Exchange; //Gamma_Exchange(m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  Multi<arma::cube> Gamma; //Gamma(m, n_alpha, n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha+1 ; n_gamma++)
      {
        arma::vec &prer = preR(m, n_alpha, n_gamma);

        arma::cube &gamma = Gamma(m, n_alpha, n_gamma);
        arma::cube &gamma_sym = Gamma(m, n_gamma, n_alpha);

        arma::cube &gammae = Gamma_Exchange(m, n_alpha, n_gamma);
        arma::cube &gammae_sym = Gamma_Exchange(m, n_gamma, n_alpha);

        arma::cube &delta = Delta(m, n_alpha, n_gamma);
        arma::cube &delta_sym = Delta(m, n_gamma, n_alpha);

        gamma = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        gamma_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

        gammae = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        gammae_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

        delta = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        delta_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
        {
          for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
          {
            for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
            {
              for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
              {
                gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(prer % Rz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                gammae(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(prer % Rz_Exchange(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                delta(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(prer % Kz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
              }//n_zgamma
            }//n_zalpha
            gamma_sym.slice(d_gamma+2*d_alpha) = gamma.slice(d_alpha+2*d_gamma).t();
            gammae_sym.slice(d_gamma+2*d_alpha) = gammae.slice(d_alpha+2*d_gamma).t();
            delta_sym.slice(d_gamma+2*d_alpha) = delta.slice(d_alpha+2*d_gamma).t();
          }//d_gamma
        }//d_alpha
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds the interaction..
  arma::mat field_direct = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  arma::mat field_exchange = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  arma::mat field_pairing = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

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
        field_direct.submat(bras.filter, kets.filter) = 2 * Gamma(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field_exchange.submat(bras.filter, kets.filter) = 2 * Gamma_Exchange(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field_pairing.submat(bras.filter, kets.filter) = 2 * Delta(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean Exchange part --.
      {
        field_direct.submat(bras.filter, kets.filter) = 2 * Gamma(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field_exchange.submat(bras.filter, kets.filter) = 2 * Gamma_Exchange(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field_pairing.submat(bras.filter, kets.filter) =  2 * Delta(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if
    }//j
  }//i

  for(INT iso: {NEUTRON, PROTON})
  {
    field(iso,DIRECT  ) = field_direct;
    field(iso,EXCHANGE) = field_exchange + field_pairing;
  }

  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangementFR::calcIz(void)
{
  DBG_ENTER;

  if (!IzDirect.empty()) DBG_LEAVE;

  //Dependencies.
  basis.calcTalmanz();

  //Useful tools.
  double p3 = parameters["p"];
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
          gphi(n_za) = d2_factor * z_factor * exp_factor * basis.zPartScalar((value_z + k_db)/sqrt(Kz), n_za) / std::pow(sqrt(Kz),n_za + 1);
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

  //Builds IzExchange.
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

void FieldRearrangementFR::calcIr(void)
{
  DBG_ENTER;

  if (!IrDirect.empty()) DBG_LEAVE;

  //Dependencies.
  basis.calcTalmanr();

  //Useful tools.
  double p3 = parameters["p"];
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
        }//mp
      }//n_gamma
    }//m
  }//index_r

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangementFR::calcWaveZ(void)
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

void FieldRearrangementFR::calcWaveMatZ(void)
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

void FieldRearrangementFR::calcRVals(void)
{
  DBG_ENTER;

  if (!rVals.empty()) DBG_LEAVE;

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

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangementFR::calcPreR(void)
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

void FieldRearrangementFR::calcPreZ(void)
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


//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangementFR::calcIndexZ(void)
{
  DBG_ENTER;

  if (!index_Z_pp.empty()) DBG_LEAVE;

  ////////////////////////////////////////Calc ++.
  INT length_pp = 0;
  for(INT mp = 0 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      length_pp +=1;
    }//n
  }//mp

  index_Z_pp = arma::zeros<IVEC>(length_pp+1);
  index_Z_pp(0) = 0;
  INT length = 0;
  INT lenZ = 0;
  for(INT mp = 0 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      lenZ += basis.dMax * basis.n_zMax(mp,n);
      index_Z_pp(length + 1) = lenZ;
      length += 1;
    }//n
  }//mp
  fullIndex_pp = lenZ;

  ////////////////////////////////////////Calc --.
  INT length_mm = 0;
  for(INT mp = 1 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      length_mm +=1;
    }//n
  }//mp

  index_Z_mm = arma::zeros<IVEC>(length_mm+1);
  index_Z_mm(0) = 0;
  length = 0;
  lenZ = 0;
  for(INT mp = 1 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      lenZ += basis.dMax*basis.n_zMax(mp,n);
      index_Z_mm(length + 1) = lenZ;
      length += 1;
    }//n
  }//mp
  fullIndex_mm = lenZ;

  ////////////////////////////////////////Calc -+.
  INT length_mp = 0;
  for(INT mp = 1 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp-1) ; n++)
    {
      length_mp +=1;
    }//n
  }//mp

  index_Z_mp = arma::zeros<IVEC>(length_mp+1);
  index_Z_mp(0) = 0;
  length = 0;
  lenZ = 0;
  for(INT mp = 1 ; mp < basis.mMax ; mp++)
  {
    for(INT n = 0 ; n < basis.nMax(mp-1) ; n++)
    {
      lenZ += basis.dMax*basis.n_zMax(mp-1,n);
      index_Z_mp(length + 1) = lenZ;
      length += 1;
    }//n
  }//mp
  fullIndex_mp = lenZ;


  DBG_LEAVE;
}
