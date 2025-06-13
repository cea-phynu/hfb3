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

#include "field_coulomb_exact.h"
#include "axis.h"
#include "interaction.h"
#include "multi.h"
#include "tools.h"

/** \file
 *  \brief Methods of the FieldCoulombExact class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldCoulombExact::FieldCoulombExact(Field::Parameters fp, State *_state) :
  Field(fp, _state),
  nQuad(10),
  Izx({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE})
{
  DBG_ENTER;

  name      = "coulomb (exact)";
  shortName = "CoulEx";

  nGLE      = 40;
  nGLA      = 40;
  nGHE      = 150;

  Axis axis_gl(Axis::GAUSS_LEGENDRE, nQuad, -1.0, 1.0);
  axis_exact = axis_gl;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldCoulombExact::calcField(void)
{
  DBG_ENTER;

  calculatingLength = -1.0;

  calcDirectField();

  if (!field(NEUTRON, EXCHANGE).empty() && !field(PROTON, EXCHANGE).empty()) DBG_LEAVE;

  double startTime = Tools::clock();

  //Dependencies
  calcIzx();
  calcIrx();

  INT dMax = basis.dMax;
  INT DMax = (dMax - 1) * 3 + 1;
  INT M = basis.mMax;
  INT N_z = basis.n_zGlobalMax;
  // INT N_p = basis.nGlobalMax;
  double fac_global = 2 * HBARC * ALPHA / sqrt(PI);

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Kappa and Rho Multis with explicit quantum numbers.
  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  HOqn.calcBlocks({0, 1, 3, 4});
  //Rho(m,mp,n,np,d+2*dp,s+2*sp)(n_z,n_zp)
  FMulti<arma::mat> Rho({MPRANGE, MPRANGE, NRANGE, NRANGE, DDRANGE, SSRANGE});

  //Kappa(m,mp,n,np,d+2*dp,s+2*sp)(n_z,n_zp)
  FMulti<arma::mat> Kappa({MPRANGE, MPRANGE, NRANGE, NRANGE, DDRANGE, SSRANGE});


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
        Rho(m, mp, n, np, d + 2 * dp, s + 2 * sp) = fac_global * rho(PROTON ).submat(bras.filter, kets.filter);
        Kappa(m, mp, n, np, d + 2 * dp, s + 2 * sp) = fac_global * kappa(PROTON ).submat(bras.filter, kets.filter);
      }//if
    }//j
  }//i

  // Multi<arma::cube> Rmat_0; //Rmat for the mean coulomb exchange term when m_alpha = m_gamma.    Rmat_0(m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)
  // Multi<arma::cube> Rmat_1; //Rmat for the mean coulomb exchange term when m_alpha = m_gamma +1. Rmat_1(m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)
  // Multi<arma::cube> Kmat_0; //Kmat for the coulomb pairing term when m_alpha = m_gamma.    Kmat_0(m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)
  // Multi<arma::cube> Kmat_1; //Kmat for the coulomb pairing term when m_alpha = m_gamma +1. Kmat_1(m_delta,n_delta,n_beta)(n_zdelta,n_zbeta,d_delta+2*d_beta)

  FMulti<arma::cube> Rmat_0({MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Rmat_1({MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Kmat_0({MRANGE, NRANGE, NRANGE});
  FMulti<arma::cube> Kmat_1({MRANGE, NRANGE, NRANGE});

  //Builds Rmat_0 and Kmat_0.
  for (INT n_delta = 0 ; n_delta < basis.nMax(0) ; n_delta++)
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(0) ; n_beta++)
    {
      arma::cube &rmat = Rmat_0(0, n_delta, n_beta);
      rmat = arma::zeros(basis.n_zMax(0, n_delta), basis.n_zMax(0, n_beta), DMax);

      arma::cube &kmat = Kmat_0(0, n_delta, n_beta);
      kmat = arma::zeros(basis.n_zMax(0, n_delta), basis.n_zMax(0, n_beta), DMax);

      for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
      {
        for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
        {
          rmat.slice(d_delta + 2 * d_beta) = -Rho(0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0);
          kmat.slice(d_delta + 2 * d_beta) = Kappa(0, 0, n_delta, n_beta, d_delta + 2 * d_beta, 0);
        }//d_delta
      }//d_beta
    }//n_delta
  }//n_beta

  for (INT mp = 1 ; mp < M ; mp++)
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(mp) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp); n_delta++)
      {
        arma::cube &rmat   = Rmat_0( mp, n_delta, n_beta);
        rmat   = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);
        arma::cube &rmat_m = Rmat_0(-mp, n_delta, n_beta);
        rmat_m = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);

        arma::cube &kmat   = Kmat_0( mp, n_delta, n_beta);
        kmat   = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);
        arma::cube &kmat_m = Kmat_0(-mp, n_delta, n_beta);
        kmat_m = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp, n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            rmat.slice(d_delta + 2 * d_beta)   = -Rho(mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0);
            rmat_m.slice(d_delta + 2 * d_beta) = -Rho(mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3);
            kmat.slice(d_delta + 2 * d_beta)   = Kappa(mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 0);
            kmat_m.slice(d_delta + 2 * d_beta) = Kappa(mp, mp, n_delta, n_beta, d_delta + 2 * d_beta, 3);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp

  //Builds Rmat_1 and Kmat_1.
  for (INT mp = -M + 2 ; mp < 1 ; mp++) //mp <= 0.
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
      {
        arma::cube &rmat = Rmat_1(mp, n_delta, n_beta);
        rmat = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);
        arma::cube &kmat = Kmat_1(mp, n_delta, n_beta);
        kmat = arma::zeros(basis.n_zMax(abs(mp), n_delta), basis.n_zMax(abs(mp - 1), n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            rmat.slice(d_delta + 2 * d_beta) = Rho(abs(mp), abs(mp - 1), n_delta, n_beta, d_delta + 2 * d_beta, 2);
            kmat.slice(d_delta + 2 * d_beta) = -Kappa(abs(mp), abs(mp - 1), n_delta, n_beta, d_delta + 2 * d_beta, 2);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp <= 0

  for (INT mp = 1 ; mp < M ; mp++) //mp > 0.
  {
    for (INT n_beta = 0 ; n_beta < basis.nMax(mp - 1) ; n_beta++)
    {
      for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
      {
        arma::cube &rmat = Rmat_1(mp, n_delta, n_beta);
        rmat = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);

        arma::cube &kmat = Kmat_1(mp, n_delta, n_beta);
        kmat = arma::zeros(basis.n_zMax(mp, n_delta), basis.n_zMax(mp - 1, n_beta), DMax);

        for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
        {
          for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
          {
            rmat.slice(d_delta + 2 * d_beta) = -Rho(mp, mp - 1, n_delta, n_beta, d_delta + 2 * d_beta, 1);
            kmat.slice(d_delta + 2 * d_beta) = Kappa(mp, mp - 1, n_delta, n_beta, d_delta + 2 * d_beta, 1);
          }//d_delta
        }//d_beta
      }//n_delta
    }//n_beta
  }//mp > 0

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Calculates IzR and IzK merging Iz with Rmat and Kmat.
  // Multi<arma::vec>  IzR_0; //RmatExchange_0 with IzExchange. IzRExchange_0(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_beta + n_delta*basis.nMax(mp) + indexM(mp));
  // Multi<arma::vec>  IzR_1; //RmatExchange_1 with IzExchange. IzRExchange_1(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index3M(mp+M-2));
  // Multi<arma::vec>  IzR_3; //RmatExchange_0 with IzExchange in a reverse order. IzRExchange_3(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index2M(-mp+M-1));
  // Multi<arma::vec>  IzK_0; //Kmat_0 with IzExchange. IzK_0(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_beta + n_delta*basis.nMax(mp) + indexM(mp));
  // Multi<arma::vec>  IzK_1; //Kmat_1 with IzExchange. Izk_1(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index3M(mp+M-2));
  // Multi<arma::vec>  IzK_3; //Kmat_0 with IzExchange in a reverse order. IzK_3(x,d_alpha+2*d_gamma,n_zalpha,n_zgamma)(n_delta + n_beta*basis.nMax(abs(mp)) + index2M(-mp+M-1));

  FMulti<arma::vec>  IzR_0({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzR_1({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzR_3({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_0({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_1({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});
  FMulti<arma::vec>  IzK_3({{0, nQuad}, DDRANGE, NZRANGE, NZRANGE});

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

  arma::cube IzcubeE_0;
  arma::cube IzcubeE_1;

  //Calculates IzRExchange_0, IzRExchange_3, IzK_0 and IzK_3.
  for (INT x = 0; x < nQuad; x++)
  {
    for (INT n_zalpha = 0; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
    {
      for (INT n_zgamma = 0; n_zgamma < n_zalpha + 1 ; n_zgamma++) //One uses the symmetry gamma <-> alpha.
      {
        for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            INT D = d_alpha + 2 * d_gamma;
            INT D_sym = d_gamma + 2 * d_alpha;

            arma::vec &izr_0 = IzR_0(x, D, n_zalpha, n_zgamma);
            arma::vec &izr_0_sym = IzR_0(x, D_sym, n_zgamma, n_zalpha);
            izr_0 = arma::zeros(sizeInit2);
            izr_0_sym = arma::zeros(sizeInit2);

            arma::vec &izr_3 = IzR_3(x, D, n_zalpha, n_zgamma);
            arma::vec &izr_3_sym = IzR_3(x, D_sym, n_zgamma, n_zalpha);
            izr_3 = arma::zeros(sizeInit2);
            izr_3_sym = arma::zeros(sizeInit2);


            arma::vec &izk_0 = IzK_0(x, D, n_zalpha, n_zgamma);
            arma::vec &izk_0_sym = IzK_0(x, D_sym, n_zgamma, n_zalpha);
            izk_0 = arma::zeros(sizeInit2);
            izk_0_sym = arma::zeros(sizeInit2);

            arma::vec &izk_3 = IzK_3(x, D, n_zalpha, n_zgamma);
            arma::vec &izk_3_sym = IzK_3(x, D_sym, n_zgamma, n_zalpha);
            izk_3 = arma::zeros(sizeInit2);
            izk_3_sym = arma::zeros(sizeInit2);

            for (INT mp = -M + 1 ; mp < M ; mp++)
            {
              for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
              {
                for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
                {
                  IzcubeE_0 = Izx(x, d_alpha + 2 * d_gamma, n_zalpha, n_zgamma).subcube(0, 0, 0, basis.n_zMax(abs(mp), n_delta) - 1, basis.n_zMax(abs(mp), n_beta) - 1, DMax - 1);
                  INT N = n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1);
                  INT N_sym = n_beta + n_delta * basis.nMax(abs(mp)) + index2M(mp + M - 1);
                  izr_0(N) = arma::accu(IzcubeE_0 % Rmat_0(mp, n_delta, n_beta));
                  izr_0_sym(N_sym) = izr_0(N);

                  INT Np = n_delta + n_beta * basis.nMax(abs(mp)) + index2M(-mp + M - 1);
                  INT Np_sym = n_beta + n_delta * basis.nMax(abs(mp)) + index2M(-mp + M - 1);
                  izr_3(Np) = izr_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                  izr_3_sym(Np_sym) = izr_3(Np);

                  izk_0(N) = arma::accu(IzcubeE_0 % Kmat_0(mp, n_delta, n_beta));
                  izk_0_sym(N_sym) = izk_0(N);

                  izk_3(Np) = izk_0(n_delta + n_beta * basis.nMax(abs(mp)) + index2M(mp + M - 1));
                  izk_3_sym(Np_sym) = izk_3(Np);
                }//n_beta
              }//n_delta
            }//mp
          }//d_gamma
        }//d_alpha
      }//n_zgamma
    }//n_zalpha
  }//x

  //Calculates IzRExchange_1 and IzK_1.
  for (INT x = 0; x < nQuad; x++)
  {
    for (INT n_zalpha = 0; n_zalpha < N_z ; n_zalpha++)
    {
      for (INT n_zgamma = 0; n_zgamma < N_z ; n_zgamma++)
      {
        for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            INT D = d_alpha + 2 * d_gamma;

            arma::vec &izk_1 = IzK_1(x, D, n_zalpha, n_zgamma);
            izk_1 = arma::zeros(sizeInit3);

            arma::vec &izr_1 = IzR_1(x, D, n_zalpha, n_zgamma);
            izr_1 = arma::zeros(sizeInit3);

            for (INT mp = -M + 2 ; mp < M ; mp++)
            {
              for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
              {
                for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
                {
                  INT N = n_delta + n_beta * basis.nMax(abs(mp)) + index3M(mp + M - 2);
                  IzcubeE_1 = Izx(x, d_alpha + 2 * d_gamma, n_zalpha, n_zgamma).subcube(0, 0, 0, basis.n_zMax(abs(mp), n_delta) - 1, basis.n_zMax(abs(mp - 1), n_beta) - 1, DMax - 1);

                  izk_1(N) = arma::accu(IzcubeE_1 % Kmat_1(mp, n_delta, n_beta));

                  izr_1(N) = arma::accu(IzcubeE_1 % Rmat_1( mp, n_delta, n_beta));
                }//np
              }//n
            }//mp
          }//d_gamma
        }//d_alpha
      }//n_zgamma
    }//n_zalpha
  }//x

  ////////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Gamma_Direct, Gamma_Exchange and Delta from merging Ir with IzR and IzK.
  // Multi<arma::cube> Gamma; //GammaDirect(m,n_alpha,n_gamma,s_alpha+2*s_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).
  // Multi<arma::cube> Delta; // Delta(m,n_alpha,n_gamma,s_alpha+2*s_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma).

  FMulti<arma::cube> Gamma({MPRANGE, NZRANGE, NZRANGE, SSRANGE});
  FMulti<arma::cube> Delta({MPRANGE, NZRANGE, NZRANGE, SSRANGE});

  for (INT m = 0 ; m < M ; m++)
  {
    for (INT n_alpha = 0; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0; n_gamma < n_alpha + 1 ; n_gamma++) //One uses the symmetry gamma <-> alpha.
      {
        arma::cube &gamma_0 = Gamma(m, n_alpha, n_gamma, 0);
        arma::cube &gamma_3 = Gamma(m, n_alpha, n_gamma, 3);
        arma::cube &gamma_0_sym = Gamma(m, n_gamma, n_alpha, 0);
        arma::cube &gamma_3_sym = Gamma(m, n_gamma, n_alpha, 3);
        gamma_0 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        gamma_3 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        gamma_0_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        gamma_3_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);


        arma::cube &delta_0 = Delta( m, n_alpha, n_gamma, 0);
        arma::cube &delta_3 = Delta( m, n_alpha, n_gamma, 3);
        arma::cube &delta_0_sym = Delta( m, n_gamma, n_alpha, 0);
        arma::cube &delta_3_sym = Delta( m, n_gamma, n_alpha, 3);
        delta_0 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        delta_3 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        delta_0_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        delta_3_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);

        for (INT x = 0; x < nQuad ; x++)
        {
          arma::vec &irx = Irx(x, 0, m, n_alpha, n_gamma);

          for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
            {
              INT D = d_alpha + 2 * d_gamma;

              for (INT n_zalpha = 0; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
                {
                  gamma_0 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzR_0(x, D, n_zalpha, n_zgamma));
                  gamma_3 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzR_3(x, D, n_zalpha, n_zgamma));

                  delta_0 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzK_0(x, D, n_zalpha, n_zgamma));
                  delta_3 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzK_3(x, D, n_zalpha, n_zgamma));
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//x

        for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
          {
            INT D = d_alpha + 2 * d_gamma;
            INT D_bar = d_gamma + 2 * d_alpha;

            gamma_0_sym.slice(D_bar) = (gamma_0.slice(D)).t();
            gamma_3_sym.slice(D_bar) = (gamma_3.slice(D)).t();

            delta_0_sym.slice(D_bar)  = (delta_0.slice(D)).t();
            delta_3_sym.slice(D_bar)  = (delta_3.slice(D)).t();
          }//d_gamma
        }//d_alpha
      }//n_gamma
    }//n_alpha
  }//m

  for (INT m = 1 ; m < M ; m++)
  {
    for (INT n_alpha = 0; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0; n_gamma < basis.nMax(m - 1) ; n_gamma++)
      {
        arma::cube &gamma_1 = Gamma(m, n_alpha, n_gamma, 1);
        gamma_1 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);

        arma::cube &delta_1 = Delta(m, n_alpha, n_gamma, 1);
        delta_1 = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m - 1, n_gamma), DMax);

        for (INT x = 0; x < nQuad ; x++)
        {
          arma::vec &irx = Irx(x, 1, m, n_alpha, n_gamma);

          for (INT d_alpha = 0; d_alpha < dMax ; d_alpha++)
          {
            for (INT d_gamma = 0; d_gamma < dMax ; d_gamma++)
            {
              for (INT n_zalpha = 0; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
              {
                for (INT n_zgamma = 0; n_zgamma < basis.n_zMax(m - 1, n_gamma) ; n_zgamma++)
                {
                  INT D = d_alpha + 2 * d_gamma;
                  gamma_1 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzR_1(x, D, n_zalpha, n_zgamma));

                  delta_1 (n_zalpha, n_zgamma, D) += arma::accu(irx % IzK_1(x, D, n_zalpha, n_zgamma));
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//x
      }//n_gamma
    }//n_alpha
  }//m


  ////////////////////////////////////////////////////////////////////////////////////////////////////////Builds the interaction.

  field(NEUTRON, EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, EXCHANGE) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(NEUTRON, PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, PAIRING ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

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

      if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0) //Mean Exchange part ++ and Pairing part ++.
      {
        arma::mat &Gmat = Gamma(m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);

        field(PROTON, EXCHANGE).submat(bras.filter, kets.filter) = Gmat;

        arma::mat &Dmat = Delta(m_alpha, n_alpha, n_gamma, 0).slice(d_alpha + 2 * d_gamma);
        field(PROTON, PAIRING ).submat(bras.filter, kets.filter) = Dmat;
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean exchange part -- and Pairing part --.
      {
        arma::mat &Gmat = Gamma( m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
        field(PROTON, EXCHANGE).submat(bras.filter, kets.filter) = Gmat;

        arma::mat &Dmat = Delta( m_alpha, n_alpha, n_gamma, 3).slice(d_alpha + 2 * d_gamma);
        field(PROTON, PAIRING ).submat(bras.filter, kets.filter) = Dmat;
      }//if

      if (m_alpha == m_gamma + 1 && m_alpha > 0 && s_alpha == 1 && s_gamma == 0) //Mean exchange part -+ and Pairing part -+.
      {
        arma::mat &Gmat = Gamma( m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
        field(PROTON, EXCHANGE).submat(bras.filter, kets.filter) = Gmat;

        arma::mat &Dmat = Delta(m_alpha, n_alpha, n_gamma, 1).slice(d_alpha + 2 * d_gamma);
        field(PROTON, PAIRING ).submat(bras.filter, kets.filter) = Dmat;
      }//if

      if (m_alpha == m_gamma - 1 && m_alpha < M - 1 && s_alpha == 0 && s_gamma == 1) //Mean exchange part +- and Pairing part +-.
      {
        arma::mat &Gmat = Gamma( m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);
        field(PROTON, EXCHANGE).submat(bras.filter, kets.filter) = Gmat.t();

        arma::mat &Dmat = Delta( m_alpha + 1, n_gamma, n_alpha, 1).slice(d_gamma + 2 * d_alpha);
        field(PROTON, PAIRING ).submat(bras.filter, kets.filter) = Dmat.t();
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

/** Calculates Irx.
 */

void FieldCoulombExact::calcIrx(void)
{
  DBG_ENTER;

  // dependencies
  basis.calcTalmanr();
  basis.calcMoshinskyr();

  if (!Irx.empty()) DBG_LEAVE;

  INT M = basis.mMax;
  INT maxn = basis.Nmaxr;
  INT nl = 0.0;
  double fac = 0.0;

  /////////////////////////////////////////////////////////////Calculates the Gaussian integral factor and merges it with a Moshinsky.
  FMulti<arma::vec> mosGauss({{0, nQuad}, {-2*M+2, 2*M-1}, {0, maxn}});

  for (INT x = 0; x < nQuad; x++)
  {
    //fac =  pow(axis_exact.p(x) / basis.b_r, 2);
    double mu_x = (1 + axis_exact.p(x)) / (1 - axis_exact.p(x));
    fac =  pow(mu_x / basis.b_r, 2);

    for (INT ma = M * -2 + 2; ma < 2 * M - 1; ma++)
    {
      for (INT nb = 0; nb < maxn; nb++)
      {
        mosGauss(x, ma, nb) = arma::zeros(maxn);
        arma::vec &mos = mosGauss(x, ma, nb);

        for (INT na = 0; na < maxn; na++)
        {
          nl = na + nb + abs(ma);
          mos(na) = basis.moshinskyr(ma, nb)(na) * fac * (pow(2.0, nl) / pow(fac + 2.0, nl + 1));
        }//na
      }//nb
    }//ma
  }//x

  /////////////////////////////////////////////////////////////Calculates Jr.
  Multi<arma::vec> Jr; //Jr(x,m_alpha,n_alpha,m_beta,n_beta)(n_b).

  for (INT x = 0; x < nQuad; x++)
  {
    double fac_x = 2 * axis_exact.w(x) / std::pow(1 + axis_exact.p(x), 2);

    for (INT m_alpha = 0; m_alpha < M; m_alpha++)
    {
      for (INT n_alpha = 0; n_alpha < basis.nMax(m_alpha); n_alpha++)
      {
        for (INT m_beta = -M + 1; m_beta < M; m_beta++)
        {
          for (INT n_beta = 0; n_beta < basis.nMax(abs(m_beta)); n_beta++)
          {
            Jr(x, m_alpha, n_alpha, m_beta, n_beta) = arma::zeros(maxn);
            Jr(x, -m_alpha, n_alpha, -m_beta, n_beta) = arma::zeros(maxn);

            arma::vec &talvec  = basis.talmanr(m_alpha, n_alpha, m_beta, n_beta);
            arma::mat &jvec  = Jr(x, m_alpha, n_alpha, m_beta, n_beta);
            arma::mat &jvec_min = Jr(x, -m_alpha, n_alpha, -m_beta, n_beta);

            for (INT nb = 0 ; nb < maxn ; nb++)
            {
              jvec(nb) = fac_x * arma::accu(talvec % mosGauss(x, m_beta - m_alpha, nb));
              jvec_min(nb) = jvec(nb);
            }//nb
          }//n_beta
        }//m_beta
      }//n_alpha
    }//m_alpha
  }//x

  /////////////////////////////////////////////////////////////////////////////////////////////////////Calculates indices for the exchange part.
  INT dimExchange_0 = 0;
  UVEC indicesExchange_0 = arma::zeros<UVEC>(2 * M);

  for (INT m = -M + 1 ; m < M ; m++)
  {
    dimExchange_0 += basis.nMax(abs(m)) * basis.nMax(abs(m));
    indicesExchange_0(m + M) = dimExchange_0;
  }

  INT dimExchange_1 = 0;
  UVEC indicesExchange_1 = arma::zeros<UVEC>(2 * M - 1);

  for (INT m = -M + 2 ; m < M ; m++)
  {
    dimExchange_1 += basis.nMax(abs(m)) * basis.nMax(abs(m - 1));
    indicesExchange_1(m + M - 1) = dimExchange_1;
  }

  arma::mat MatExchange;
  arma::mat MatExchange_sym;

  Irx = FMulti<arma::vec>({{0, nQuad}, {0, 2}, MPRANGE, NRANGE, NRANGE});

  ////////////////////////////////////////////////////////////////////////////////////////Calculates IrExchange for m_alpha = m_gamma.
  for (INT x = 0; x < nQuad; x++)
  {
    for (INT m = 0 ; m < M  ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //One uses the symmetry gamma <-> alpha.
        {
          arma::vec &IrEvec =  Irx(x, 0, m, n_alpha, n_gamma);
          arma::vec &IrEvec_sym =  Irx(x, 0, m, n_gamma, n_alpha);
          IrEvec = arma::zeros(dimExchange_0);
          IrEvec_sym = arma::zeros(dimExchange_0);

          for (INT mp = -M + 1 ; mp < M; mp++)
          {
            MatExchange = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp)));
            MatExchange_sym = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp)));

            for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
            {
              arma::vec &jvec = Jr(x, m, n_alpha, mp, n_delta);

              for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp)) ; n_beta++)
              {
                MatExchange(n_delta, n_beta) = arma::accu(basis.talmanr(mp, n_beta, m, n_gamma) % jvec);
                MatExchange_sym(n_beta, n_delta) = MatExchange(n_delta, n_beta);
              }// n_belta
            }// n_delta

            IrEvec.subvec(indicesExchange_0(mp - 1 + M), indicesExchange_0(mp + M) - 1) = arma::vectorise(MatExchange);
            IrEvec_sym.subvec(indicesExchange_0(mp - 1 + M), indicesExchange_0(mp + M) - 1) = arma::vectorise(MatExchange_sym);
          }//mp
        }//n_gamma
      }//n_alpha
    }//m
  }//x

  //////////////////////////////////////////////////////////////////////////////////Calculates IrExchange for m_alpha = m_gamma + 1.
  for (INT x = 0; x < nQuad; x++)
  {
    for (INT m = 1 ; m < M  ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < basis.nMax(m - 1) ; n_gamma++)
        {
          arma::vec &IrEvec =  Irx(x, 1, m, n_alpha, n_gamma);
          IrEvec = arma::zeros(dimExchange_1);

          for (INT mp = -M + 2 ; mp < M; mp++)
          {
            MatExchange = arma::zeros(basis.nMax(abs(mp)), basis.nMax(abs(mp - 1)));

            for (INT n_delta = 0 ; n_delta < basis.nMax(abs(mp)) ; n_delta++)
            {
              arma::vec &jvec = Jr(x, m, n_alpha, mp, n_delta);

              for (INT n_beta = 0 ; n_beta < basis.nMax(abs(mp - 1)) ; n_beta++)
              {
                MatExchange(n_delta, n_beta) = arma::accu(basis.talmanr(mp - 1, n_beta, m - 1, n_gamma) % jvec);
              }//n_belta
            }//n_delta

            IrEvec.subvec(indicesExchange_1(mp - 2 + M), indicesExchange_1(mp + M - 1) - 1) = arma::vectorise(MatExchange);
          }//mp
        }//n_gamma
      }//n_alpha
    }//m
  }//x

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Izx.
 */

void FieldCoulombExact::calcIzx(void)
{
  DBG_ENTER;

  // dependencies
  basis.calcTalmanz();
  basis.calcMoshinskyz();

  if (!Izx.empty()) DBG_LEAVE;

  INT dMax = basis.dMax;
  INT n_zGlobalMax = basis.n_zGlobalMax;

  /////////////////////////////////////////////////////////////Gaussian integral factor's calculation.
  FMulti<arma::vec> gaussianIntegralz({{0, nQuad}, DRANGE, DRANGE, DRANGE, DRANGE});

  double fac0 = sqrt(sqrt(PI) * basis.b_z);

  for (INT x = 0; x < nQuad; x++)
  {
    double mu_x = (1 + axis_exact.p(x)) / (1 - axis_exact.p(x));
    double fac1 = mu_x / basis.b_z;

    for (INT d_alpha = 0; d_alpha < basis.dMax; d_alpha++)
    {
      double zd_alpha = (0.5 - d_alpha) * basis.d_0;

      for (INT d_beta = 0; d_beta < basis.dMax; d_beta++)
      {
        double zd_beta = (0.5 - d_beta) * basis.d_0;

        for (INT d_gamma = 0; d_gamma < basis.dMax; d_gamma++)
        {
          double zd_gamma = (0.5 - d_gamma) * basis.d_0;
          double d_alga = (zd_alpha + zd_gamma) / 2.0;

          for (INT d_delta = 0; d_delta < basis.dMax; d_delta++)
          {
            double zd_delta = (0.5 - d_delta) * basis.d_0;
            double d_bede = (zd_beta + zd_delta) / 2.0;
            gaussianIntegralz(x, d_alpha, d_beta, d_gamma, d_delta) = arma::zeros(basis.n_zGlobalMax * 4);
            arma::vec &gauvec = gaussianIntegralz(x, d_alpha, d_beta, d_gamma, d_delta);

            for (INT n_z = 0; n_z < basis.n_zGlobalMax * 4 - 3; n_z++)
            {
              double fac2 = sqrt(pow(2.0, n_z) / pow(2.0 + pow(fac1, 2), (n_z + 1)));
              double sigma = (d_bede - d_alga) / sqrt(2.0 + pow(fac1, 2));
              gauvec(n_z) =  fac0 * fac1 * fac2 * exp(-0.5 * pow(sigma / basis.b_z, 2)) * basis.zPartScalar(sigma, n_z);
            }//n_z
          }//d_delta
        }//d_gamma
      }//d_beta
    }//d_alpha
  }//x

  //////////////////////////////////////////////////////////////////////Merges a Moshinsky with the Gaussian integral factor.
  FMulti<arma::vec> mosvec({{0, nQuad}, DDRANGE, DDRANGE, {0, 2*basis.n_zGlobalMax}});

  for (INT x = 0; x < nQuad; x++)
  {
    for (INT d_alpha = 0; d_alpha < dMax; d_alpha++)
    {
      for (INT d_beta = 0; d_beta < dMax; d_beta++)
      {
        INT D_albe = d_alpha + 2 * d_beta;

        for (INT d_gamma = 0; d_gamma < dMax; d_gamma++)
        {
          for (INT d_delta = 0; d_delta < dMax; d_delta++)
          {
            INT D_gade = d_gamma + 2 * d_delta;
            arma::vec &gausVec = gaussianIntegralz(x, d_alpha, d_beta, d_gamma, d_delta);

            for (INT n_za = 0; n_za < 2 * n_zGlobalMax; n_za ++)
            {
              mosvec(x, D_albe, D_gade, n_za) = arma::zeros(2 * n_zGlobalMax);
              arma::vec &mos = mosvec(x, D_albe, D_gade, n_za);

              for (INT n_zb = 0; n_zb < 2 * n_zGlobalMax; n_zb ++)
              {
                mos(n_zb) = basis.moshinskyz(n_za, n_zb) * gausVec(n_za + n_zb);
              }//n_zb
            }//n_za
          }//d_delta
        }//d_gamma
      }//d_beta
    }//d_alpha
  }//x

  //////////////////////////////////////////////////////////////////Calculates Jz.

  //Jz(x,d_alpha+2*d_beta,d_gamma+2*d_delta,n_zbeta,n_zdelta)(n_z)
  FMulti<arma::vec> Jz({{0, nQuad}, DDRANGE, DDRANGE, NZRANGE, NZRANGE});

  for (INT d_beta = 0; d_beta < dMax; d_beta++)
  {
    for (INT d_delta = 0; d_delta < dMax; d_delta++)
    {
      for (INT n_zbeta = 0; n_zbeta < n_zGlobalMax; n_zbeta++)
      {
        for (INT n_zdelta = 0; n_zdelta < n_zbeta + 1; n_zdelta++) //One uses the symmetry beta <-> delta.
        {
          arma::vec talvec = basis.talmanz(n_zbeta, d_beta, n_zdelta, d_delta).subvec(0, 2 * n_zGlobalMax - 1);

          for (INT x = 0; x < nQuad ; x++)
          {
            for (INT d_alpha = 0; d_alpha < dMax; d_alpha++)
            {
              INT D_albe = d_alpha + 2 * d_beta;
              INT D_albe_sym = d_alpha + 2 * d_delta;

              for (INT d_gamma = 0; d_gamma < dMax; d_gamma++)
              {
                INT D_gade = d_gamma + 2 * d_delta;
                INT D_gade_sym = d_gamma + 2 * d_beta;
                Jz(x, D_albe, D_gade, n_zbeta, n_zdelta) = arma::zeros(basis.n_zGlobalMax * 2);
                Jz(x, D_albe_sym, D_gade_sym, n_zdelta, n_zbeta) = arma::zeros(basis.n_zGlobalMax * 2);
                arma::vec &jzvec_up = Jz(x, D_albe, D_gade, n_zbeta, n_zdelta);
                arma::vec &jzvec_low = Jz(x, D_albe_sym, D_gade_sym, n_zdelta, n_zbeta);

                for (INT n_za = 0; n_za < n_zGlobalMax * 2; n_za++)
                {
                  jzvec_up(n_za) = arma::accu(talvec % mosvec(x, D_albe, D_gade, n_za));
                  jzvec_low(n_za) = jzvec_up(n_za);
                }//n_za
              }//d_gamma
            }//d_alpha
          }//x
        }//n_zdelta
      }//n_zbeta
    }//d_delta
  }//d_beta

  INT DMax = (dMax - 1) * 3 + 1;;
  INT n_zGMax = basis.n_zGlobalMax;

  /////////////////////////////////////////////////////////////Calculates IzDirect.
  FMulti<arma::cube> IzDirect({{0, nQuad}, DRANGE, DRANGE, NZRANGE, NZRANGE});

  for (INT x = 0; x < nQuad ; x++)
  {
    for (INT d_alpha = 0; d_alpha < dMax; d_alpha++)
    {
      for (INT d_gamma = 0; d_gamma < dMax; d_gamma++)
      {
        for (INT n_zalpha = 0; n_zalpha < n_zGMax; n_zalpha++)
        {
          for (INT n_zgamma = 0; n_zgamma < n_zalpha + 1; n_zgamma++) //One uses the symmetry gamma <-> alpha.
          {
            IzDirect(x, d_alpha, d_gamma, n_zalpha, n_zgamma) = arma::zeros(n_zGMax, n_zGMax, DMax);
            IzDirect(x, d_gamma, d_alpha, n_zgamma, n_zalpha) = arma::zeros(n_zGMax, n_zGMax, DMax);
            arma::cube &IzCubeD_up = IzDirect(x, d_alpha, d_gamma, n_zalpha, n_zgamma);
            arma::cube &IzCubeD_low = IzDirect(x, d_gamma, d_alpha, n_zgamma, n_zalpha);
            arma::vec talvec = basis.talmanz(n_zalpha, d_alpha, n_zgamma, d_gamma).subvec(0, 2 * n_zGMax - 1);

            for (INT d_delta = 0; d_delta < dMax; d_delta++)
            {
              INT D_gade = d_gamma + 2 * d_delta;

              for (INT n_zdelta = 0; n_zdelta < n_zGMax; n_zdelta++)
              {
                for (INT d_beta = 0; d_beta < dMax; d_beta++)
                {
                  INT Dd = 2 * d_beta + d_delta;
                  //INT De = 2 * d_beta + d_gamma;
                  INT Dd_sym = 2 * d_delta + d_beta;
                  INT D_albe = d_alpha + 2 * d_beta;

                  for (INT n_zbeta = 0; n_zbeta < n_zdelta + 1; n_zbeta++)
                  {
                    IzCubeD_up(n_zdelta, n_zbeta, Dd) = arma::accu(talvec % Jz(x, D_albe, D_gade, n_zbeta, n_zdelta));
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
  }//x

  /////////////////////////////////////////////////////////////Calculates IzExchange.
  for (INT x = 0; x < nQuad; x++)
  {
    for (INT d_alpha = 0; d_alpha < dMax; d_alpha++)
    {
      for (INT d_delta = 0; d_delta < dMax; d_delta++)
      {
        for (INT n_zalpha = 0; n_zalpha < n_zGMax; n_zalpha++)
        {
          for (INT n_zdelta = 0; n_zdelta < n_zalpha + 1; n_zdelta++) //One uses the symmetry delta <-> alpha.
          {
            Izx(x, d_alpha + 2 * d_delta, n_zalpha, n_zdelta) = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax, DMax);
            Izx(x, d_delta + 2 * d_alpha, n_zdelta, n_zalpha) = arma::zeros(basis.n_zGlobalMax, basis.n_zGlobalMax, DMax);
            arma::cube &IzCubeE = Izx(x, d_alpha + 2 * d_delta, n_zalpha, n_zdelta);
            arma::cube &IzCubeE_sym = Izx(x, d_delta + 2 * d_alpha, n_zdelta, n_zalpha);

            for (INT d_gamma = 0; d_gamma < dMax; d_gamma++)
            {
              for (INT d_beta = 0; d_beta < dMax; d_beta++)
              {
                INT Dd     = 2 * d_beta + d_delta;
                INT De     = 2 * d_beta + d_gamma;
                INT De_sym = 2 * d_gamma + d_beta;

                for (INT n_zgamma = 0; n_zgamma < n_zGMax; n_zgamma++)
                {
                  for (INT n_zbeta = 0; n_zbeta < n_zGMax; n_zbeta++)
                  {
                    IzCubeE(n_zgamma, n_zbeta, De) = IzDirect(x, d_alpha, d_gamma, n_zalpha, n_zgamma)(n_zdelta, n_zbeta, Dd);
                    IzCubeE_sym(n_zbeta, n_zgamma, De_sym) = IzCubeE(n_zgamma, n_zbeta, De);
                  }//n_zbeta
                }//n_zgamma
              }//d_beta
            }//d_gamma
          }//n_zdelta
        }//n_zalpha
      }//d_delta
    }//d_alpha
  }//x

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldCoulombExact::calcDirectField(void)
{
  DBG_ENTER;

  calculatingLength = -1.0;

  if (!field(NEUTRON, DIRECT).empty() && !field(PROTON, DIRECT).empty()) DBG_LEAVE;

  double startTime = Tools::clock();

  //Dependencies
  calcJcoul();

  //Usefull quantities
  INT N_z = 2 * basis.n_zGlobalMax; //maximum for n_za and n_zb.
  INT dMax = basis.dMax;
  INT DMax = (dMax - 1) * 3 + 1;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Rmat.
  FMulti<arma::mat> Rmat({MPRANGE}); //Rmat(mp)(vectorized(n_delta,n_zdelta,d_delta),vectorized(n_delta,n_zdelta,d_delta))

  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  HOqn.calcBlocks({0, 4});

  for (UINT i = 0; i < nbBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m = bras(0, 0);
    INT s = bras(4, 0);

    for (UINT j = 0; j < nbBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT mp = kets(0, 0);
      INT sp = kets(4, 0);

      if (mp == m && s == sp)
      {
        if (s == 0)
        {
          Rmat(mp) =  rho(PROTON ).submat(bras.filter, kets.filter);
        }//if

        if (s == 1)
        {
          Rmat(mp) += rho(PROTON ).submat(bras.filter, kets.filter);
        }//if
      }//if
    }//j
  }//i

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges Rmat with Jcoul.
  arma::cube Rj = arma::zeros(2 * basis.n_zGlobalMax, basis.Nmaxr, DMax); //Rj(n_za,n_a,d_alpha+2*d_gamma)

  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < d_alpha + 1 ; d_gamma++) //Symmetry d_alpha <=> d_gamma
    {
      for (INT n_za = 0 ; n_za < N_z ; n_za++)
      {
        for (INT n_a = 0 ; n_a < basis.Nmaxr ; n_a++)
        {
          double sum_m = 0;

          for (INT mp = 0 ; mp < basis.mMax ; mp++)
          {
            arma::cube &jcoul = Jcoul(n_za, d_gamma + 2 * d_alpha, mp);
            sum_m += arma::accu(jcoul.slice(n_a) % Rmat(mp));
          }//mp

          Rj(n_za, n_a, d_alpha + 2 * d_gamma) = sum_m;
          Rj(n_za, n_a, d_gamma + 2 * d_alpha) = sum_m;
        }//n_a
      }//n_za
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges Rj with Talmanr.

  //Test(m,n_alpha,n_gamma)(n_za,d_alpha+2*d_gamma)
  FMulti<arma::mat> TRj({MPRANGE, NRANGE, NRANGE});

  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //Symmetry n_alpha <=> n_gamma
      {
        arma::vec &talmanr = basis.talmanr(m, n_alpha, m, n_gamma);
        arma::mat &trj     = TRj(m, n_alpha, n_gamma);
        trj = arma::zeros(N_z, DMax);

        for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
          {
            trj.col(d_alpha + 2 * d_gamma) = Rj.slice(d_alpha + 2 * d_gamma) * talmanr;
          }//d_gamma
        }//d_alpha

        TRj(m, n_gamma, n_alpha) = trj;
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges TRj with Talmanz.

  //Gamma(m,n_alpha,n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma)
  FMulti<arma::cube> Gamma({MPRANGE, NRANGE, NRANGE}) ;

  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //Symmetry alpha <=> gamma
      {
        arma::cube &gamma = Gamma(m, n_alpha, n_gamma);
        arma::cube &gamma_sym = Gamma(m, n_gamma, n_alpha);
        gamma = arma::zeros(basis.n_zMax(m, n_alpha), basis.n_zMax(m, n_gamma), DMax);
        gamma_sym = arma::zeros(basis.n_zMax(m, n_gamma), basis.n_zMax(m, n_alpha), DMax);
        arma::mat  &trj   = TRj(m, n_alpha, n_gamma);

        for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
        {
          for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
          {
            for (INT n_zalpha = 0 ; n_zalpha < basis.n_zMax(m, n_alpha) ; n_zalpha++)
            {
              for (INT n_zgamma = 0 ; n_zgamma < basis.n_zMax(m, n_gamma) ; n_zgamma++)
              {
                arma::vec talvecz = basis.talmanz(n_zalpha, d_alpha, n_zgamma, d_gamma).subvec(0, N_z - 1);
                gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::accu(trj.col(d_alpha + 2 * d_gamma) % talvecz);
              }//n_zgamma
            }//n_zalpha

            gamma_sym.slice(d_gamma + 2 * d_alpha) = (gamma.slice(d_alpha + 2 * d_gamma)).t();
          }//d_gamma
        }//d_alpha
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds the interaction.
  field(PROTON, DIRECT) = arma::zeros(HOqn.nb, HOqn.nb);
  field(NEUTRON, DIRECT) = arma::zeros(HOqn.nb, HOqn.nb);

  UINT nbMBlocks_2 =  HOqn.calcBlocks({0, 1, 3, 4});

  for (UINT i = 0; i < nbMBlocks_2; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT n_alpha =   bras(1, 0);
    INT d_alpha =   bras(3, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks_2; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT n_gamma = kets(1, 0);
      INT d_gamma = kets(3, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma && s_alpha == s_gamma)//Mean Direct part.
      {
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
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

/** Calculates the coulomb integral via the Gauss-Legendre quadrature.
 */

void FieldCoulombExact::calcIcoul(void)
{
  DBG_ENTER;

  if (!Icoul.empty()) DBG_LEAVE;


  //Dependencies
  basis.calcTalmanz();
  basis.calcMoshinskyz();
  basis.calcTalmanr();
  basis.calcMoshinskyr();

  //Usefull quantities
  double gamma = std::pow(basis.b_r / basis.b_z, 2.0) - 1.0;
  double factor_b = basis.b_r / basis.b_z;
  INT N_d  = 2 * basis.Nmaxr; //maximum for n_d.
  INT N_zd = 4 * basis.n_zGlobalMax; //maximum for n_zd.
  arma::vec zerosVec = arma::zeros(nGLE);

  Icoul = FMulti<arma::vec>({DDRANGE, DDRANGE, {0, N_d}});

  //Creates quadrature quantities for GaussLegendre [-1,1].
  //Since the integrand is an even function, [0,1] is directly related to [-1,1].
  Axis axis_gl(Axis::GAUSS_LEGENDRE, nGLE, -1.0, 1.0);
  arma::vec pre_den = arma::pow(axis_gl.p, 2) * gamma + arma::ones(nGLE);

  //Initializes Icoul
  for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
      {
        for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
        {
          for (INT n_d = 0 ; n_d < N_d ; n_d++)
          {
            Icoul(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_d) = arma::zeros(N_zd);
          }//n_d
        }//d_beta
      }//d_delta
    }//d_gamma
  }//d_alpha

  for (INT n_d = 0 ; n_d < N_d ; n_d++)
  {
    arma::vec &icoul = Icoul(0, 0, n_d);
    double factor_global = 0.5 * std::pow(factor_b, 2 * n_d) / basis.b_z;
    arma::vec den = 1.0 / arma::pow(pre_den, n_d + 1);

    for (INT n_zd = 0 ; n_zd < N_zd ; n_zd++)
    {
      arma::vec vecPhi = factor_global * basis.zPartNorm(zerosVec, n_zd);
      arma::vec num = arma::pow(axis_gl.p, 2 * n_d + n_zd);
      icoul(n_zd) = arma::accu(vecPhi % num % den % axis_gl.w);
    }//n_zd
  }//n_d

  if (basis.dMax == 2) //2ct
  {
    for (INT n_d = 0 ; n_d < N_d ; n_d++)
    {
      arma::vec &icoul_1p = Icoul(1, 0, n_d);
      arma::vec &icoul_1m = Icoul(0, 1, n_d);
      arma::vec &icoul_2p = Icoul(3, 0, n_d);
      arma::vec &icoul_2m = Icoul(0, 3, n_d);

      double factor_global = 0.5 * std::pow(factor_b, 2 * n_d) / basis.b_z;
      arma::vec den = 1.0 / arma::pow(pre_den, n_d + 1);

      for (INT n_zd = 0 ; n_zd < N_zd ; n_zd++)
      {
        arma::vec num = arma::pow(axis_gl.p, 2 * n_d + n_zd);

        //if d_tot = d_0/(2*sqrt(2)) or - d_0/(2*sqrt(2))
        double d_tot_1 = basis.d_0 / (2 * sqrt(2));
        arma::vec vecPhi_1 = factor_global * basis.zPartNorm(d_tot_1 * axis_gl.p / basis.b_z, n_zd);
        arma::vec fac_exp_1 = exp(-0.5 * arma::pow(d_tot_1 / basis.b_z * axis_gl.p, 2 ) );

        icoul_1p(n_zd) = arma::accu(vecPhi_1 % num % den % axis_gl.w % fac_exp_1);
        icoul_1m(n_zd) = std::pow(-1.0, n_zd) * icoul_1p(n_zd);

        //if d_tot = d_0/(sqrt(2)) or -d_0/(sqrt(2))
        double d_tot_2 = basis.d_0 / (sqrt(2));
        arma::vec vecPhi_2 = factor_global * basis.zPartNorm(d_tot_2 * axis_gl.p / basis.b_z, n_zd);
        arma::vec fac_exp_2 = exp(-0.5 * arma::pow(d_tot_2 / basis.b_z * axis_gl.p, 2 ) );

        icoul_2p(n_zd) = arma::accu(vecPhi_2 % num % den % axis_gl.w % fac_exp_2);
        icoul_2m(n_zd) = std::pow(-1.0, n_zd) * icoul_2p(n_zd);
      }//n_zd

      Icoul(1, 1, n_d) = Icoul(0, 0, n_d);
      Icoul(2, 2, n_d) = Icoul(0, 0, n_d);
      Icoul(3, 3, n_d) = Icoul(0, 0, n_d);
      Icoul(1, 2, n_d) = Icoul(0, 0, n_d);
      Icoul(2, 1, n_d) = Icoul(0, 0, n_d);

      Icoul(2, 0, n_d) = Icoul(1, 0, n_d);
      Icoul(3, 1, n_d) = Icoul(1, 0, n_d);
      Icoul(3, 2, n_d) = Icoul(1, 0, n_d);

      Icoul(0, 2, n_d) = Icoul(0, 1, n_d);
      Icoul(1, 3, n_d) = Icoul(0, 1, n_d);
      Icoul(2, 3, n_d) = Icoul(0, 1, n_d);
    }//n_d
  }//if 2ct

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Jcoulomb.
 */

void FieldCoulombExact::calcJcoul(void)
{
  DBG_ENTER;

  if (!Jcoul.empty()) DBG_LEAVE;

  //Dependencies
  calcIcoul();

  double factor_global = sqrt(2. / sqrt(PI)) * 2.0 * HBARC * ALPHA;

  if (general.compatibility == General::COMPAT_BERGER) factor_global = sqrt(2. / sqrt(PI)) * 2.88;

  INT N_z = 2 * basis.n_zGlobalMax; //maximum for n_za and n_zb.
  INT n_zMax = basis.n_zGlobalMax;
  INT dMax = basis.dMax;
  INT N_d  = 2 * basis.Nmaxr;
  IVEC Jdim = arma::zeros<IVEC>(basis.mMax);
  FMulti<ICUBE> Jindex({MPRANGE});

  Jcoul = FMulti<arma::cube>({{0, N_z}, DDRANGE, MPRANGE});

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges Icoul with Moshinskyz and Talmanz.
  Multi<arma::mat> Jz; //Jz(d_alpha+2*d_gamma,d_delta+2*d_beta,n_zdelta,n_zbeta)(n_za,n_d)

  //Initializes Jz
  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          for (INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < n_zMax ; n_zbeta++)
            {
              Jz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta) = arma::zeros(N_z, N_d);
            }//n_zbeta
          }//n_zdelta
        }//d_beta
      }//d_delta
    }//d_gamma
  }//d_alpha

  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          for (INT n_d = 0 ; n_d < N_d ; n_d++)
          {
            //Sets corrected moshinskyz
            arma::mat mosz = factor_global * basis.moshinskyz.submat(0, 0, N_z - 1, N_z - 1);
            arma::vec &icoul = Icoul(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_d);

            for (INT n_za = 0 ; n_za < N_z ; n_za++)
            {
              for (INT n_zb = 0 ; n_zb < N_z ; n_zb++)
              {
                mosz(n_za, n_zb) = mosz(n_za, n_zb) * icoul(n_za + n_zb);
              }//n_zb
            }//n_za

            for (INT n_zdelta = 0 ; n_zdelta < n_zMax ; n_zdelta++)
            {
              for (INT n_zbeta = 0 ; n_zbeta < n_zdelta + 1 ; n_zbeta++) //symmetry delta <=> beta
              {
                Jz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta).col(n_d) = mosz * basis.talmanz(n_zdelta, d_delta, n_zbeta, d_beta).subvec(0, N_z - 1);
                Jz(d_alpha + 2 * d_gamma, d_beta + 2 * d_delta, n_zbeta, n_zdelta).col(n_d) = Jz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta).col(n_d);
              }//n_zbeta
            }//n_zdelta
          }//n_d
        }//d_beta
      }//d_delta
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges Jz with Moshinskyr.

  Multi<arma::cube> MJz; //Jzr(d_alpha+2*d_gamma,d_delta+2*d_beta,n_zdelta,n_zbeta) (n_a,n_b,n_za)

  for (INT d_delta = 0 ; d_delta < basis.dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < basis.dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
        {
          for (INT n_zdelta = 0 ; n_zdelta < basis.n_zGlobalMax ; n_zdelta++)
          {
            for (INT n_zbeta = 0 ; n_zbeta < n_zdelta + 1 ; n_zbeta++) //symmetry delta <=> beta
            {
              arma::cube &mjz = MJz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta);
              mjz = arma::zeros(basis.Nmaxr, basis.Nmaxr, N_z);
              arma::mat &jz = Jz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta);

              for (INT n_za = 0 ; n_za < 2 * basis.n_zGlobalMax ; n_za++)
              {
                for (INT n_a = 0 ; n_a < basis.Nmaxr ; n_a++)
                {
                  for (INT n_b = 0 ; n_b < n_a + 1 ; n_b++) //symmetry n_a <=> n_b
                  {
                    mjz(n_a, n_b, n_za) = basis.moshinskyr(0, n_a)(n_b) * jz(n_za, n_a + n_b);
                    mjz(n_b, n_a, n_za) = mjz(n_a, n_b, n_za);
                  }//n_b
                }//n_a

                MJz(d_alpha + 2 * d_gamma, d_beta + 2 * d_delta, n_zbeta, n_zdelta) = mjz;
              }//n_za
            }//n_zbeta
          }//n_zdelta
        }//d_gamma
      }//d_alpha
    }//d_beta
  }//d_delta


  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges MJz with Talmanr.
  //Stores J's dimensions
  for (INT mp = 0 ; mp < basis.mMax ; mp++)
  {
    INT i = 0;

    for (INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      for (INT n_z = 0 ; n_z < basis.n_zMax(mp, n) ; n_z++)
      {
        for (INT d = 0 ; d < dMax ; d++)
        {
          i += 1;
        }//d
      }//n_z
    }//n

    Jdim(mp) = i;
  }//mp

  //Stores J's index
  for (INT mp = 0 ; mp < basis.mMax ; mp++)
  {
    INT i = 0;
    ICUBE &jindex = Jindex(mp);
    jindex = arma::zeros<ICUBE>(basis.nMax(mp), basis.n_zMax(mp, 0), dMax);

    for (INT n = 0 ; n < basis.nMax(mp) ; n++)
    {
      for (INT n_z = 0 ; n_z < basis.n_zMax(mp, n) ; n_z++)
      {
        for (INT d = 0 ; d < dMax ; d++)
        {
          jindex(n, n_z, d) = i;
          i += 1;
        }//d
      }//n_z
    }//n
  }//mp

  for (INT mp = 0 ; mp < basis.mMax ; mp++)
  {
    ICUBE &jindex = Jindex(mp);

    for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < d_alpha + 1 ; d_gamma++) //Symmetry d_alpha <=> d_gamma
      {
        for (INT n_za = 0 ; n_za < N_z ; n_za++)
        {
          arma::cube &jcube = Jcoul(n_za, d_alpha + 2 * d_gamma, mp);
          jcube = arma::zeros(Jdim(mp), Jdim(mp), basis.Nmaxr);

          for (INT n_delta = 0 ; n_delta < basis.nMax(mp) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < n_delta + 1 ; n_beta++) //Symmetry delta <=> beta
            {
              arma::vec talvec = basis.talmanr(mp, n_delta, mp, n_beta);

              for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
              {
                for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
                {
                  for (INT n_zdelta = 0 ; n_zdelta < basis.n_zMax(mp, n_delta) ; n_zdelta++)
                  {
                    INT index_delta = jindex(n_delta, n_zdelta, d_delta);

                    for (INT n_zbeta = 0 ; n_zbeta < basis.n_zMax(mp, n_beta) ; n_zbeta++)
                    {
                      INT index_beta = jindex(n_beta, n_zbeta, d_beta);

                      arma::cube &mjz = MJz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta, n_zdelta, n_zbeta);
                      jcube.tube(index_delta, index_beta) = mjz.slice(n_za) * talvec;
                      jcube.tube(index_beta, index_delta) = jcube.tube(index_delta, index_beta);

                    }//n_zbeta
                  }//n_zdelta
                }//d_beta
              }//d_delta
            }//n_bet
          }//n_delta

          Jcoul(n_za, d_gamma + 2 * d_alpha, mp) = jcube;
        }//n_za
      }//d_gamma
    }//d_alpha
  }//mp

  DBG_LEAVE;
}


