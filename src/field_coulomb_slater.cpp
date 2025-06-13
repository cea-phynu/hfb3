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

#include "field_coulomb_slater.h"
#include "axis.h"
#include "interaction.h"
#include "global.h"
#include "multi.h"
#include "tools.h"
#include "mesh.h"
#include "state.h"

/** \file
 *  \brief Methods of the FieldCoulombSlater class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldCoulombSlater::FieldCoulombSlater(Field::Parameters fp, State *_state) :
  Field(fp, _state),
  discrete(&(state.basis)),
  discrete0(&(state.basis))
{
  DBG_ENTER;

  name      = "coulomb (Slater)";
  shortName = "CoulSl";

  nGLE      = 40;
  nGLA      = 40;
  nGHE      = 150;

  discrete.mesh = Mesh::gaussLaguerreHermite(nGLA, nGHE);
  discrete0 = discrete;
  discrete0.mesh.ax.p = arma::sqrt(discrete.mesh.ax.p);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldCoulombSlater::calcField(void)
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
  Multi<arma::mat> Rmat; //Rmat(mp)(vectorized(n_delta,n_zdelta,d_delta),vectorized(n_delta,n_zdelta,d_delta))
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
  Multi<arma::mat> TRj ; //Test(m,n_alpha,n_gamma)(n_za,d_alpha+2*d_gamma)

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
  Multi<arma::cube> Gamma ; //Gamma(m,n_alpha,n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma)

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

void FieldCoulombSlater::calcIcoul(void)
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

void FieldCoulombSlater::calcJcoul(void)
{
  DBG_ENTER;

  if (!Jcoul.empty()) DBG_LEAVE;

  //Dependencies
  calcIcoul();

  double factor_global = sqrt(2. / sqrt(PI)) * 2.0 * HBARC * ALPHA;

  if (general.compatibility == General::COMPAT_HFBTHO) factor_global = sqrt(2. / sqrt(PI)) * 2.0 * 1.4399784085965135;
  if (general.compatibility == General::COMPAT_BERGER) factor_global = sqrt(2. / sqrt(PI)) * 2.88;

  INT N_z = 2 * basis.n_zGlobalMax; //maximum for n_za and n_zb.
  INT n_zMax = basis.n_zGlobalMax;
  INT dMax = basis.dMax;
  INT N_d  = 2 * basis.Nmaxr;
  IVEC Jdim = arma::zeros<IVEC>(basis.mMax);
  Multi<ICUBE> Jindex;

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

//==============================================================================
//==============================================================================
//==============================================================================
/** Calculate the Energy contribution.
 */

void FieldCoulombSlater::calcEnergy(void)
{
  DBG_ENTER;

  // calculate non-specialized energies
  Field::calcEnergy();

  double fac_slater = -(3.0 / 4.0) * HBARC * ALPHA * std::pow(3.0 / PI, 1.0 / 3.0); // HBARC*ALPHA = e^2

  if (general.compatibility == General::COMPAT_HFBTHO) fac_slater = - 0.75 * 1.4399784085965135 * std::pow(3.0 / PI, 1.0 / 3.0);

  double fac_integration  = PI;
  arma::mat Rmat = arma::pow(discrete0.getLocalXZ(rho(PROTON ), true), 4. / 3);
  arma::vec Rz = Rmat * discrete.mesh.az.we;
  energy(PROTON, EXCHANGE) = fac_slater * fac_integration * arma::accu(Rz % discrete.mesh.ax.we);
  energy(NEUTRON, EXCHANGE) = 0.0;

  DBG_LEAVE;
}


