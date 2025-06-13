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

#include "field_density.h"
#include "tools.h"
#include "basis.h"
#include "interaction.h"
#include "plot.h"
#include "quadratures.h"

/** \file
 *  \brief Methods of the FieldDensity class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldDensity::FieldDensity(Field::Parameters fp, State *_state) : Field(fp, _state),
  discrete(&basis),
  discrete0(&basis),
  discrete1(&basis),
  discrete2(&basis)

{
  DBG_ENTER;
  name = "density";
  shortName = "Densi.";
  nGLA = 40; //26 Berger / Initially 40
  nGHE = 100; //54 Berger / Initially 100

  discrete.mesh = Mesh::gaussLaguerreHermite(nGLA, nGHE);
  discrete0 = discrete;
  discrete0.mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
  discrete0.mesh.az.p = basis.b_z * discrete.mesh.az.p;

  if (basis.dMax == 2) //2ct
  {
    discrete1 = discrete0;
    discrete2 = discrete0;
    discrete1.mesh.az.p +=  0.5 * basis.d_0;
    discrete2.mesh.az.p += -0.5 * basis.d_0;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensity::calcField(void)
{
  DBG_ENTER;

  calculatingLength = -1.0;

  if (!field(NEUTRON, DIRECT).empty() && !field(PROTON, DIRECT).empty()) DBG_LEAVE;

  double startTime = Tools::clock();

  //dependencies
  calcPreZ();
  calcPreR();

  //usefull quantities
  double fac_global = PI * basis.b_z * std::pow(basis.b_r, 2);
  INT DMax = (basis.dMax - 1) * 3 + 1;

  double t = parameters["t"];
  double x = parameters["x"];
  double a = parameters["a"];

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Rmat.
  Multi<arma::mat> localRho; //localRho(isospin,dcase)
  Multi<arma::mat> Rmat; //Rmat(isospin,dcase)

  localRho(NEUTRON, 0) = discrete0.getLocalXZ(rho(NEUTRON), true);
  localRho(PROTON, 0) = discrete0.getLocalXZ(rho(PROTON ), true);
  localRho(TOTAL, 0)  = localRho(NEUTRON, 0) + localRho(PROTON, 0);
  localRho(3, 0)       = arma::pow(localRho(TOTAL, 0), a);
  Rmat(NEUTRON, 0) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 0) - (x + 0.5) * localRho(NEUTRON, 0)) % localRho(3, 0);
  Rmat(PROTON, 0) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 0) - (x + 0.5) * localRho(PROTON, 0)) % localRho(3, 0);

  if (basis.dMax == 2) //2ct
  {
    localRho(NEUTRON, 1) = discrete1.getLocalXZ(rho(NEUTRON), true);
    localRho(PROTON, 1) = discrete1.getLocalXZ(rho(PROTON ), true);
    localRho(TOTAL, 1)  = localRho(NEUTRON, 1) + localRho(PROTON, 1);
    localRho(3, 1)       = arma::pow(localRho(TOTAL, 1), a);
    Rmat(NEUTRON, 1) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 1) - (x + 0.5) * localRho(NEUTRON, 1)) % localRho(3, 1);
    Rmat(PROTON, 1) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 1) - (x + 0.5) * localRho(PROTON, 1)) % localRho(3, 1);


    localRho(NEUTRON, 2) = discrete2.getLocalXZ(rho(NEUTRON), true);
    localRho(PROTON, 2) = discrete2.getLocalXZ(rho(PROTON ), true);
    localRho(TOTAL, 2)  = localRho(NEUTRON, 2) + localRho(PROTON, 2);
    localRho(3, 2)       = arma::pow(localRho(TOTAL, 2), a);
    Rmat(NEUTRON, 2) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 2) - (x + 0.5) * localRho(NEUTRON, 2)) % localRho(3, 2);
    Rmat(PROTON, 2) = t * ((1.0 + x / 2.0) * localRho(TOTAL, 2) - (x + 0.5) * localRho(PROTON, 2)) % localRho(3, 2);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with zPart.
  Multi<arma::vec> Rz; //Rz(iso,n_zalpha,n_gamma,d_alpha+2*d_gamma)(r)

  if (basis.dMax == 1 ) //1ct
  {
    for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
    {
      for (INT n_zgamma = 0 ; n_zgamma < n_zalpha + 1 ; n_zgamma++) //Symmetry
      {
        arma::vec &prez = preZ(n_zalpha, n_zgamma, 0);

        for (INT iso: {NEUTRON, PROTON})
        {
          Rz(iso, n_zalpha, n_zgamma, 0) = Rmat(iso, 0) * prez;
          Rz(iso, n_zgamma, n_zalpha, 0) = Rmat(iso, 0) * prez;
        }//iso
      }//n_zgamma
    }//n_zalpha
  }//if

  if (basis.dMax == 2) //2ct
  {
    for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
    {
      for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
      {
        INT dcase = 0;

        if (d_alpha == 0 && d_gamma == 0)
        {
          dcase = 1;
        }

        if (d_alpha == 1 && d_gamma == 1)
        {
          dcase = 2;
        }

        for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
        {
          for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
          {
            arma::vec &prez = preZ(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);

            for (INT iso: {NEUTRON, PROTON})
            {
              Rz(iso, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = Rmat(iso, dcase) * prez;
            }//iso
          }//n_zgamma
        }//n_zalpha
      }//d_gamma
    }//d_alpha
  }//if

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with rPart.
  Multi<arma::cube> Gamma; //Gamma(iso,m,n_alpha,n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma)

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT m = 0 ; m < basis.mMax ; m++)
    {
      for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //Symmetry
        {
          arma::vec &prer = preR(m, n_alpha, n_gamma);
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
                  gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_global * arma::accu(prer % Rz(iso, n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                  gamma_sym(n_zgamma, n_zalpha, d_gamma + 2 * d_alpha) = gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
                }//n_zgamma
              }//n_zalpha
            }//d_gamma
          }//d_alpha
        }//n_gamma
      }//n_alpha
    }//m
  }//iso

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Gamma.
  field(NEUTRON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  Qnumbers &HOqn = basis.HOqn;
  UINT nbBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});

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
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean Exchange part --.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(NEUTRON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
        field(PROTON, DIRECT).submat(bras.filter, kets.filter) = Gamma(PROTON, m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
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

void FieldDensity::calcPreZ(void)
{
  DBG_ENTER;

  if (!preZ.empty()) DBG_LEAVE;

  //usefull quantities;
  double fac_bz  = std::pow(1.0 / sqrt(basis.b_z), 2);
  double fac_exp = 1.0;

  for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
  {
    double dv_alpha = (d_alpha - 0.5) * basis.d_0;

    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      double dv_gamma = (d_gamma - 0.5) * basis.d_0;
      double k_ag = (dv_alpha - dv_gamma) / std::sqrt(2);
      fac_exp = exp(-0.5 * std::pow(k_ag / basis.b_z, 2));

      arma::vec nodes_alpha = discrete.mesh.az.p + (dv_alpha - dv_gamma) / (2 * basis.b_z);
      arma::vec nodes_gamma = discrete.mesh.az.p + (dv_gamma - dv_alpha) / (2 * basis.b_z);

      for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
        {
          preZ(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_exp * fac_bz * basis.zPartNormReduced(nodes_alpha, n_zalpha) % basis.zPartNormReduced(nodes_gamma, n_zgamma) % discrete.mesh.az.w;
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldDensity::calcPreR(void)
{
  DBG_ENTER;

  if (!preR.empty()) DBG_LEAVE;

  //usefull quantities
  double fac_br = 1.0 / (basis.b_r * sqrt(PI));
  fac_br = fac_br * fac_br;

  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    arma::vec nodes_m = arma::pow(discrete.mesh.ax.p, m);

    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      double fac_alpha = sqrt(fact[n_alpha] / fact[n_alpha + m]);

      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        double fac_gamma = sqrt(fact[n_gamma] / fact[n_gamma + m]);
        preR(m, n_alpha, n_gamma) = fac_br * fac_alpha * fac_gamma * basis.laguerre(m, n_alpha, discrete.mesh.ax.p) % basis.laguerre(m, n_gamma, discrete.mesh.ax.p) % nodes_m % discrete.mesh.ax.w;
      }//n_gamma
    }//n_alpha
  }//m

  DBG_LEAVE;
}

