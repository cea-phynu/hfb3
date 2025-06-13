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

#include "field_rearrangement.h"
#include "tools.h"
#include "state.h"
#include "discrete.h"
#include "interaction.h"

/** \file
 *  \brief Methods of the FieldRearrangement class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldRearrangement::FieldRearrangement(Field::Parameters fp, State *_state) :
  Field(fp, _state),
  discrete(&basis),
  discrete0(&basis)
{
  DBG_ENTER;

  name = "rearrangement";
  shortName = "Rearr.";

  nGLA = 20;//20
  nGHE = 70;//

  contributeToEnergy = false;

  discrete.mesh = Mesh::gaussLaguerreHermite(nGLA, nGHE);
  discrete0 = discrete;
  discrete0.mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
  discrete0.mesh.az.p = basis.b_z * discrete.mesh.az.p;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangement::calcField(void)
{
  DBG_ENTER;

  //Useful tools.
  double t = parameters["t"];
  double x = parameters["x"];
  double a = parameters["a"];

  if (!field(NEUTRON, DIRECT).empty() && !field(PROTON, DIRECT).empty()) DBG_LEAVE;

  calculatingLength = -1.0;
  double startTime = Tools::clock();

  //dependencies
  calcPreZ();
  calcPreR();

  //usefull quantities
  double fac_global = 2 * PI * basis.b_z * std::pow(basis.b_r, 2);
  INT DMax = (basis.dMax - 1) * 3 + 1;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Rmat.
  Multi<arma::mat> localRho; //localRho(isospin)
  arma::mat Rmat; //(r,z)

  localRho(NEUTRON) = discrete0.getLocalXZ(rho(NEUTRON), true);
  localRho(PROTON ) = discrete0.getLocalXZ(rho(PROTON ), true);
  localRho(TOTAL  ) = localRho(NEUTRON) + localRho(PROTON);
  localRho(3)       = arma::pow(localRho(TOTAL), a - 1);
  localRho(4)       = arma::pow(localRho(TOTAL), 2);
  localRho(5)       = arma::pow(localRho(NEUTRON), 2) + arma::pow(localRho(PROTON), 2);
  Rmat = 0.25 * a * t * ((1.0 + x / 2.0) * localRho(4) - (x + 0.5) * localRho(5)) % localRho(3);

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with zPart.
  Multi<arma::vec> Rz; //Rz(n_zalpha,n_gamma,d_alpha+2*d_gamma)(r)

  for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
        {
          arma::vec &prez = preZ(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
          Rz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = Rmat * prez;
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  //////////////////////////////////////////////////////////////////////////////////////////////////////Merges with rPart.
  Multi<arma::cube> Gamma; //Gamma(m,n_alpha,n_gamma)(n_zalpha,n_zgamma,d_alpha+2*d_gamma)

  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < n_alpha + 1 ; n_gamma++) //Symmetry
      {
        arma::vec &prer = preR(m, n_alpha, n_gamma);
        arma::cube &gamma = Gamma(m, n_alpha, n_gamma);
        arma::cube &gamma_sym = Gamma(m, n_gamma, n_alpha);
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
                gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_global * arma::accu(prer % Rz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma));
                gamma_sym(n_zgamma, n_zalpha, d_gamma + 2 * d_alpha) = gamma(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma);
              }//n_zgamma
            }//n_zalpha
          }//d_gamma
        }//d_alpha
      }//n_gamma
    }//n_alpha
  }//m

  //////////////////////////////////////////////////////////////////////////////////////////////////////Builds Gamma.
  field(NEUTRON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

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
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1) //Mean Exchange part --.
      {
        field(NEUTRON, DIRECT).submat(bras.filter, kets.filter) = Gamma(m_alpha, n_alpha, n_gamma).slice(d_alpha + 2 * d_gamma);
      }//if
    }//j
  }//i

  field(PROTON, DIRECT) = field(NEUTRON, DIRECT);

  //Stores the calculating length
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangement::calcPreZ(void)
{
  DBG_ENTER;

  if (!preZ.empty()) DBG_LEAVE;

  //usefull quantities;
  double fac_bz  = std::pow(1.0 / sqrt(basis.b_z), 2);

  for (INT d_alpha = 0 ; d_alpha < basis.dMax ; d_alpha++)
  {
    double dv_alpha = (d_alpha - 0.5) * basis.d_0;

    for (INT d_gamma = 0 ; d_gamma < basis.dMax ; d_gamma++)
    {
      double dv_gamma = (d_gamma - 0.5) * basis.d_0;

      arma::vec nodes_alpha = discrete.mesh.az.p + dv_alpha / basis.b_z;
      arma::vec nodes_gamma = discrete.mesh.az.p + dv_gamma / basis.b_z;

      for (INT n_zalpha = 0 ; n_zalpha < basis.n_zGlobalMax ; n_zalpha++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < basis.n_zGlobalMax ; n_zgamma++)
        {
          preZ(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = fac_bz * basis.zPartNorm(nodes_alpha, n_zalpha) % basis.zPartNorm(nodes_gamma, n_zgamma) % discrete.mesh.az.we;
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldRearrangement::calcPreR(void)
{
  DBG_ENTER;

  if (!preR.empty()) DBG_LEAVE;

  //usefull quantities
  double fac_br = 1.0 / (basis.b_r * sqrt(PI));
  fac_br = fac_br * fac_br;

  for (INT m = 0 ; m < basis.mMax ; m++)
  {
    for (INT n_alpha = 0 ; n_alpha < basis.nMax(m) ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < basis.nMax(m) ; n_gamma++)
      {
        preR(m, n_alpha, n_gamma) = fac_br * basis.rPartNorm(discrete.mesh.ax.p, m, n_alpha) % basis.rPartNorm(discrete.mesh.ax.p, m, n_gamma) % discrete.mesh.ax.we;
      }//n_gamma
    }//n_alpha
  }//m

  DBG_LEAVE;
}

