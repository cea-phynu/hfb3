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

#include "field_cdm2.h"
#include "global.h"
#include "tools.h"
#include "state.h"
#include "interaction.h"

/** \file
 *  \brief Methods of the FieldCDM2 class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldCDM2::FieldCDM2(Field::Parameters _parameters, State *_state) : Field(_parameters, _state)
{
  DBG_ENTER;

  name = "2-body COM cor.";
  shortName = "COMas2";

  if (general.compatibility == General::COMPAT_BERGER)
  {
    Tools::warning("Removing CDM2 pairing contribution ! (general.compatibility with BERGER2CT solver)");
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the 2-body center of mass field.
 */

void FieldCDM2::calcField(void)
{
  DBG_ENTER;

  // No CDM2 field in HFBTHO;
  if (general.compatibility == General::COMPAT_HFBTHO)
  {
    Tools::warning("deactivating CDM2 field (HFBTHO general.compatibility)");
    DBG_LEAVE;
  }

  if (!field(NEUTRON, DIRECT ).empty() && !field(PROTON, DIRECT ).empty() &&
      !field(NEUTRON, PAIRING).empty() && !field(PROTON, PAIRING).empty()) DBG_LEAVE;

  calculatingLength = -1.0;
  double startTime = Tools::clock();

  //dependencies
  calcPblocks();

  field(NEUTRON, PAIRING) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, PAIRING) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  field(NEUTRON, DIRECT) =  arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, DIRECT) =  arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates block 0.
  for (INT omega = 0 ; omega < basis.mMax ; omega++)
  {
    arma::mat P0_omega = P(0).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
    field(NEUTRON, DIRECT).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += - P0_omega * rho(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * P0_omega;
    field(PROTON, DIRECT).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += - P0_omega * rho(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * P0_omega;
    field(NEUTRON, PAIRING).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += - P0_omega * kappa(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * P0_omega;
    field(PROTON, PAIRING).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += - P0_omega * kappa(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * P0_omega;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates block time even/odd.
  INT omega = 0;
  arma::mat p_omega = p.submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
  field(NEUTRON, DIRECT).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += p_omega * rho(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * p_omega;
  field(PROTON, DIRECT).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += p_omega * rho(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * p_omega;
  field(NEUTRON, PAIRING).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += p_omega * kappa(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * p_omega;
  field(PROTON, PAIRING).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += p_omega * kappa(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * p_omega;

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates +- blocks.
  UINT nbMBlocks =  basis.HOqn.calcBlocks({0, 4});
  Multi<arma::mat> gamman_mp = SizeInit;
  Multi<arma::mat> gammap_mp = SizeInit;
  Multi<arma::mat> kappan_mp = SizeInit;
  Multi<arma::mat> kappap_mp = SizeInit;

  Qnumbers &HOqn = basis.HOqn;

  nbMBlocks =  basis.HOqn.calcBlocks({0, 4});

  for (UINT i = 0; i < nbMBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha < basis.mMax - 1 && m_gamma < basis.mMax - 1)
      {
        if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0)
        {
          gamman_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          gammap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappan_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
        }

        if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1)
        {
          gamman_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          gammap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappan_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
        }

        if (m_alpha == m_gamma + 1 && s_alpha == 1 && s_gamma == 0)
        {
          gamman_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          gammap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappan_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma).t();
          kappap_mp(m_alpha + 1, m_gamma + 1, s_alpha, s_gamma) += Pblocks(m_alpha) * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma).t();
        }
      }//if

      if (m_alpha > 1 || m_gamma > 1)
      {
        if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0)
        {
          gamman_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          gammap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappan_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        }

        if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1)
        {
          gamman_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          gammap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappan_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        }

        if (m_alpha == m_gamma + 1 && s_alpha == 1 && s_gamma == 0)
        {
          gamman_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          gammap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappan_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
          kappap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        }
      }

      if (m_alpha == 1 && m_gamma == 1 && s_alpha == s_gamma && s_alpha == 0)
      {
        gamman_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        gammap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * rho(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        kappan_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(NEUTRON)(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
        kappap_mp(m_alpha - 1, m_gamma - 1, s_alpha, s_gamma) += Pblocks(m_alpha - 1).t() * kappa(PROTON )(bras.filter, kets.filter) * Pblocks(m_gamma - 1);
      }//if
    }//j
  }//i

  //Adds the calculations to the interaction.
  for (UINT i = 0; i < nbMBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma && s_alpha == 0 && s_gamma == 0)
      {
        field(NEUTRON, DIRECT ).submat(bras.filter, kets.filter) += gamman_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, DIRECT ).submat(bras.filter, kets.filter) += gammap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) += kappan_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, PAIRING).submat(bras.filter, kets.filter) += kappap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
      }

      if (m_alpha == m_gamma && s_alpha == 1 && s_gamma == 1)
      {
        field(NEUTRON, DIRECT ).submat(bras.filter, kets.filter) += gamman_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, DIRECT ).submat(bras.filter, kets.filter) += gammap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) += kappan_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, PAIRING).submat(bras.filter, kets.filter) += kappap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
      }

      if (m_alpha == m_gamma + 1 && s_alpha == 1 && s_gamma == 0)
      {
        field(NEUTRON, DIRECT ).submat(bras.filter, kets.filter) += gamman_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, DIRECT ).submat(bras.filter, kets.filter) += gammap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) += kappan_mp(m_alpha, m_gamma, s_alpha, s_gamma);
        field(PROTON, PAIRING).submat(bras.filter, kets.filter) += kappap_mp(m_alpha, m_gamma, s_alpha, s_gamma);
      }

      if (m_alpha + 1 == m_gamma && s_alpha == 0 && s_gamma == 1) //Interaction symmetry.
      {
        field(NEUTRON, DIRECT ).submat(bras.filter, kets.filter) += gamman_mp(m_gamma, m_alpha, s_gamma, s_alpha).t();
        field(PROTON, DIRECT ).submat(bras.filter, kets.filter) += gammap_mp(m_gamma, m_alpha, s_gamma, s_alpha).t();
        field(NEUTRON, PAIRING).submat(bras.filter, kets.filter) += kappan_mp(m_gamma, m_alpha, s_gamma, s_alpha).t();
        field(PROTON, PAIRING).submat(bras.filter, kets.filter) += kappap_mp(m_gamma, m_alpha, s_gamma, s_alpha).t();
      }
    }//j
  }//i

  //Factors.
  field(NEUTRON, DIRECT ) = field(NEUTRON, DIRECT ) * HBARC2 / (NUCLEON_MASS * (double)(state.sys.nNeut + state.sys.nProt));
  field(PROTON , DIRECT ) = field(PROTON , DIRECT ) * HBARC2 / (NUCLEON_MASS * (double)(state.sys.nNeut + state.sys.nProt));
  field(NEUTRON, PAIRING) = field(NEUTRON, PAIRING) * HBARC2 / (NUCLEON_MASS * (double)(state.sys.nNeut + state.sys.nProt));
  field(PROTON , PAIRING) = field(PROTON , PAIRING) * HBARC2 / (NUCLEON_MASS * (double)(state.sys.nNeut + state.sys.nProt));

  if (general.compatibility == General::COMPAT_BERGER)
  {
    field(NEUTRON, PAIRING) *= 0.0;
    field(PROTON , PAIRING) *= 0.0;
  }

  // Store the calculating length
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates Pblocks.
 */

void FieldCDM2::calcPblocks(void)
{
  DBG_ENTER;

  if (!Pblocks.empty()) DBG_LEAVE;

  // dependencies
  calcP();

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates special blocks for P.
  Qnumbers &HOqn = basis.HOqn;
  UINT nbMBlocks =  basis.HOqn.calcBlocks({0, 4});

  for (UINT i = 0; i < nbMBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma + 1 && s_alpha == s_gamma && s_alpha == 0)
      {
        Pblocks(m_gamma) = P(1)(bras.filter, kets.filter);
      }//if
    }//j
  }//i

  //////////////////////////////////////////////////////////////////////////////////////////////////////Size initializer for blocks calculations.
  arma::mat SizeMat = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  for (UINT i = 0; i < nbMBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m_alpha =   bras(0, 0);
    INT s_alpha =   bras(4, 0);

    for (UINT j = 0; j < nbMBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT m_gamma = kets(0, 0);
      INT s_gamma = kets(4, 0);

      if (m_alpha == m_gamma)
      {
        arma::mat sizemat = SizeMat(bras.filter, kets.filter);
        SizeInit(m_alpha, m_gamma, s_alpha, s_gamma) = arma::zeros(sizemat.n_rows, sizemat.n_cols);
      }//if

      if (m_alpha == m_gamma + 1)
      {
        arma::mat sizemat = SizeMat(bras.filter, kets.filter);
        SizeInit(m_alpha, m_gamma, s_alpha, s_gamma) = arma::zeros(sizemat.n_rows, sizemat.n_cols);
        SizeInit(m_gamma, m_alpha, s_gamma, s_alpha) = arma::zeros(sizemat.n_cols, sizemat.n_rows);
      }//if
    }//j
  }//i

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates P.
 */

void FieldCDM2::calcP(void)
{
  DBG_ENTER;

  if (!P.empty() && !p.empty()) DBG_LEAVE;

  // dependencies
  calcIz();
  calcIr();

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates the 0 part.
  P(0) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  Qnumbers &HOqn = basis.HOqn;
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

      if (m_alpha == m_gamma && n_alpha == n_gamma && s_alpha == s_gamma)
      {
        P(0).submat(bras.filter, kets.filter) = (Iz_0.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_alpha, n_alpha) - 1);
      }//if
    }//j
  }//i

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates the +- part.
  P(1) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  nbMBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});

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

      if (m_alpha - 1 == m_gamma && n_alpha == n_gamma  && s_alpha == s_gamma)
      {
        P(1)(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_alpha)(n_alpha, n_gamma);
      }//if

      if (m_alpha - 1 == m_gamma && n_alpha + 1 == n_gamma  && s_alpha == s_gamma)
      {
        P(1)(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_alpha)(n_alpha, n_gamma);
      }//if
    }//j
  }//i

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates p
  p = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  nbMBlocks =  basis.HOqn.calcBlocks({0, 1, 3, 4});

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

      if (m_alpha == 0 && m_gamma == 1 && s_alpha == 0 && s_gamma == 1)
      {
        if (n_alpha == n_gamma)
        {
          p(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_gamma)(n_gamma, n_alpha);
        }

        if (n_alpha == n_gamma + 1)
        {
          p(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_gamma)(n_gamma, n_alpha);
        }
      }//if

      if (m_alpha == 1 && m_gamma == 0 && s_alpha == 1 && s_gamma == 0)
      {
        if (n_alpha == n_gamma)
        {
          p(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_alpha)(n_alpha, n_gamma);
        }

        if (n_alpha + 1 == n_gamma)
        {
          p(bras.filter, kets.filter) = (Iz_mp.slice(d_alpha + 2 * d_gamma)).submat(0, 0, basis.n_zMax(m_alpha, n_alpha) - 1, basis.n_zMax(m_gamma, n_gamma) - 1) * Ir_plus(m_alpha)(n_alpha, n_gamma);
        }
      }//if
    }//j
  }//i

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates z spatial part.
 */

void FieldCDM2::calcIz(void)
{
  DBG_ENTER;

  if (!Iz_0.empty() && !Iz_mp.empty()) DBG_LEAVE;

  // dependencies
  basis.calcTab();

  //Usefull quantities
  INT n_zGlobalMax = basis.n_zGlobalMax;
  INT dMax = basis.dMax;
  INT DMax = (dMax - 1) * 3 + 1;
  double fac_bz = 1 / (basis.b_z * sqrt(2));

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Iz_0.
  Iz_0 = arma::zeros(n_zGlobalMax, n_zGlobalMax, DMax);

  for (INT n_zbeta = 0; n_zbeta < n_zGlobalMax ; n_zbeta++)
  {
    for (INT n_zdelta = 0; n_zdelta < n_zbeta + 1 ; n_zdelta++) //symmetry  beta<=>delta x -1
    {
      for (INT d_beta = 0; d_beta < dMax ; d_beta++)
      {
        for (INT d_delta = 0; d_delta < dMax ; d_delta++)
        {
          Iz_0(n_zbeta, n_zdelta, d_beta + 2 * d_delta) = fac_bz * (sqrt(n_zdelta) * basis.tabzd(n_zbeta, n_zdelta - 1, d_beta, d_delta) - sqrt(n_zdelta + 1) * basis.tabzd(n_zbeta, n_zdelta + 1, d_beta, d_delta));
          Iz_0(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = -Iz_0(n_zbeta, n_zdelta, d_beta + 2 * d_delta);
        }//d_delta
      }//d_beta
    }//n_zdelta
  }//n_zbeta

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Iz_mp.
  Iz_mp = arma::zeros(n_zGlobalMax, n_zGlobalMax, DMax);

  for (INT n_zbeta = 0; n_zbeta < n_zGlobalMax ; n_zbeta++)
  {
    for (INT n_zdelta = 0; n_zdelta < n_zGlobalMax ; n_zdelta++)//symmetry  beta<=>delta
    {
      for (INT d_beta = 0; d_beta < dMax ; d_beta++)
      {
        for (INT d_delta = 0; d_delta < dMax ; d_delta++)
        {
          Iz_mp(n_zbeta, n_zdelta, d_beta + 2 * d_delta) = basis.tabzd(n_zbeta, n_zdelta, d_beta, d_delta);
          Iz_mp(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = Iz_mp(n_zbeta, n_zdelta, d_beta + 2 * d_delta);
        }//d_delta
      }//d_beta
    }//n_zdelta
  }//n_zbeta

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates the r spatial part.
 */

void FieldCDM2::calcIr(void)
{
  DBG_ENTER;

  if (!Ir_plus.empty()) DBG_LEAVE;

  //Usefull quantities
  INT M = basis.mMax;
  double fac_br = 1 / (basis.b_r * sqrt(2));

  //////////////////////////////////////////////////////////////////////////////////////////////////////Calculates Ir_plus.
  for (INT m_beta = 1;  m_beta < M ; m_beta ++)
  {
    Ir_plus(m_beta) = arma::zeros(basis.nMax(m_beta), basis.nMax(m_beta - 1));
    Ir_plus(m_beta)(basis.nMax(m_beta) - 1, basis.nMax(m_beta) - 1) = fac_br * sqrt(m_beta + INT(basis.nMax(m_beta)) - 1);

    if (basis.nMax(m_beta) < basis.nMax(m_beta - 1))
    {
      Ir_plus(m_beta)(basis.nMax(m_beta) - 1, basis.nMax(m_beta)) = fac_br * sqrt(INT(basis.nMax(m_beta)));
    }

    for (INT n_beta = 0; n_beta < basis.nMax(m_beta) - 1 ; n_beta++)
    {
      Ir_plus(m_beta)(n_beta, n_beta)   = fac_br * sqrt(m_beta + n_beta);
      Ir_plus(m_beta)(n_beta, n_beta + 1) = fac_br * sqrt(n_beta + 1);
    }//n_beta
  }//m_beta

  DBG_LEAVE;
}

