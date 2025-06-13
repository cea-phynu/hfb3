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

#include "io_amedee.h"
#include "tools.h"
#include "solver_hfb_broyden.h"

/** \file
 *  \brief Methods of the IOamedee class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

IOamedee::IOamedee(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Read basis parameters and \f$\rho_n\f$, \f$\rho_p\f$ matrices, then create a Basis instance.
 */

void IOamedee::readBasisRho(const std::string &content)
{
  DBG_ENTER;

  if (!checkFileType(content))
  {
    Tools::warning("wrong format in IOamedee::readBasisRho()");
    DBG_LEAVE;
  }

  std::stringstream stream;
  stream << content;
  std::string inputStr;
  std::getline(stream, inputStr); // '#amedee_dpfp9'
  std::getline(stream, inputStr); // Z
  Z = atoi(inputStr.c_str());
  std::getline(stream, inputStr); // N
  N = atoi(inputStr.c_str());
  std::getline(stream, inputStr); // nOscil
  INT nOscil = INT(atol(inputStr.c_str()));
  std::getline(stream, inputStr); // alpha
  double alpha = atof(inputStr.c_str());
  std::getline(stream, inputStr); // beta
  double beta = atof(inputStr.c_str());
  std::getline(stream, inputStr); // hbar omega
  std::getline(stream, inputStr); // Q (deformation for the basis truncation)
  double gq_amedee = atof(inputStr.c_str());
  std::getline(stream, inputStr); // q (used for the basis truncation ?)
  std::getline(stream, inputStr); // number of gaussian nodes [z-axis]
  std::getline(stream, inputStr); // number of gaussian nodes [r-axis]
  double br = 1.0 / sqrt(beta);
  double bz = 1.0 / sqrt(alpha);

  if (fabs(gq_amedee - 1.0) > 1e-9) Tools::warning("GQ read from AMEDEE file is different from 1.0 !");

  basis = Basis(-1.0, br, bz, nOscil, 24, 1.0 / pow(gq_amedee, 2));
  rho(NEUTRON)   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  rho(PROTON )   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappa(NEUTRON) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappa(PROTON ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  customFieldn   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  customFieldp   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  while (!stream.eof())
  {
    std::getline(stream, inputStr);

    if (inputStr.length() == 0) continue;

    if (inputStr[0] == '#')
    {
      std::vector<std::string> elems = Tools::stringSplit(Tools::trim(inputStr), ':');

      if (elems.size() == 2)
      {
        if (elems[0] == "#kinetic_n"          ) energiesNeut["Kinetic"         ] = atof(elems[1].c_str());

        if (elems[0] == "#kinetic_p"          ) energiesProt["Kinetic"         ] = atof(elems[1].c_str());

        if (elems[0] == "#coulomb_direct_p"   ) energiesProt["Coulomb_Direct"  ] = atof(elems[1].c_str());

        if (elems[0] == "#central_direct_n"   ) energiesNeut["Central_Direct"  ] = atof(elems[1].c_str());

        if (elems[0] == "#central_direct_p"   ) energiesProt["Central_Direct"  ] = atof(elems[1].c_str());

        if (elems[0] == "#central_exchange_n" ) energiesNeut["Central_Exchange"] = atof(elems[1].c_str());

        if (elems[0] == "#central_exchange_p" ) energiesProt["Central_Exchange"] = atof(elems[1].c_str());

        if (elems[0] == "#spin-orbit_n"       ) energiesNeut["Spin-orbit"      ] = atof(elems[1].c_str());

        if (elems[0] == "#spin-orbit_p"       ) energiesProt["Spin-orbit"      ] = atof(elems[1].c_str());

        if (elems[0] == "#density_n"          ) energiesNeut["Density_D2"      ] = atof(elems[1].c_str());

        if (elems[0] == "#density_p"          ) energiesProt["Density_D2"      ] = atof(elems[1].c_str());

        if (elems[0] == "#rearrangement_n"    ) energiesNeut["Rearrangement_D2"] = atof(elems[1].c_str());

        if (elems[0] == "#rearrangement_p"    ) energiesProt["Rearrangement_D2"] = atof(elems[1].c_str());
      }

      continue;
    }

    double vrho, vkappa;
    INT ma, na, nza, sa;
    INT mc, nc, nzc, sc;
    INT itype, iso;
    std::stringstream ss(inputStr);
    ss >> itype;
    ss >> iso;
    ss >> ma;
    ss >> na;
    ss >> nza;
    ss >> sa;
    ss >> mc;
    ss >> nc;
    ss >> nzc;
    ss >> sc;
    ss >> vrho;
    ss >> vkappa;
    INT ia, ic;
    ia = basis.HOqn.find({ma, na, nza, 0, sa});

    ASSERT(ia >= 0, "state not found: m=%d n=%d n_z=%d s=%d", ma, na, nza, sa);

    ic = basis.HOqn.find({mc, nc, nzc, 0, sc});

    ASSERT(ic >= 0, "state not found: m=%d n=%d n_z=%d s=%d", mc, nc, nzc, sc);

    if (iso == 1)
    {
      if (itype == 0)
      {
        rho(PROTON)(ia, ic) = vrho;
        rho(PROTON)(ic, ia) = vrho;
        kappa(PROTON)(ia, ic) = vkappa;
        kappa(PROTON)(ic, ia) = vkappa;
      }
    }
    else
    {
      if (itype == 0)
      {
        rho(NEUTRON)(ia, ic) = vrho;
        rho(NEUTRON)(ic, ia) = vrho;
        kappa(NEUTRON)(ia, ic) = vkappa;
        kappa(NEUTRON)(ic, ia) = vkappa;
      }
      else if (itype == 1)
      {
        customFieldn(ia, ic) = vrho;
        customFieldn(ic, ia) = vrho;
        customFieldp(ia, ic) = vkappa;
        customFieldp(ic, ia) = vkappa;
      }
    }
  }

  // Tools::info("customFieldn", customFieldn, true);
  // Tools::info("customFieldp", customFieldp, true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance from AMEDEE formatted content.
 */

DataTree IOamedee::fromContent(const std::string &content)
{
  DBG_ENTER;
  DataTree result;

  if (!checkFileType(content))
  {
    DBG_RETURN(result);
  }

  readBasisRho(content);

  result.merge(basis.getDataTree());
  result.set("state/rho", rho);
  result.set("state/kappa", kappa);
  // result.set("state/customFieldNeut", customFieldn);
  // result.set("state/customFieldProt", customFieldp);
  result.set("system/nNeut", N);
  result.set("system/nProt", Z);
  result.set("interaction/name", "D1S");

  // for (auto &ene: energiesNeut)
  // {
  //   result.set("referenceEnergiesNeut/" + ene.first, ene.second);
  // }
  // for (auto &ene: energiesProt)
  // {
  //   result.set("referenceEnergiesProt/" + ene.first, ene.second);
  // }


  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check file type
 */

bool IOamedee::checkFileType(const std::string &content)
{
  DBG_ENTER;

  std::stringstream stream;
  stream << content;
  std::string inputStr;
  std::getline(stream, inputStr);

  if (inputStr != "#amedee_dpfp9")
  {
    DBG_RETURN(false);
  }

  DBG_RETURN(true);
}


