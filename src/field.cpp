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

#include "field.h"
#include "tools.h"
#include "interaction.h"

/** \file
 *  \brief Methods of the Field class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

const std::vector<std::string> Field::typeStr = {"Direct", "Exchange", "Pairing"};

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

Field::Field(Parameters _parameters, State *_state):
  state_ptr(_state), // TODO: if possible, remove this pointer, use the reference instead
  state(*state_ptr),
  parameters(_parameters),
  basis(state.basis),
  rho(state.rho),
  kappa(state.kappa)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set the deformation parameters (used by FieldWoodsSaxon class).
 */

void Field::setDef(const Multi<double> &)
{
  DBG_ENTER;
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get the Woods-Saxon potential (used by FieldWoodsSaxon class).
 */

WSPot Field::getWSPot(void)
{
  DBG_ENTER;
  DBG_RETURN(WSPot());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the field (must be defined by the children classes).
 */

void Field::calcField(void)
{
  DBG_ENTER;

  ASSERT(!rho(NEUTRON  ).empty(), "empty rho(NEUTRON) matrix");
  ASSERT(!rho(PROTON   ).empty(), "empty rho(PROTON) matrix");
  ASSERT(!kappa(NEUTRON).empty(), "empty kappa(NEUTRON) matrix");
  ASSERT(!kappa(PROTON ).empty(), "empty kappa(PROTON) matrix");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the field (must be defined by the children classes).
 */

void Field::clear(void)
{
  DBG_ENTER;

  field.clear();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The method to calculate the energy.
 */

void Field::calcEnergy()
{
  DBG_ENTER;

  ASSERT(!rho(NEUTRON  ).empty(), "empty rho(NEUTRON) matrix");
  ASSERT(!rho(PROTON   ).empty(), "empty rho(PROTON) matrix");
  ASSERT(!kappa(NEUTRON).empty(), "empty kappa(NEUTRON) matrix");
  ASSERT(!kappa(PROTON ).empty(), "empty kappa(PROTON) matrix");

  for (INT type = 0; type < 3; type++)
  {
    for (INT iso: {NEUTRON, PROTON})
    {
      if (!field.contains(iso, type)) continue;

      ASSERT(!field(iso, type).empty(), "empty field: " + name);

      energy(iso, type) = 0.0;

      for (INT omega = 0; omega < basis.mMax; omega++)
      {
        if (type == PAIRING)
        {
          energy(iso, type) += arma::trace(field(iso, type).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega))
                                           * kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)));
        }
        else
        {
          energy(iso, type) += arma::trace(field(iso, type).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega))
                                           * rho(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)));
        }
      }

      energy(iso, type) *= energyFactor;
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Field::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  std::vector<std::string> isoStr = {"n", "p"};

  if (isShort)
  {
    std::vector<std::string> typeStr = {"D", "E", "P"};

    std::vector<std::pair<std::string, double> > pairList;

    for (INT type = 0; type < 4; type++)
    {
      for (INT iso = 0; iso < 3; iso++)
      {
        if (energy.contains(iso, type))
        {
          pairList.push_back({typeStr[type] + isoStr[iso], energy(iso, type)});
        }
      }
    }

    result += "(";

    bool ifirst = true;

    for (auto &p : pairList)
    {
      if (ifirst) ifirst = false;
      else result += ", ";

      result += PF_YELLOW(p.first) + ":" + PF_BLUE("%12.6f", p.second);
    }

    result += ")";

    DBG_RETURN(result)
  }


  // 1. Energy contribution info
  std::vector<std::string> typeStr = {"Dire", "Exch", "Pair"};
  std::vector<std::pair<std::string, std::string> > strList;

  strList.push_back({"Field", ""});
  strList.push_back({"name  ", name});
  // strList.push_back({"nGLA  ", PF("%d", nGLA)});
  // strList.push_back({"nGHE  ", PF("%d", nGHE)});
  // strList.push_back({"nGLE  ", PF("%d", nGLE)});

  for (INT iso = 0; iso < 3; iso++)
  {
    for (INT type = 0; type < 4; type++)
    {
      if (energy.contains(iso, type))
      {
        std::string fieldLabel = typeStr[type] + "_" + isoStr[iso];
        std::string energyStr = PF("%13.6f", energy(iso, type));

        strList.push_back({fieldLabel, energyStr});
      }
    }
  }
  result += Tools::treeStr(strList, isShort) + "\n";

  // 2. Ram Info
  result += ramInfo() + "\n";

  // 3. Timing Info
  result += timeInfo() + "\n";

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Accessor to the energy contributions.
 */

double Field::getEnergy(INT iso, INT type)
{
  // TODO: replace me with Python bindings for the Multi<> class
  return energy(iso, type);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Field::ramInfo() const
{
  std::string table;

  table +=
    TABLE_YELLOW + TABLE_LEFT + name + TABLE_TD
  + TABLE_YELLOW + TABLE_LEFT + "Size"  + TABLE_TD + TABLE_TR;

  UINT s = 0;
  for (auto &o : ramTable)
  {
    table +=

      TABLE_NORM  + o.first + TABLE_TD
    + TABLE_GREEN + Tools::humanSize(o.second) + TABLE_TD + TABLE_TR;

    s += o.second;
  }

  table +=
      TABLE_RED  + "TOTAL" + TABLE_TD
    + TABLE_RED + Tools::humanSize(s) + TABLE_TD + TABLE_TR;

  return Tools::printTable(table);

}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Field::timeInfo() const
{
  std::string table;

  table +=
    TABLE_YELLOW + TABLE_LEFT + name + TABLE_TD
  + TABLE_YELLOW + TABLE_LEFT + "Time (s)"  + TABLE_TD + TABLE_TR;

  double s = 0.;
  for (auto &o : timeTable)
  {
    table +=
      TABLE_NORM  + o.first + TABLE_TD
    + TABLE_GREEN + PF("%f", o.second) + TABLE_TD + TABLE_TR;

    s += o.second;
  }

  return Tools::printTable(table);

}

