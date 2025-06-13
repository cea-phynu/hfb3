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

#include "states.h"
#include "tools.h"

/** \file
 *  \brief Methods of the States class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from an isospin value.
 */

States::States(const std::string _title) : title(_title)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Add a state.
 */

void States::add(INT _index, double _energy, double _occupation, const std::string &_label)

{
  // DBG_ENTER;

  IndivState indivState = {_index, _energy, _occupation, _label};
  states.push_back(indivState);
  n_elem = states.size();

  // DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Sort the states.
 */

void States::sort(void)

{
  DBG_ENTER;

  arma::vec energies = arma::zeros(states.size());

  INT i = 0;
  for (auto &state: states)
  {
    energies(i) = state.energy;
    i++;
  }

  UVEC sortIndex = arma::sort_index(energies, "ascend");

  std::vector<IndivState> sortedStates;

  for (UINT i = 0; i < sortIndex.n_elem; i++)
  {
    sortedStates.push_back(states[sortIndex(i)]);
  }

  states = sortedStates;
  n_elem = states.size();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Clear the states.
 */

void States::clear(void)
{
  DBG_ENTER;

  states.clear();
  n_elem = states.size();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string States::info(INT id, INT nbShown, bool isShort) const
{
  DBG_ENTER;

  std::string result;

  if (id == -1)
  {
    result = "None";
  }

  if (id >= 0)
  {
    if (isShort)
    {
      for (auto &state: states)
      {
        if (state.index == id)
        {
          result = PF("(%d,%s,%6.3f)", state.index, state.label.c_str(), state.energy);
          /* break; */
        }
      }
    }
    else
    {
      for (auto &state: states)
      {
        if (state.index == id)
        {
          result = Tools::treeStr({{"id", PF("%d", state.index)},
              {"ene", PF("%.3f", state.energy)},
              {"v^2", PF("%.3f", state.occupation)},
              {"label", PF("%s", state.label.c_str())},
              }, true);
          break;
        }
      }
    }
  }

  if (id == -2)
  {
    std::list<std::list<std::string> > tempList;

    INT nb = 0;
    for (auto &state: states)
    {
      std::string occupStr;
      if (state.occupation > 0.4)
      {
        occupStr = TABLE_BLUE + PF("%.6f", state.occupation) + TABLE_NORM;
      }
      else
      {
        occupStr = PF("%.6f", state.occupation);
      }

      if ((nbShown > -1) && (nb > nbShown))
      {
        tempList.push_back({"...", "...", "...", "..."});
        break;
      }
      else
      {
        tempList.push_back(
            {
            PF("%03d", state.index),
            PF("%.6f", state.energy),
            occupStr,
            state.label,
            });
      }
      nb++;
    } // state

    result += Tools::valueTable(
        title,
        {"energy", "v^2", "label"},
        {"[MeV]",  "[]",  "[]"},
        tempList
        );


    double occupSum = 0.0;
    tempList.clear();
    for (auto &state: states)
    {
      std::string occupStr;
      if (state.occupation > 0.4)
      {
        occupStr = TABLE_BLUE + PF("%.6f", state.occupation) + TABLE_NORM;
        tempList.push_back(
            {
            PF("%03d", state.index),
            PF("%.6f", state.energy),
            occupStr,
            state.label,
            });
      }
      occupSum += state.occupation;
    } // state

    tempList.push_back(
        {
        TABLE_RED + "Total",
        "",
        TABLE_RED + PF("%.6f", occupSum),
        ""
        });

    result += "\n";
    result += Tools::valueTable(
        title,
        {"energy", "v^2", "label"},
        {"[MeV]",  "[]",  "[]"},
        tempList
        );
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the index of the first [unoccupied] states.
 */

const IVEC States::getFirstEmptyStates(INT nbStates) const
{
  DBG_ENTER;

  IVEC result;

  for (auto &state: states)
  {
    // if (state.occupation > 0.5) continue;

    Tools::growIVec(result, state.index);
    if (result.n_elem >= nbStates) break;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Equality operator for States instances.
 */

bool operator==(const States &s0, const States &s1)
{
  return (s0.states == s1.states);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Equality operator for IndivState instances.
 */

bool operator==(const IndivState &s0, const IndivState &s1)
{
  return ((s0.index == s1.index) && (s0.label == s1.label));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Inequality operator for States instances.
 */

bool operator!=(const States &s0, const States &s1)
{
  return (s0.states != s1.states);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Inequality operator for IndivState instances.
 */

bool operator!=(const IndivState &s0, const IndivState &s1)
{
  return ((s0.index != s1.index) || (s0.label != s1.label));
}


