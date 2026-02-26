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

void States::add(INT _index, double _energy, double _occupation, StateId _stateId)
{
  // DBG_ENTER;

  IndivState indivState = {_index, _energy, _occupation, _stateId};
  states.push_back(indivState);
  n_elem = states.size();
  // INFO("new state: " + Tools::infoStr(indivState));

  // DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Sort the states.
 */

void States::sort(const std::string &type)
{
  DBG_ENTER;

  VEC vals = arma::zeros(states.size());

  INT i = 0;
  for (auto &state: states)
  {
    if (type == "energy")     vals(i) = state.energy;
    if (type == "index")      vals(i) = double(state.index);
    if (type == "occupation") vals(i) = -1.0 * double(state.occupation);
    i++;
  }

  UVEC sortIndex = arma::sort_index(vals, "ascend");

  std::vector<IndivState> sortedStates;

  for (UINT i = 0; i < sortIndex.n_elem; i++)
  {
    sortedStates.push_back(states[sortIndex(i)]);
  }

  states = sortedStates;

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

VEC States::getEnergy(void) const
{
  VEC result = arma::zeros(n_elem);

  for (UINT i = 0; i < n_elem; i++)
  {
    result(states[i].index) = states[i].energy;
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

VEC States::getOccupation(void) const
{
  VEC result = arma::zeros(n_elem);

  for (UINT i = 0; i < n_elem; i++)
  {
    result(states[i].index) = states[i].occupation;
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

IVEC States::getIndex(void) const
{
  IVEC result = arma::zeros<IVEC>(n_elem);

  for (UINT i = 0; i < n_elem; i++)
  {
    result(states[i].index) = states[i].index;
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Find a state id from a stateId (seriously).
 */

UINT States::findState(const StateId &_stateId) const
{
  DBG_ENTER;

  for (auto &s: states)
  {
    if (s.stateId == _stateId) DBG_RETURN(s.index);
  }

  ERROR("State not found: " + Tools::color() + Tools::infoStr(_stateId));

  DBG_RETURN(0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string States::info(INT nbShown, bool onlyOccupiedStates) const
{
  DBG_ENTER;

  std::string result;
  std::list<std::list<std::string> > tempList;

  if (!onlyOccupiedStates)
  {
    INT nb = 0;
    for (auto &state: states)
    {
      std::string occupStr;
      if (state.occupation > 0.1)
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
            PF("(%d,%d/2)", state.stateId.index, state.stateId.omega * 2 + 1),
            });
      }
      nb++;
    } // state
  }
  else
  {
    double occupSum = 0.0;
    tempList.clear();
    for (auto &state: states)
    {
      std::string occupStr;
      if (state.occupation > 0.1)
      {
        occupStr = TABLE_BLUE + PF("%.6f", state.occupation) + TABLE_NORM;
        tempList.push_back(
            {
            PF("%03d", state.index),
            PF("%.6f", state.energy),
            occupStr,
            PF("(%d,%d/2)", state.stateId.index, state.stateId.omega * 2 + 1),
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

  }

  result += Tools::valueTable(
      title,
      {"energy", "v^2", "(id,omega)"},
      {"[MeV]",  "[]",  "[]"},
      tempList
      );

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string States::info(StateId _stateId, bool isShort) const
{
  DBG_ENTER;

  std::string result;

  UINT id = findState(_stateId);

  if (isShort)
  {
    for (auto &state: states)
    {
      if (state.index == id)
      {
        result = PF("(%d,%s,%6.3f)", state.index, Tools::infoStr(state.stateId).c_str(), state.energy);
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
            {"id" , PF("%s", Tools::infoStr(state.stateId).c_str())},
            }, true);
        break;
      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the index of the first [unoccupied] states.
 */

const std::list<StateId> States::getFirstEmptyStates(INT nbStates)
{
  DBG_ENTER;

  sort("energy");

  std::list<StateId> result;

  for (auto &state: states)
  {
    // if (state.occupation > 0.5) continue;

    result.push_back(state.stateId);
    if (result.size() >= nbStates) break;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check for emptiness.
 */

bool States::empty(void) const
{
  return states.empty();
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

/** Comparison operator for StateId instances.
 */

bool operator<(const StateId &s0, const StateId &s1)
{
  return ((s0.index < s1.index) || (s0.omega < s1.omega) || (s0.isospin < s1.isospin));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Equality operator for StateId instances.
 */

bool operator==(const StateId &s0, const StateId &s1)
{
  return ((s0.index == s1.index) && (s0.omega == s1.omega) && (s0.isospin == s1.isospin));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Equality operator for IndivState instances.
 */

bool operator==(const IndivState &s0, const IndivState &s1)
{
  return ((s0.index == s1.index) && (s0.stateId == s1.stateId));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Inequality operator for StateId instances.
 */

bool operator!=(const StateId &s0, const StateId &s1)
{
  return !(s0 == s1);
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
  return ((s0.index != s1.index) || (s0.stateId != s1.stateId));
}
