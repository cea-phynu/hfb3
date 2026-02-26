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

#ifndef STATES_H
#define STATES_H


/** \file
 *  \brief Headers for the States class.
 */

#include "global.h"
#include "generic.h"

//============================================================================
//============================================================================
//============================================================================

/** \brief Identification of a state.
 */

struct StateId
{
  INT index;
  INT omega;
  INT isospin;
};

//============================================================================
//============================================================================
//============================================================================

/** \brief Store one individual state.
 */

struct IndivState
{
  INT index;
  double energy;
  double occupation;
  StateId stateId;
};

//============================================================================
//============================================================================
//============================================================================

/** \brief Store properties associated to a list of states.
 *
 * This class stores the properties associated to a list of states.
 */

class States : public Generic
{
public :
  States(const std::string _title = " ");
  void sort(const std::string& type = "energy");
  void add(INT _index,
           double _energy,
           double _occupation,
           StateId _stateId);
  const std::string info(StateId _stateId, bool = false) const;
  const std::string info(INT nbShown = 25, bool onlyOccupiedStates = false) const;
  void clear(void);
  const std::list<StateId> getFirstEmptyStates(INT nbStates = 1);
  VEC getEnergy(void) const;
  IVEC getIndex(void) const;
  VEC getOccupation(void) const;
  UINT findState(const StateId &) const;
  bool empty(void) const;

  //============================================================================
  //============================================================================
  //============================================================================

  // Title of the table.
  std::string title;

  std::vector<IndivState> states;

  /// Number of states
  int n_elem = 0;

private:

};

extern bool operator==(const States &, const States &);
extern bool operator==(const IndivState &, const IndivState &);
extern bool operator!=(const States &, const States &);
extern bool operator!=(const IndivState &, const IndivState &);
extern bool operator<(const StateId &, const StateId &);
extern bool operator==(const StateId &, const StateId &);
extern bool operator!=(const StateId &, const StateId &);

#endif // STATES_H

