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

#include "global.h"
#include "general.h"
#include "datatree.h"
#include "action.h"
#include "interaction.h"
#include "solver.h"
#include "system.h"
#include "solver_alternator.h"
#include "solver_basis.h"
#include "solver_hfb_broyden.h"
#include "solver_hfb_gradient.h"
#include "solver_ws.h"
#include "state.h"

/** \file
 *  \brief Methods of the General class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > General::validKeys =
  {
    { "general/compatibility", "general.compatibility mode (can be: 'ringAndSchuck', 'berger', 'robledo', 'hfbtho')", "ringAndSchuck", "S" },
    { "general/version"      , "Version of the HFB3 library"                                                        , ""             , "S" },
    { "general/skill"        , "SKILL level (related to compilation optimization options) of the HFB3 library"      , ""             , "S" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The default constructor.
 */

General::General()
{
  version = CFG_GIT_VERSION;
  skill = SKILL;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree instance.
 */

General::General(const DataTree &dataTree)
{
  setGlobalValidKeys();

  // Set the general.compatibility mode.
  std::string compatibilityStr;
  dataTree.get(compatibilityStr, "general/compatibility", true);

  if      (compatibilityStr == "robledo")
  {
    compatibility = General::COMPAT_ROBLEDO;
    Tools::warning("Compatibility mode: 'robledo'");
  }
  else if (compatibilityStr == "berger" )
  {
    compatibility = General::COMPAT_BERGER;
    Tools::warning("Compatibility mode: 'berger'");
  }
  else if (compatibilityStr == "hfbtho" )
  {
    compatibility = General::COMPAT_HFBTHO;
    Tools::warning("Compatibility mode: 'hfbtho'");
  }

  if (dataTree.get(version, "general/version", true)) Tools::warning("Modifying the automatic version value");
  if (dataTree.get(skill,   "general/skill",   true)) Tools::warning("Modifying the automatic skill value");
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Gather the valid keys.
 */

void General::setGlobalValidKeys(void)
{
  if (globalValidKeys.size() == 0)
  {
    ADDLIST(globalValidKeys,            Action::validKeys);
    ADDLIST(globalValidKeys,             Basis::validKeys);
    ADDLIST(globalValidKeys,            Solver::validKeys);
    ADDLIST(globalValidKeys,       SolverBasis::validKeys);
    ADDLIST(globalValidKeys,  SolverAlternator::validKeys);
    ADDLIST(globalValidKeys,  SolverHFBBroyden::validKeys);
    ADDLIST(globalValidKeys, SolverHFBGradient::validKeys);
    ADDLIST(globalValidKeys,          SolverWS::validKeys);
    ADDLIST(globalValidKeys,        Constraint::validKeys);
    ADDLIST(globalValidKeys,       Interaction::validKeys);
    ADDLIST(globalValidKeys,            Mixing::validKeys);
    ADDLIST(globalValidKeys,             State::validKeys);
    ADDLIST(globalValidKeys,            System::validKeys);
    ADDLIST(globalValidKeys,           General::validKeys);
  }
}


