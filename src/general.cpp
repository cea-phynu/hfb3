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

/** \file
 *  \brief Methods of the General class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The default constructor.
 */

General::General()
{
#ifndef CFG_GIT_VERSION
    version = CFG_VERSION;
#else
    version = CFG_VERSION + std::string("(") + CFG_GIT_VERSION + std::string(")") ;
#endif

  skill = SKILL;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree instance.
 */

General::General(const DataTree &dataTree)
{

  // Set the general.compatibility mode.
  std::string compatibilityStr;
  dataTree.get(compatibilityStr, "general/compatibility", true);

  int _compatibility = General::COMPAT_NONE;
  if      (compatibilityStr == "robledo") _compatibility = General::COMPAT_ROBLEDO;
  else if (compatibilityStr == "berger" ) _compatibility = General::COMPAT_BERGER;
  else if (compatibilityStr == "hfbtho" ) _compatibility = General::COMPAT_HFBTHO;

  setCompatibility(_compatibility);

  if (dataTree.get(version, "general/version", true)) Tools::warning("Modifying the automatic version value");
  if (dataTree.get(skill,   "general/skill",   true)) Tools::warning("Modifying the automatic skill value");
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set the compatibility flag.
 */

void General::setCompatibility(int _compatibility)
{
  compatibility = _compatibility;

  switch(compatibility)
  {
    case COMPAT_NONE:
      // Tools::warning("Compatibility mode: 'ringAndSchuck'");
    break;
    case COMPAT_ROBLEDO:
      Tools::warning("Compatibility mode: 'robledo'");
    break;
    case COMPAT_BERGER:
      Tools::warning("Compatibility mode: 'berger'");
    break;
    case COMPAT_HFBTHO:
      Tools::warning("Compatibility mode: 'hfbtho'");
    break;
    default:
      ERROR("Unknown compatibility value.");
  }
}
