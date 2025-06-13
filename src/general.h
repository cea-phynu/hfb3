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

#ifndef GENERAL_H
#define GENERAL_H

/** \file
 *  \brief Headers for the General class.
 */

#include <map>
#include <string>
#include <list>

class DataTree;

/** \brief A struct to store ax key with its description, default value and type.
 */

struct KeyStruct
{
  std::string key;
  std::string description;
  std::string defaultValue;
  std::string type;
};

/** \brief General parameters.
 *
 *  This class stores general parameters.
 */

class General
{
public:

  /// Constructors
  General(void);                                                       // #TEST#
  General(const DataTree &dataTree);                                   // #TEST#

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// Possible compatibility values,
  enum
  {
    COMPAT_NONE,
    COMPAT_BERGER,
    COMPAT_ROBLEDO,
    COMPAT_HFBTHO
  };

  /// Compatibility mode.
  int compatibility = COMPAT_NONE;

  /// Version of the HFB3 library.
  std::string version = CFG_GIT_VERSION;

  /// SKILL level (related to compilation optimization options) of the HFB3 library.
  std::string skill = SKILL;
};

#endif // GENERAL_H

