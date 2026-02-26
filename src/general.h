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

//==============================================================================
//==============================================================================
//==============================================================================

class DataTree;

//==============================================================================
//==============================================================================
//==============================================================================

/// Helper macro to append to an std::list).
#define ADDLIST(A,B) A.insert(A.end(), B.begin(), B.end())

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief A struct to store a key with its description, default value and type.
 */

#define GENERAL_VALID_KEYS \
{ "general/compatibility", "general.compatibility mode (can be: 'ringAndSchuck', 'berger', 'robledo', 'hfbtho')", "ringAndSchuck", "S" }, \
{ "general/version"      , "Version of the HFB3 library"                                                        , ""             , "S" }, \
{ "general/skill"        , "SKILL level (related to compilation optimization options) of the HFB3 library"      , ""             , "S" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief General parameters.
 *
 *  This class stores general parameters.
 */

class General
{
public:

  //============================================================================
  //============================================================================
  //============================================================================

  /// Possible compatibility values,
  enum
  {
    COMPAT_NONE,
    COMPAT_BERGER,
    COMPAT_ROBLEDO,
    COMPAT_HFBTHO
  };

  //============================================================================
  //============================================================================
  //============================================================================

  /// Constructors
  General(void);                                                       // #TEST#
  General(const DataTree &dataTree);                                   // #TEST#

  void setCompatibility(int);

  //============================================================================
  //============================================================================
  //============================================================================

  /// Compatibility mode.
  int compatibility = COMPAT_NONE;

  //============================================================================
  //============================================================================
  //============================================================================

#ifdef SKILL0
/// SKILL level 0 (cf. Makefile for meaning)
#define SKILL "SKILL0 (I'm too young to die)"
#endif

#ifdef SKILL1
/// SKILL level 1 (cf. Makefile for meaning)
#define SKILL "SKILL1 (Hey, not too rough)"
#endif

#ifdef SKILL2
/// SKILL level 2 (cf. Makefile for meaning)
#define SKILL "SKILL2 (Hurt me plenty)"
#endif

#ifdef SKILL3
/// SKILL level 3 (cf. Makefile for meaning)
#define SKILL "SKILL3 (Ultra-Violence)"
#endif

#ifdef SKILL4
/// SKILL level 4 (cf. Makefile for meaning)
#define SKILL "SKILL4 (Nightmare!)"
#endif


/// Default Git version string (should be set at compile time).
#ifndef CFG_GIT_VERSION
#define CFG_GIT_VERSION ""
#endif

  /// Version of the HFB3 library.
  std::string version = CFG_GIT_VERSION;

// Default SKILL string value (should be set above).
#ifndef SKILL
#define SKILL "NO SKILL DEFINED !?"
#endif

  /// SKILL level (related to compilation optimization options) of the HFB3 library.
  std::string skill = SKILL;

  /// Flag to display the compatibility warning once.
  bool compatWarningDisplayed = false;
};


#endif // GENERAL_H

