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

#ifndef CLIPARSER_H
#define CLIPARSER_H

/** \file
 *  \brief Headers for the CliParser class.
 */

#include "global.h"
#include "generic.h"

/** \brief CLI parser.
 *
 * This class parses the CLI options.
 */

class CliParser : public Generic
{
public:

  // Constructor
  CliParser(int argc, char **argv);                                    // #TEST#

  static void help(void);
  static void version(void);
  static void printKeys(void);
  static void printInteractions(void);

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// The files to load
  std::vector<std::string> fileToLoad;

  /// The filename to save the DataTree.
  std::string saveTo;

  /// Will show the HFB3 logo at the start of the program if true.
  bool showLogo = false;
};

#endif // CLIPARSER_H
