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

#ifndef IO_JSON_H
#define IO_JSON_H

/** \file
 *  \brief Headers for the IOjson class.
 */

#include "global.h"
class DataTree;

/** \brief Methods to load pseudo-JSON content.
 *
 *  This class provides methos to load pseudo-JSON content into a DataTree instance.
 */

class IOjson
{
public:

  IOjson(void);                                                        // #TEST#
  DataTree fromContent(const std::string &content) const;              // #TEST#
  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  void updateDataTree(DataTree &, const std::string key, const std::string val, const std::string type) const;
};

#endif // IO_JSON_H
