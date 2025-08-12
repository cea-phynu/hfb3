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

#ifndef IO_HFB3_H
#define IO_HFB3_H

/** \file
 *  \brief Headers for the IOhfb3 class.
 */

#include "global.h"
#include "base64.h"
#include "gzstream.h"
#include "datatree.h"

/** \brief Methods to load .hfb3 files.
 *
 *  This class provides methos to load .hfb3 files.
 */

class IOhfb3
{
public :
  IOhfb3(void);                                                       // #TEST#
  DataTree fromContent(const std::string &content) const;             // #TEST#
  static void updateDataTree(DataTree &, const std::string key, const std::string val, std::string type);

  //============================================================================
  //============================================================================
  //============================================================================

};

#endif // IO_HFB3_H
