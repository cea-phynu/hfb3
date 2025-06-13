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

/// @cond BASE64

#ifndef BASE64_H
#define BASE64_H

#include <string>

/** \file
 *  \brief Headers for the Base64 namespace.
 */

/** \brief A namespace containing routines to encode and decode in base64 format.
 *
 * This namespace contains routines to encode and decode in base64 format.
 */

namespace Base64
{
  std::string encode(char const *, long unsigned int len);             // #TEST#
  std::string encode(std::string const &s);                            // #TEST#
  std::string decode(std::string const &s);                            // #TEST#
}

#endif // BASE64_H

/// @endcond
