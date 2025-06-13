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

#ifndef FRAGMENTS_H
#define FRAGMENTS_H

/** \file
 *  \brief Headers for the Fragments class.
 */

#include "global.h"
#include "generic.h"
#include "state.h"
#include "geometry.h"

/** \brief Compute some fragment properties.
 *
 * This class allows to identify fragments and calculate some of their properties.
 */

class Fragments : public Generic
{
public :
  Fragments();
  Fragments(State &_state);                        // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string getNiceInfo(const std::string what = "");          // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// The id of the neck.
  INT izNeck = -1;

  std::vector<Geometry> geom = std::vector<Geometry>(3);

private:
};

#endif // FRAGMENTS_H
