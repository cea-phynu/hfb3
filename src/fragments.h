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

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define FRAGMENTS_VALID_KEYS \
{ "fragments/neckPos" , "neck position [fm]"                           , "", "D" }, \
{ "fragments/neckDens", "neck density [fm-3]"                          , "", "D" }, \
{ "fragments/q_neck"  , "<Qneck> [fm-3]"                               , "", "D" }, \
{ "fragments/z_left"  , "protons in the left fragment"                 , "", "D" }, \
{ "fragments/n_left"  , "neutrons in the left fragment"                , "", "D" }, \
{ "fragments/a_left"  , "nucleons in the left fragment"                , "", "D" }, \
{ "fragments/z_right" , "protons in the right fragment"                , "", "D" }, \
{ "fragments/n_right" , "neutrons in the right fragment"               , "", "D" }, \
{ "fragments/a_right" , "nucleons in the right fragment"               , "", "D" }, \
{ "fragments/z_radius" , "mean radius of the proton local density [fm]", "", "D" }, \
{ "fragments/z_rms"    , "RMS of the proton local density [fm]"        , "", "D" }, \
{ "fragments/n_radius" , "mean radius of the proton local density [fm]", "", "D" }, \
{ "fragments/n_rms"    , "RMS of the proton local density [fm]"        , "", "D" }, \
{ "fragments/chargeRms", "charge RMS of the proton local density [fm]" , "", "D" }, \
{ "fragments/a_radius" , "mean radius of the proton local density [fm]", "", "D" }, \
{ "fragments/a_rms"    , "RMS of the proton local density [fm]"        , "", "D" }

//==============================================================================
//==============================================================================
//==============================================================================

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
  const DataTree getDataTree(void) const;                              // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// The id of the neck.
  INT izNeck = -1;

  std::vector<Geometry> geom = std::vector<Geometry>(3);

private:
};

#endif // FRAGMENTS_H
