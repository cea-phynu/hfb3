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

#include "fragments.h"
#include "discrete.h"

/** \file
 *  \brief Methods of the Fragments class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 * \param state Reference to a State object.
 */

Fragments::Fragments(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

Fragments::Fragments(State &state)
{
  DBG_ENTER;

  if (!state.rho(NEUTRON).empty())
  {
    Discrete discrete(state.basis, Mesh::regular(0, 0, -20, 20, 0, 20, 401, 1, 401));
    arma::mat densn = discrete.getLocalXZ(state.rho(NEUTRON), true);
    arma::mat densp = discrete.getLocalXZ(state.rho(PROTON ), true);
    arma::mat denst = densn + densp;

    // Numerical calculation of the mass multipole moments and of the fagment properties
    Geometry geomt(discrete.mesh, denst, state.sys);
    izNeck = geomt.izNeck;

    geom[NEUTRON] = Geometry(discrete.mesh, densn, state.sys, izNeck);
    geom[PROTON ] = Geometry(discrete.mesh, densp, state.sys, izNeck);
    geom[TOTAL  ] = Geometry(discrete.mesh, denst, state.sys, izNeck);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Fragments::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  std::vector<std::pair<std::string, std::string> > list;

  if (isShort)
  {
    DBG_RETURN(PF_RED("empty"));

    // TODO
  }
  else
  {
    list.push_back({"Fragments", ""});

    // TODO
  }

  result += Tools::treeStr(list, isShort);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice table with the fragment properties.
 */

const std::string Fragments::getNiceInfo(const std::string what)
{
  DBG_ENTER;

  std::string result;

  if (what == "radii")
  {
    // NUCELAR RADII
    result += Tools::valueTable( "Radii",
        {"RMS", "<r>", "ch. RMS"},
        {"[fm]", "[fm]", "[fm]"},
        {
        {"neutron", PF("%.3f", geom[NEUTRON].rms), PF("%.3f", geom[NEUTRON].radius), " "},
        {"proton",  PF("%.3f", geom[PROTON ].rms), PF("%.3f", geom[PROTON ].radius), PF("%.3f", geom[PROTON].chargeRMS)},
        {"total",   PF("%.3f", geom[TOTAL  ].rms), PF("%.3f", geom[TOTAL  ].radius), " "},
        }
        );
  }
  else
  {
    if (izNeck != -1)
    {
      result += Tools::valueTable( "Fragments",
      {"N. pos.", "N. dens.", "<Qneck>", "L. Mass", "R. Mass", "T. Mass"},
      {"[fm]", "[n.fm-3]", "[n]", "[n]", "[n]", "[n]"    },
      {
        {
          "neutron",
          " ",
          PF("%.3f", geom[NEUTRON].neckDens),
          PF("%.3f", geom[NEUTRON].q_neck  ),
          PF("%.3f", geom[NEUTRON].intLeft ),
          PF("%.3f", geom[NEUTRON].intRight),
          PF("%.3f", geom[NEUTRON].intTotal),

        },
        {
          "proton",
          " ",
          PF("%.3f", geom[PROTON].neckDens),
          PF("%.3f", geom[PROTON].q_neck  ),
          PF("%.3f", geom[PROTON].intLeft ),
          PF("%.3f", geom[PROTON].intRight),
          PF("%.3f", geom[PROTON].intTotal),

        },
        {
          "total",
          PF("%.3f", geom[TOTAL].neckPos ),
          PF("%.3f", geom[TOTAL].neckDens),
          PF("%.3f", geom[TOTAL].q_neck  ),
          PF("%.3f", geom[TOTAL].intLeft ),
          PF("%.3f", geom[TOTAL].intRight),
          PF("%.3f", geom[TOTAL].intTotal),

        },
      }
                                 );
    }
    else
    {
      result += Tools::valueTable( "Fragments",
      {"N. pos."},
      {"[fm]"},
      {
        {"neutron", " "},
        {"proton", " "},
        {"total", TABLE_RED + "no neck found"}
      }
                                 );
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Save the results in a DataTree instance.
 */

const DataTree Fragments::getDataTree(void) const
{
  DBG_ENTER;

  DataTree dt;

  dt.set("fragments/neckPos"  , geom[TOTAL  ].neckPos);
  dt.set("fragments/neckDens" , geom[TOTAL  ].neckDens);
  dt.set("fragments/q_neck"   , geom[TOTAL  ].q_neck);
  dt.set("fragments/z_left"   , geom[PROTON ].intLeft);
  dt.set("fragments/n_left"   , geom[NEUTRON].intLeft);
  dt.set("fragments/a_left"   , geom[TOTAL  ].intLeft);
  dt.set("fragments/z_right"  , geom[PROTON ].intRight);
  dt.set("fragments/n_right"  , geom[NEUTRON].intRight);
  dt.set("fragments/a_right"  , geom[TOTAL  ].intRight);

  dt.set("fragments/z_radius" , geom[PROTON ].radius);
  dt.set("fragments/z_rms"    , geom[PROTON ].rms);
  dt.set("fragments/n_radius" , geom[NEUTRON].radius);
  dt.set("fragments/n_rms"    , geom[NEUTRON].rms);
  dt.set("fragments/chargeRms", geom[PROTON ].chargeRMS);
  dt.set("fragments/a_radius" , geom[TOTAL  ].radius);
  dt.set("fragments/a_rms"    , geom[TOTAL  ].rms);

  DBG_RETURN(dt);
}
