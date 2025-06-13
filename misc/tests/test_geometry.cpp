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

/** \file
 *  \brief Test suite for the Geometry class.
 */

#include "gtest/gtest.h"
#include "geometry.h"
#include "state.h"
#include "discrete.h"

//==============================================================================

TEST(Geometry, Geometry)
{
  {
    // 1ct state
    State s("examples/42Ca_deformed_1x11.msg.gz");
    Discrete discrete(&s.basis, Mesh::regular(0, 0, -15, 15, 0, 15, 101, 1, 301));
    arma::mat rhot = s.rho(NEUTRON) + s.rho(PROTON );
    arma::mat denst = discrete.getLocalXZ(rhot, true);
    Geometry geomt(discrete.mesh, denst, s.sys);
    ASSERT_NEAR(geomt.intTotal, 42.0, 3e-2);

    // This works with an orthonormal basis (not with 2ct)
    double traceRho = arma::trace(s.rho(NEUTRON)) + arma::trace(s.rho(PROTON ));
    ASSERT_NEAR(traceRho, 21.0, 1e-4);
  }

  {
    // 2ct state
    State s("examples/42Ca_deformed_2x9.msg.gz");
    Discrete discrete(&s.basis, Mesh::regular(0, 0, -15, 15, 0, 15, 101, 1, 301));
    arma::mat rhot = s.rho(NEUTRON) + s.rho(PROTON );
    arma::mat denst = discrete.getLocalXZ(rhot, true);
    Geometry geomt(discrete.mesh, denst, s.sys);
    ASSERT_NEAR(geomt.intTotal, 42.0, 5e-2);
  }
}

//==============================================================================

TEST(Geometry, Geometry_)
{
  {
    // 1ct state
    State state("examples/42Ca_deformed_1x11.msg.gz");
    Geometry geomt(state);
    ASSERT_NEAR(geomt.intTotal, 42.0, 3e-2);
  }

  {
    // 2ct state
    State state("examples/42Ca_deformed_2x9.msg.gz");
    Geometry geomt(state);
    ASSERT_NEAR(geomt.intTotal, 42.0, 5e-2);
  }
}

//==============================================================================

/*
TEST(Geometry, info)
{
}
*/
