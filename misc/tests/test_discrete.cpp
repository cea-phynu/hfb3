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
 *  \brief Test suite for the Discrete class.
 */

#include "gtest/gtest.h"
#include "discrete.h"
#include "global.h"
#include "state.h"

//==============================================================================

/*
TEST(Discrete, Discrete)
{
}
*/

//==============================================================================

TEST(Discrete, getLocalXZ)
{
  {
    // 1ct state
    State state("examples/42Ca_deformed_1x11.msg.gz");
    Discrete d(&state.basis);
    arma::mat densn = d.getLocalXZ(state.rho(NEUTRON), true);
    arma::mat densp = d.getLocalXZ(state.rho(PROTON ), true);

    double precision = 1e-12;

    // ALL_OUT;
    // INFO("%.6e", densn(0,  0));
    // INFO("%.6e", densn(5, 17));
    // INFO("%.6e", densp(0,  0));
    // INFO("%.6e", densp(5, 17));

    ASSERT_NEAR(densn(0,  0), 6.030257e-10, precision);
    ASSERT_NEAR(densn(5, 17), 3.970771e-06, precision);
    ASSERT_NEAR(densp(0,  0), 1.076740e-10, precision);
    ASSERT_NEAR(densp(5, 17), 1.066791e-06, precision);
  }

  {
    // 2ct state
    State state("examples/42Ca_deformed_2x9.msg.gz");
    Discrete d(&state.basis);
    arma::mat densn = d.getLocalXZ(state.rho(NEUTRON), true);
    arma::mat densp = d.getLocalXZ(state.rho(PROTON ), true);

    double precision = 1e-12;

    // ALL_OUT;
    // INFO("%.6e", densn(0,  0));
    // INFO("%.6e", densn(5, 17));
    // INFO("%.6e", densp(0,  0));
    // INFO("%.6e", densp(5, 17));

    ASSERT_NEAR(densn(0,  0), 1.172235e-07, precision);
    ASSERT_NEAR(densn(5, 17), 3.358336e-06, precision);
    ASSERT_NEAR(densp(0,  0), 1.786813e-08, precision);
    ASSERT_NEAR(densp(5, 17), 5.842514e-07, precision);
  }
}

//==============================================================================

/*
TEST(Discrete, getWaveFunctionXYZ)
{
}
*/

//==============================================================================

/*
TEST(Discrete, getHOXZ)
{
}
*/

//==============================================================================

TEST(Discrete, getHFmXZ)
{
  // // 2ct state
  // state state("examples/42Ca_deformed_2x9.msg.gz");
  // Discrete d(&state.basis);
  //
  // Discrete discrete(&state.basis, Mesh::regular(0, 0, -20, 10, 0, 20, 101, 1, 201));
  // //Discrete discrete(&state.basis, Mesh::gaussLaguerreHermite(120, 120));
  //
  // double ratio = 5; // test only 1 out of ratio values
  // arma::cube cn;
  //
  // for (UINT idHF = 0; idHF < state.HOtoHFn.n_rows; idHF++)
  // {
  //   if (Tools::irand(ratio) != 0) continue;
  //
  //   // NEUTRONS
  //   cn = discrete.getHFmXZ(state.HOtoHFn, idHF);
  //   arma::mat densP = arma::zeros(discrete.mesh.ax.nb, discrete.mesh.az.nb);
  //   arma::mat densM = arma::zeros(discrete.mesh.ax.nb, discrete.mesh.az.nb);
  //
  //   for (UINT i = 0; i < cn.n_slices; i++)
  //   {
  //     if (i % 2 == 0)
  //       densP += cn.slice(i);
  //     else
  //       densM += cn.slice(i);
  //   }
  //
  //   arma::mat funcP = arma::square(densP)
  //                     % (discrete.mesh.ax.we * discrete.mesh.az.we.t())
  //                     % Tools::matFromCol(discrete.mesh.ax.p, discrete.mesh.az.nb)
  //                     * 2.0 * PI;
  //   arma::mat funcM = arma::square(densM)
  //                     % (discrete.mesh.ax.we * discrete.mesh.az.we.t())
  //                     % Tools::matFromCol(discrete.mesh.ax.p, discrete.mesh.az.nb)
  //                     * 2.0 * PI;
  //   ASSERT_NEAR(arma::accu(funcP + funcM), 1.0, 3e-2); // a bit high ?
  // }
}

//==============================================================================

/*
TEST(Discrete, getLocalXYZ)
{
}
*/

//==============================================================================

/*
TEST(Discrete, calcDensit)
{
}
*/

//==============================================================================

/*
TEST(Discrete, calcWaveFunctions)
{
}
*/

//==============================================================================

/*
TEST(Discrete, getWaveFunctionXZ)
{
}
*/

//==============================================================================


TEST(Discrete, clear)
{
  Discrete discrete;

  Multi<arma::vec> mv;
  mv(1,2,3) = arma::randu(4);
  mv(1,8,9,7) = arma::randu(10);
  mv(11) = arma::randu(2);
  mv(1,33) = arma::vec();
  Multi<arma::mat> mv2;
  mv2(0) = arma::mat(5,6);
  Multi<arma::mat> mv3;
  mv3(11) = arma::randu(22,3);
  mv3(6,3) = arma::vec(9);

  discrete.zVals = mv;
  discrete.rVals = mv;
  discrete.rpVals = mv2;
  discrete.Densit = mv2;
  discrete.SODensit = mv3;
  discrete.BDensit = mv2;
  discrete.BSODensit = mv3;
  discrete.BSODensit1 = mv2;

  ASSERT_TRUE(!discrete.zVals.empty());
  ASSERT_TRUE(!discrete.rVals.empty());
  ASSERT_TRUE(!discrete.rpVals.empty());
  ASSERT_TRUE(!discrete.Densit.empty());
  ASSERT_TRUE(!discrete.SODensit.empty());
  ASSERT_TRUE(!discrete.BDensit.empty());
  ASSERT_TRUE(!discrete.BSODensit.empty());
  ASSERT_TRUE(!discrete.BSODensit1.empty());

  discrete.clear();

  ASSERT_TRUE(discrete.zVals.empty());
  ASSERT_TRUE(discrete.rVals.empty());
  ASSERT_TRUE(discrete.rpVals.empty());
  ASSERT_TRUE(discrete.Densit.empty());
  ASSERT_TRUE(discrete.SODensit.empty());
  ASSERT_TRUE(discrete.BDensit.empty());
  ASSERT_TRUE(discrete.BSODensit.empty());
  ASSERT_TRUE(discrete.BSODensit1.empty());
}

//==============================================================================

/*
TEST(Discrete, info)
{
}
*/

