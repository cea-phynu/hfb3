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
 *  \brief Test suite for the State class.
 */

#include "gtest/gtest.h"
#include "state.h"
#include "constraint.h"
#include "multipole_operators.h"

//==============================================================================

void ASSERT_REL_ERROR(double a, double b, double c)
{
  if (fabs(a) > 1e-8)
    ASSERT_NEAR((a - b) / a, 0.0, c);
  else
    ASSERT_NEAR(a, b, c);
}

//==============================================================================

TEST(State, State)
{
  // Test State::constraints

  State s;
  s.constraints["q10t"] = Constraint("q10t", 1.1);
  s.constraints["q20t"] = Constraint("q20t", 2.2);
  s.constraints["q30t"] = Constraint("q30t", 3.3);

  ASSERT_EQ(s.constraints.size(), 3);
  ASSERT_EQ(s.constraints["q10t"].typeStr, "q10t");
  ASSERT_EQ(s.constraints["q20t"].typeStr, "q20t");
  ASSERT_EQ(s.constraints["q30t"].typeStr, "q30t");

  ASSERT_NEAR(s.constraints["q10t"].val, 1.1, 1e-12);
  ASSERT_NEAR(s.constraints["q20t"].val, 2.2, 1e-12);
  ASSERT_NEAR(s.constraints["q30t"].val, 3.3, 1e-12);

  s.constraints["q10t"] = Constraint("q10t", -1.1);
  s.constraints["q20t"] = Constraint("q20t", -2.2);
  s.constraints["q30t"] = Constraint("q30t", -3.3);
  s.constraints.erase("q20t");

  ASSERT_EQ(s.constraints.size(), 2);
  ASSERT_NEAR(s.constraints["q10t"].val, -1.1, 1e-12);
  ASSERT_NEAR(s.constraints["q30t"].val, -3.3, 1e-12);

  ASSERT_FALSE(s.converged);
}

//==============================================================================

TEST(State, State_)
{
  {
    // 1ct berger2ct state
    State s("examples/42Ca_deformed_1x11.msg.gz");

    Basis &basis = s.basis;

    ASSERT_TRUE(s.converged);
    ASSERT_NEAR(basis.b_r, 2.4144410198519912e+00, 1e-10);
    ASSERT_NEAR(basis.b_z, 3.1170299534289376e+00, 1e-10);
    ASSERT_NEAR(arma::trace(s.rho(NEUTRON)) * 2.0, 22.0, 1e-4);
    ASSERT_NEAR(arma::trace(s.rho(PROTON )) * 2.0, 20.0, 1e-5);
  }

  {
    // 2ct berger2ct state
    State s("examples/42Ca_deformed_2x9.msg.gz");

    Basis &basis = s.basis;

    ASSERT_NEAR(basis.b_r, 2.3878704919130063e+00, 1e-10);
    ASSERT_NEAR(basis.b_z, 3.1087425890539850e+00, 1e-10);
  }

  {
    // 1ct msg state
    State s("examples/42Ca_deformed_1x11.msg.gz");

    Basis &basis = s.basis;

    ASSERT_NEAR(basis.b_r, 2.4144410198519912e+00, 1e-10);
    ASSERT_NEAR(basis.b_z, 3.1170299534289376e+00, 1e-10);
    ASSERT_NEAR(arma::trace(s.rho(NEUTRON)) * 2.0, 22.0, 1e-4);
    ASSERT_NEAR(arma::trace(s.rho(PROTON )) * 2.0, 20.0, 1e-5);
  }

  {
    // 2ct msg state
    State s("examples/42Ca_deformed_2x9.msg.gz");

    Basis &basis = s.basis;

    ASSERT_NEAR(basis.b_r, 2.3878704919130063e+00, 1e-10);
    ASSERT_NEAR(basis.b_z, 3.1087425890539850e+00, 1e-10);
  }
}

//==============================================================================

TEST(State, State__)
{
  // With some values
  {
    // Construction of the datatree
    DataTree d;
    d.strict_mode = false;

    Multi<arma::mat> U;
    arma::mat u0 = arma::randu(3,3);
    arma::mat u1 = arma::randu(3,3);
    U(NEUTRON) = u0;
    U(PROTON ) = u1;
    Multi<arma::mat> V;
    arma::mat v0 = arma::randu(3,3);
    arma::mat v1 = arma::randu(3,3);
    V(NEUTRON) = v0;
    V(PROTON ) = v1;
    d.set("test/U", U);
    d.set("test/V", V);

    // Creation of the state
    State b(d);

    ASSERT_FALSE(b.rho.empty());
    ASSERT_FALSE(b.kappa.empty());
  }

  // Without values
  {
    DataTree dataTree = DataTree();
    State b(dataTree);

    ASSERT_TRUE(b.rho(NEUTRON).empty());
    ASSERT_TRUE(b.rho(PROTON).empty());
    ASSERT_TRUE(b.kappa(NEUTRON).empty());
    ASSERT_TRUE(b.kappa(PROTON).empty());
    ASSERT_TRUE(b.U(NEUTRON).empty());
    ASSERT_TRUE(b.U(PROTON).empty());
    ASSERT_TRUE(b.V(NEUTRON).empty());
    ASSERT_TRUE(b.V(PROTON).empty());

    ASSERT_TRUE(b.vecOc(NEUTRON).empty());
    ASSERT_TRUE(b.vecOc(PROTON ).empty());
    ASSERT_TRUE(b.eneQP(NEUTRON).empty());
    ASSERT_TRUE(b.eneQP(PROTON).empty());
    ASSERT_EQ(b.sys.nProt,  94);
    ASSERT_EQ(b.sys.nNeut, 146);
    ASSERT_EQ(b.totalEnergy, 1e99);
  }
}

//==============================================================================

TEST(State, RhoKappa)
{
  // // 1ct berger2ct state
  // State s("examples/42Ca_deformed_1x11.msg.gz");
  //
  // arma::mat rhon   = s.V(NEUTRON) * s.V(NEUTRON).t();
  // arma::mat rhop   = s.V(PROTON ) * s.V(PROTON ).t();
  // arma::mat kappan = s.V(NEUTRON) * s.U(NEUTRON).t();
  // arma::mat kappap = s.V(PROTON ) * s.U(PROTON ).t();
  //
  // ASSERT_TRUE(rhon == s.rho(NEUTRON));
  // ASSERT_TRUE(rhop == s.rho(PROTON ));
  // ASSERT_TRUE(kappan == s.kappa(NEUTRON));
  // ASSERT_TRUE(kappap == s.kappa(PROTON ));
}

//==============================================================================

TEST(State, convertFrom)
{
  {
    // 1ct - identical bases
    State sol1("examples/42Ca_deformed_1x11.msg.gz");
    Basis targetBasis("examples/42Ca_deformed_1x11.msg.gz");



    State sol2 = sol1;
    sol2.basis = targetBasis;
    sol2.convertFrom(sol1.basis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-14);
  }

  {
    // 1ct - bigger basis with identical deformation
    State sol1("examples/42Ca_deformed_1x11.msg.gz");

    double d_0 = sol1.basis.d_0;
    double b_r = sol1.basis.b_r;
    double b_z = sol1.basis.b_z;
    INT nOscil = sol1.basis.nOscil + 4;
    INT n_zMaxImposed = sol1.basis.n_zMaxImposed;
    double g_q = sol1.basis.g_q;
    Basis targetBasis(d_0, b_r, b_z, nOscil, n_zMaxImposed, g_q);

    State sol2 = sol1;
    sol2.basis = targetBasis;
    sol2.convertFrom(sol1.basis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-13);
  }

  {
    // 2ct - identical bases
    State sol1("examples/42Ca_deformed_2x9.msg.gz");
    Basis targetBasis("examples/42Ca_deformed_2x9.msg.gz");

    State sol2 = sol1;
    sol2.basis = targetBasis;
    sol2.convertFrom(sol1.basis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-14);
  }

}

//==============================================================================

TEST(State, convertTo)
{
  {
    // 1ct identical transformation
    State sol1("examples/42Ca_deformed_1x11.msg.gz");
    Basis targetBasis("examples/42Ca_deformed_1x11.msg.gz");

    State sol2 = sol1.convertTo(targetBasis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-14);
  }

  {
    // 1ct - bigger basis with identical deformation
    State sol1("examples/42Ca_deformed_1x11.msg.gz");

    double d_0 = sol1.basis.d_0;
    double b_r = sol1.basis.b_r;
    double b_z = sol1.basis.b_z;
    INT nOscil = sol1.basis.nOscil + 2;
    INT n_zMaxImposed = sol1.basis.n_zMaxImposed;
    double g_q = sol1.basis.g_q;
    Basis targetBasis(d_0, b_r, b_z, nOscil, n_zMaxImposed, g_q);

    State sol2 = sol1.convertTo(targetBasis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-13);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-13);
  }

  {
    // 2ct identical transformation
    State sol1("examples/42Ca_deformed_2x9.msg.gz");
    Basis targetBasis("examples/42Ca_deformed_2x9.msg.gz");

    State sol2 = sol1.convertTo(targetBasis);

    MultipoleOperators multipoleOperators1(sol1);
    MultipoleOperators multipoleOperators2(sol2);

    ASSERT_REL_ERROR(multipoleOperators1.qlm(1, 0), multipoleOperators2.qlm(1, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(2, 0), multipoleOperators2.qlm(2, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(3, 0), multipoleOperators2.qlm(3, 0), 1e-14);
    ASSERT_REL_ERROR(multipoleOperators1.qlm(4, 0), multipoleOperators2.qlm(4, 0), 1e-14);
  }
}

//==============================================================================


TEST(State, getDataTree)
{
  State bogo;
  bogo.sys.nProt = 90;
  bogo.sys.nNeut = 140;
  Basis basis;
  DataTree d = bogo.getDataTree();

  DataTree result;
  result.merge(basis.getDataTree("state/"));


  result.set("system/nProt", 90);
  result.set("system/nNeut", 140);
  result.set("state/rho", Multi<arma::mat>());
  result.set("state/kappa", Multi<arma::mat>());
  result.set("state/individualStates/occupation", Multi<arma::vec>());
  result.set("state/individualStates/energy", Multi<arma::vec>());
  result.set("state/chemicalPotential", (arma::vec)arma::zeros(2));
  result.set("state/totalEnergy", 1e99);
  result.set("state/converged", false);
  result.set("state/calculationLength", 0.0);
  result.set("state/nbIter", 0);

  result.validate();

  ASSERT_TRUE(d == result);

  auto dt = basis.getDataTree("state/");
  dt.set("state/basis/d_0", 9.8);
  result.merge(dt);
  ASSERT_FALSE(d == result);
}

//==============================================================================

/*
TEST(State, getDensityDistance)
{
}
*/

//==============================================================================

/*
TEST(State, info)
{
}
*/
