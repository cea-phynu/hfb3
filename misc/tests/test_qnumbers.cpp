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
 *  \brief Test suite for the Qnumbers class.
 */

#include "gtest/gtest.h"
#include "qnumbers.h"
#include "basis.h"
#include "state.h"

//==============================================================================


TEST(Qnumbers, Qnumbers)
{
  Qnumbers q1(0);
  ASSERT_TRUE(q1.nbq == 0);
  ASSERT_TRUE(q1.nb == 0);
  ASSERT_TRUE(q1.name == std::vector<std::string>(0));

  Qnumbers q2(22);
  ASSERT_TRUE(q2.nbq == 22);
  ASSERT_TRUE(q2.nb == 0);
  ASSERT_TRUE(q2.name == std::vector<std::string>(22));
}

//==============================================================================


TEST(Qnumbers, Qnumbers_)
{
  Qnumbers q1("n_z", 2);
  IMAT mat(1, 2);
  for (INT i = 0; i < 2; i++)
  {
    mat(0, i) = i;
  }

  ASSERT_TRUE(q1.nbq == 1);
  ASSERT_TRUE(q1.nb == 2);
  ASSERT_TRUE(q1.name[0] == "n_z");
  ASSERT_TRUE(q1.mat == mat);
}

//==============================================================================


TEST(Qnumbers, appendName)
{
  Qnumbers q(2);
  q.appendName("q2");

  ASSERT_TRUE(q.name[2] == "q2");
  ASSERT_TRUE(q.nbq == 3);
}

//==============================================================================

/*
TEST(Qnumbers, sub)
{
}
*/

//==============================================================================

/*
TEST(Qnumbers, reorder)
{
}
*/

//==============================================================================

/*
TEST(Qnumbers, extract)
{
}
*/

//==============================================================================

/*
TEST(Qnumbers, operator_parenthesis)
{
  Qnumbers q("n_z", 2);
  ASSERT_TRUE(q(0,0) == 0);
  ASSERT_TRUE(q(0,1) == 1);
}
*/

//==============================================================================

/*
TEST(Qnumbers, operator_parenthesis_)
{
}
*/

//==============================================================================

TEST(Qnumbers, calcBlocks)
{
  State state("examples/42Ca_deformed_2x9.msg.gz");
  Qnumbers &HOqn = state.basis.HOqn;
  arma::mat m = arma::zeros(HOqn.nb, HOqn.nb);
  UINT nbBlocks = state.basis.HOqn.calcBlocks({0});

  for (UINT i = 0; i < nbBlocks; i++)
  {
    for (UINT j = 0; j < nbBlocks; j++)
    {
      m.submat(HOqn.blocks[i].filter, HOqn.blocks[j].filter) = state.rho(NEUTRON).submat(HOqn.blocks[i].filter, HOqn.blocks[j].filter);
    }
  }

  ASSERT_NEAR(arma::norm(m - state.rho(NEUTRON), "inf"), 0.0, 1e-16);
}

//==============================================================================

/*
TEST(Qnumbers, setNames)
{
}
*/

//==============================================================================

/*
TEST(Qnumbers, append)
{
}
*/

//==============================================================================

TEST(Qnumbers, find)
{
  {
    // 1ct cylindrical basis
    Basis b(-1, 1.3, 1.2, 12, 20, 1.1);
    ASSERT_EQ(b.HOqn.nb, 465);

    for (UINT i = 0; i < b.HOqn.nb; i++)
    {
      INT a1 = b.HOqn(0, i);
      INT a2 = b.HOqn(1, i);
      INT a3 = b.HOqn(2, i);
      INT a4 = b.HOqn(3, i);
      INT a5 = b.HOqn(4, i);
      INT id = b.HOqn.find({a1, a2, a3, a4, a5});
      ASSERT_EQ(i, id);
    }
  }

  {
    // 2ct cylindrical basis
    Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);
    ASSERT_EQ(b.HOqn.nb, 930);

    for (UINT i = 0; i < b.HOqn.nb; i++)
    {
      INT a1 = b.HOqn(0, i);
      INT a2 = b.HOqn(1, i);
      INT a3 = b.HOqn(2, i);
      INT a4 = b.HOqn(3, i);
      INT a5 = b.HOqn(4, i);
      INT id = b.HOqn.find({a1, a2, a3, a4, a5});
      ASSERT_EQ(i, id);
    }
  }
}

//==============================================================================

/*
TEST(Qnumbers, checkOmegaSym)
{
}
*/

//==============================================================================


TEST(Qnumbers, clear)
{
  Qnumbers q("n_z", 2);
  IMAT mat(1, 2);
  for (INT i = 0; i < 2; i++)
  {
    mat(0, i) = i;
  }
  ASSERT_TRUE(q.nbq == 1);
  ASSERT_TRUE(q.nb == 2);
  ASSERT_TRUE(q.name[0] == "n_z");
  ASSERT_TRUE(q.mat == mat);

  q.clear();
  ASSERT_TRUE(q.nbq == 0);
  ASSERT_TRUE(q.nb == 0);
  ASSERT_TRUE(q.name.empty());
  ASSERT_TRUE(q.mat.empty());
} 
