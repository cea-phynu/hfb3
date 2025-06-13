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
 *  \brief Test suite for the Mesh class.
 */

#include "gtest/gtest.h"
#include "mesh.h"

//==============================================================================


TEST(Mesh, Mesh)
{
  Mesh m;
  Mesh m2 = Mesh::regular(0, 0, -20, 10, 0, 20, 41, 1, 161);
  ASSERT_EQ(m, m2);
  ASSERT_TRUE(m == m2);

  Mesh copy = m;
  ASSERT_EQ(m, copy);
}

//==============================================================================


TEST(Mesh, operator_equality)
{
  Mesh a;
  Mesh b;
  ASSERT_TRUE(a == b);
  ASSERT_EQ(a, b);

  Mesh c = Mesh::regular(0., 0., -10., 20., 0., 10., 41, 1, 81);
  ASSERT_FALSE(a == c);

  Mesh copy = a;
  ASSERT_TRUE(copy == a);
  ASSERT_EQ(copy, a);
}

//==============================================================================


TEST(Mesh, operator_inequality)
{
  Mesh a;
  Mesh b;
  ASSERT_FALSE(a != b);

  Mesh c = Mesh::regular(0., 0., -10., 20., 0., 10., 41, 1, 81);
  ASSERT_TRUE(a != c);

  Mesh copy = a;
  ASSERT_FALSE(copy != a);
}

//==============================================================================


TEST(Mesh, df3)
{
  Mesh m = Mesh::regular(-10, -10, -20, 10, 10, 20, 32, 32, 64);
  Mesh df3 = Mesh::df3();
  
  ASSERT_EQ(m, df3);
  ASSERT_TRUE(m == df3);
}

//==============================================================================

TEST(Mesh, regular)
{
  Mesh m = Mesh::regular(2.1, 3.2, 4.3, 8.4, 9.5, 10.6, 7, 8, 9);
  ASSERT_EQ(m.ax.nb, 7);
  ASSERT_EQ(m.ay.nb, 8);
  ASSERT_EQ(m.az.nb, 9);
}

//==============================================================================

/*
TEST(Mesh, info)
{
}
*/

//==============================================================================

/*
TEST(Mesh, gaussLaguerreHermite)
{
}
*/

