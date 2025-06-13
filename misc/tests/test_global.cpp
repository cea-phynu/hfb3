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
 *  \brief Test suite for the Global class.
 */

#include "gtest/gtest.h"
#include "global.h"

//==============================================================================

TEST(Global, operatorEquals)
{
  {
    // uvec
    UVEC v1 = UVEC({1,2,3});
    UVEC v2 = UVEC({1,3,2});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    UVEC v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // vec
    arma::vec v1 = arma::vec({1,2,3});
    arma::vec v2 = arma::vec({1,3,2});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    arma::vec v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // ivec
    IVEC v1 = IVEC({1,2,3});
    IVEC v2 = IVEC({1,3,2});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    IVEC v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // umat
    arma::umat v1 = arma::umat({{1,2,3},{4,5,6}});
    arma::umat v2 = arma::umat({{1,3,2},{5,4,6}});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    arma::umat v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // mat
    arma::mat v1 = arma::mat({{1,2,3},{4,5,6}});
    arma::mat v2 = arma::mat({{1,3,2},{5,4,6}});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    arma::mat v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // imat
    IMAT v1 = IMAT({{1,2,3},{4,5,6}});
    IMAT v2 = IMAT({{1,3,2},{5,4,6}});

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    IMAT v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // ucube
    UCUBE v1 = UCUBE(3,3,3);
    UCUBE v2 = UCUBE(3,3,3);
    v2(2,0,1) = 4;

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    UCUBE v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // cube
    arma::cube v1 = arma::cube(3,3,3);
    arma::cube v2 = arma::cube(3,3,3);
    v2(2,0,1) = 4;

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    arma::cube v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }

  {
    // icube
    ICUBE v1 = ICUBE(3,3,3);
    ICUBE v2 = ICUBE(3,3,3);
    v2(2,0,1) = 4;

    ASSERT_TRUE(v1 == v1);
    ASSERT_FALSE(v1 == v2);

    ICUBE v1_copy = v1;
    ASSERT_TRUE(v1 == v1_copy);
  }
}
