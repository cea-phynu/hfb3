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
 *  \brief Test suite for the Multi and FMulti classes.
 */

#include "gtest/gtest.h"
#include "multi.h"

//==============================================================================

TEST(Multi, Multi)
{
  Multi<std::string> toto;

  toto(2, 3, 4) = "plop";
  toto(2, 3, 1) = "plip";

  ASSERT_EQ(toto(2, 3, 4), "plop");
  ASSERT_EQ(toto(2, 3, 1), "plip");

  ASSERT_EQ(toto.size(), 2);

  Multi<std::string> titi = toto;
  ASSERT_EQ(titi, toto);
}

//==============================================================================

TEST(Multi, operator_parenthesis)
{
  // Non const operator
  Multi<arma::vec> multiVec;

  arma::vec v1 = {1, 2, 3};
  arma::vec v2 = {9, 8, 7};

  multiVec(2, 3, 4) = v1;
  multiVec(5, 6, 7, 8) = v2;

  ASSERT_TRUE(arma::approx_equal(multiVec(2, 3, 4), v1, "absdiff", 1e-8));
  ASSERT_TRUE(arma::approx_equal(multiVec(5, 6, 7, 8), v2, "absdiff", 1e-8));
}

//==============================================================================

TEST(Multi, empty)
{
  Multi<arma::vec> multi;
  ASSERT_TRUE(multi.empty());

  multi(2,3) = arma::vec({1,2,3});
  ASSERT_FALSE(multi.empty());
}

//==============================================================================

TEST(Multi, clear)
{
  Multi<arma::vec> multi;
  multi(2,3) = arma::vec({1,2,3});
  ASSERT_FALSE(multi.empty());

  multi.clear();
  ASSERT_TRUE(multi.empty());
}

//==============================================================================

TEST(Multi, getKeys)
{
  Multi<arma::vec> multi;

  multi(2,3) = arma::vec({1,2,3});
  multi(2,3,4) = arma::vec({7,2});
  multi(7,8,9,10) = arma::vec({7,8,9,10});
  multi(4) = arma::vec({11});

  std::vector<std::vector<INT> > keys = {{4}, {7,8,9,10}, {2,3, 4}, {2,3}};
  ASSERT_NE(multi.getKeys(), keys);

  std::sort(keys.begin(), keys.end());

  ASSERT_EQ(multi.getKeys(), keys);
}

//==============================================================================

TEST(Multi, getValues)
{
  Multi<arma::vec> multi;

  arma::vec v1 = {1,2,3};
  arma::vec v2 = {7,2};
  arma::vec v3 = {7,8,9,10};
  arma::vec v4 = {11};

  multi(2,3) = v1;
  multi(2,3,4) = v2;
  multi(7,8,9,10) = v3;
  multi(4) = v4;

  int sum0 = arma::accu(v1) + arma::accu(v2) + arma::accu(v3) + arma::accu(v4);
  int sum1 = 0;

  for (arma::vec val : multi.getValues())
    sum1 += arma::accu(val);

  ASSERT_EQ(sum0, sum1);
}

//==============================================================================

TEST(Multi, keyToStr)
{
  Multi<std::string> multi;
  std::vector<INT> k = {2,3};

  ASSERT_EQ(multi.keyToStr(k), "[2,3]");
}

//==============================================================================

TEST(Multi, size)
{
  Multi<std::string> multi;

  multi(2,3) = "plip";
  multi(2,3,4) = "plop";
  multi(7,8,9,10) = "plouf";
  multi(4) = "paf";

  ASSERT_EQ(multi.size(), 4);
}

//==============================================================================

TEST(Multi, contains)
{
  Multi<std::string> multi;

  multi(2,3) = "plip";
  multi(2,3,4) = "plop";
  multi(7,8,9,10) = "plouf";
  multi(4) = "paf";

  ASSERT_TRUE(multi.contains(4));
  ASSERT_FALSE(multi.contains(2,4));
  ASSERT_FALSE(multi.contains(2,4,6,8,0));
}

//==============================================================================

TEST(Multi, operatorEquals)
{
  Multi<arma::vec> multi;
  multi(2,3) = arma::randu(50);
  multi(2,6,9,8) = arma::randu(100);

  Multi<arma::vec> multi2;
  multi2(4,9,4) = arma::randu(50);

  Multi<arma::vec> multi3;
  std::vector<INT> key = {2,6,9,8};
  multi3(key) = multi(key);
  multi3(2,3) = multi(2,3);

  // Same element
  ASSERT_TRUE(multi == multi);
  // Different keys and values
  ASSERT_FALSE(multi == multi2);
  // Same keys but different values
  ASSERT_FALSE(multi3 == multi2);
  // Same object but different order of the elements
  ASSERT_TRUE(multi3 == multi);

  multi3(7,2,0) = arma::vec();
  ASSERT_FALSE(multi3 == multi);
}

