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
 *  \brief Test suite for the Tools class.
 */

#include "gtest/gtest.h"
#include "tools.h"

//==============================================================================

TEST(Tools, checkSymmetry)
{
  bool result = true;

  arma::mat asymMat = arma::randu(12, 12);
  result = Tools::checkSymmetry(asymMat, "this warning message should be printed");
  ASSERT_FALSE(result);

  arma::mat symMat = arma::symmatu(arma::randu(12, 12));
  result = Tools::checkSymmetry(symMat, "this warning message should not be printed");
  ASSERT_TRUE(result);
}

//==============================================================================

TEST(Tools, growIVec)
{
  IVEC m = {2, 3, 4};
  Tools::growIVec(m, 5);

  ASSERT_EQ(m(3), 5);
}

//==============================================================================

TEST(Tools, cubeTimesMat)
{
  arma::mat  m  = {{1, 2}, {1, 2}, {1, 2}};

  arma::cube c  = arma::zeros(3, 2, 2);
  c.slice(0)    = {{1, 2}, {3,  4}, { 5, 6.5}};
  c.slice(1)    = {{7, 8}, {9, 10}, {11 , 12}};

  arma::cube r  = Tools::cubeTimesMat(c, m);

  arma::cube r2 = arma::zeros(3, 2, 2);
  r2.slice(0)   = {{ 1,  4}, { 3,  8}, { 5, 13}};
  r2.slice(1)   = {{ 7, 16}, { 9, 20}, {11, 24}};

  ASSERT_TRUE(r == r2); // will be using operator==() defined in global.h
}

//==============================================================================

TEST(Tools, cubeTimesVec)
{
  arma::vec v   = {1, 2};

  arma::cube c  = arma::zeros(3, 2, 2);
  c.slice(0)    = {{1, 2}, {3, 4}, {5, 6.5}};
  c.slice(1)    = {{2, 4}, {6, 8}, {10, 13}};

  arma::cube r  = Tools::cubeTimesVec(c, v);

  arma::cube r2 = arma::zeros(3, 2, 2);
  r2.slice(0)   = {{1, 2}, {3, 4}, {5, 6.5}};
  r2.slice(1)   = {{4, 8}, {12, 16}, {20, 26}};

  ASSERT_TRUE(r == r2); // should use `operator==()` defined in global.h
}

//==============================================================================

TEST(Tools, cubeToMat)
{
  arma::cube c = arma::zeros(3, 2, 2);
  c.slice(0)   = {{1, 2}, {3, 4}, {5, 6.5}};
  c.slice(1)   = {{2, 4}, {6, 8}, {10, 13}};

  arma::mat r  = Tools::cubeToMat(c, 0, 0);
  arma::mat r2 = {{1, 2}, {2, 4}};

  ASSERT_TRUE(r == r2); // should use `operator==()` defined in global.h
                        //
  r  = Tools::cubeToMat(c, 1, 0);
  r2 = {{1, 2}, {3, 6}, {5, 10}};

  ASSERT_TRUE(r == r2); // should use `operator==()` defined in global.h

  r  = Tools::cubeToMat(c, 2, 1);
  r2 = {{2, 4}, {6, 8}, {10, 13}};

  ASSERT_TRUE(r == r2); // should use `operator==()` defined in global.h
}

//==============================================================================

TEST(Tools, matFromCol)
{
  arma::vec v  = {2, 3, 4.5};
  arma::mat r  = Tools::matFromCol(v, 3);
  arma::mat r2 = {{2, 2, 2}, {3, 3, 3}, {4.5, 4.5, 4.5}};

  ASSERT_NEAR(arma::norm(r - r2, "inf"), 0.0, 1e-20);
}

//==============================================================================

TEST(Tools, matFromRow)
{
  arma::rowvec v  = {2, 3, 4.5};
  arma::mat r  = Tools::matFromRow(v, 3);
  arma::mat r2 = {{2, 3, 4.5}, {2, 3, 4.5}, {2, 3, 4.5}};

  ASSERT_NEAR(arma::norm(r - r2, "inf"), 0.0, 1e-20);
}

//==============================================================================

TEST(Tools, matTimesVec)
{
  arma::vec v   = {1, 2};
  arma::mat m   = {{1, 2}, {3, 4}, {5, 6.5}};
  arma::cube r  = Tools::matTimesVec(m, v);

  arma::cube r2 = arma::zeros(3, 2, 2);
  r2.slice(0)   = {{1, 2}, {3, 4}, {5, 6.5}};
  r2.slice(1)   = {{2, 4}, {6, 8}, {10, 13}};

  ASSERT_NEAR(arma::abs(r - r2).max(), 0.0, 1e-20);
}

//==============================================================================

TEST(Tools, matToCsv)
{
  arma::mat m = {{1, 2.5}, {3, 4}, {5, 6}};
  std::string r = Tools::matToCsv(m);

  ASSERT_EQ(r, "0,1.000000e+00,2.500000e+00\n1,3.000000e+00,4.000000e+00\n2,5.000000e+00,6.000000e+00\n");
}

//==============================================================================

TEST(Tools, matToTxt)
{
  arma::mat m = {{1, 2.5}, {3, 4}, {5, 6}};
  std::string r = Tools::matToTxt(m);

  ASSERT_EQ(r, "1 2.5 \n3 4 \n5 6 \n");
}

//==============================================================================

TEST(Tools, trim)
{
  std::string s = "  test string with 2 leading spaces and 3 trailing spaces   ";
  std::string r = Tools::trim(s);
  ASSERT_EQ(r, "test string with 2 leading spaces and 3 trailing spaces");

  s = "lookMaNoSpaces";
  r = Tools::trim(s);
  ASSERT_EQ(r, "lookMaNoSpaces");
}

//==============================================================================

TEST(Tools, trim_b)
{
  std::string s = "  test string with 2 leading spaces and 3 trailing spaces   ";
  std::string r = Tools::trim_b(s);

  ASSERT_EQ(r, "test string with 2 leading spaces and 3 trailing spaces   ");

  s = "lookMaNoSpaces";
  r = Tools::trim_b(s);
  ASSERT_EQ(r, "lookMaNoSpaces");
}

//==============================================================================

TEST(Tools, trim_e)
{
  std::string s = "  test string with 2 leading spaces and 3 trailing spaces   ";
  std::string r = Tools::trim_e(s);

  ASSERT_EQ(r, "  test string with 2 leading spaces and 3 trailing spaces");

  s = "lookMaNoSpaces";
  r = Tools::trim_e(s);
  ASSERT_EQ(r, "lookMaNoSpaces");
}

//==============================================================================

TEST(Tools, strIsospin)
{
  ASSERT_EQ(Tools::strIsospin(NEUTRON), "neutron");
  ASSERT_EQ(Tools::strIsospin(PROTON ), "proton" );
  ASSERT_EQ(Tools::strIsospin(TOTAL  ), "total"  );
  ASSERT_EQ(Tools::strIsospin(8      ), "unknown isospin !?");
}

//==============================================================================

TEST(Tools, stringLen)
{
  ASSERT_EQ(Tools::stringLen("123456789012345"), 15);
}

