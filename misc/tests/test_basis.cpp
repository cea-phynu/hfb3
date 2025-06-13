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
 *  \brief Test suite for the Basis class.
 */

#include "gtest/gtest.h"
#include "basis.h"
#include "discrete.h"
#include "state.h"
#include "maxima_targets.h"
#include "tools.h"

//==============================================================================

TEST(Basis, Basis_)
{
  std::string filename = "misc/data/test_2ct_256Fm.hfb3";
  DataTree dataTree(filename);
  Basis basis(dataTree);
  ASSERT_EQ(basis.nOscil, 11);
  ASSERT_EQ(basis.dMax, 2);
  ASSERT_EQ(basis.n_zMaxImposed, 24);
  ASSERT_NEAR(basis.d_0, 5.5, 1e-15);
  ASSERT_NEAR(basis.b_r, 2.1, 1e-15);
  ASSERT_NEAR(basis.b_z, 2.4, 1e-15);
  ASSERT_NEAR(basis.g_q, 1.3, 1e-15);
}

//==============================================================================

TEST(Basis, Basis__)
{
  std::string filename = "misc/data/test_2ct_256Fm.hfb3";
  DataTree dataTree(filename);
  Basis basis(dataTree);
  ASSERT_EQ(basis.nOscil, 11);
  ASSERT_EQ(basis.dMax, 2);
  ASSERT_EQ(basis.n_zMaxImposed, 24);
  ASSERT_NEAR(basis.d_0, 5.5, 1e-15);
  ASSERT_NEAR(basis.b_r, 2.1, 1e-15);
  ASSERT_NEAR(basis.b_z, 2.4, 1e-15);
  ASSERT_NEAR(basis.g_q, 1.3, 1e-15);
}

//==============================================================================

TEST(Basis, cylindrical)
{
  Basis basis(-1.0, 2.6108016647931507, 2.8296839564912175, 14, 24, 1.3);
  ASSERT_EQ(basis.nOscil, 14);
  ASSERT_EQ(basis.dMax, 1);
  ASSERT_EQ(basis.n_zMaxImposed, 24);
  ASSERT_NEAR(basis.d_0, 0.0, 1e-15);
  ASSERT_NEAR(basis.b_r, 2.6108016647931507, 1e-15);
  ASSERT_NEAR(basis.b_z, 2.8296839564912175, 1e-15);
  ASSERT_NEAR(basis.g_q, 1.3, 1e-15);
}

//==============================================================================

TEST(Basis, berger2ct)
{
  Basis basis = Basis::fromBerger2ct(14, 1.3, 12.0, 1.2, -1.0, 24);
  ASSERT_NEAR(basis.d_0, 0.0, 1e-15);
  ASSERT_NEAR(basis.b_r, 1.9163431074141117, 1e-15);
  ASSERT_NEAR(basis.b_z, 1.7493739130837576, 1e-15);
  ASSERT_NEAR(basis.g_q, 1.3, 1e-15);
  ASSERT_EQ(basis.nOscil, 14);
  ASSERT_EQ(basis.dMax, 1);
  ASSERT_EQ(basis.n_zMaxImposed, 24);
}

//==============================================================================

TEST(Basis, cylTruncate)
{
  Basis b1(-1, 1.3, 1.2, 12, 20, 1.1);
  ASSERT_EQ(b1.HOqn.nb, 465u);
  ASSERT_EQ(b1.HOqn(0, 301), 4);
  ASSERT_EQ(b1.HOqn(1, 301), 1);
  ASSERT_EQ(b1.HOqn(2, 301), 4);
  ASSERT_EQ(b1.HOqn(3, 301), 0);
  ASSERT_EQ(b1.HOqn(4, 301), 0);

  Basis b2(12.3, 1.3, 1.2, 12, 20, 1.1);
  ASSERT_EQ(b2.HOqn.nb, 930u);
  ASSERT_EQ(b2.HOqn(0, 831), 7);
  ASSERT_EQ(b2.HOqn(1, 831), 1);
  ASSERT_EQ(b2.HOqn(2, 831), 3);
  ASSERT_EQ(b2.HOqn(3, 831), 0);
  ASSERT_EQ(b2.HOqn(4, 831), 1);
}

//==============================================================================

TEST(Basis, hermite)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-09;

  double x;
  x = -2.238;
  ASSERT_NEAR(b.hermite(0, x), 1, precision);
  ASSERT_NEAR(b.hermite(1, x), 2 * x, precision);
  ASSERT_NEAR(b.hermite(2, x), 4 * pow(x, 2) - 2, precision);
  ASSERT_NEAR(b.hermite(3, x), 8 * pow(x, 3) - 12 * x, precision);
  ASSERT_NEAR(b.hermite(7, x), 128 * pow(x, 7) - 1344 * pow(x, 5) + 3360 * pow(x, 3) - 1680 * x, precision);
  x = 4.123;
  ASSERT_NEAR(b.hermite(0, x), 1, precision);
  ASSERT_NEAR(b.hermite(1, x), 2 * x, precision);
  ASSERT_NEAR(b.hermite(2, x), 4 * pow(x, 2) - 2, precision);
  ASSERT_NEAR(b.hermite(3, x), 8 * pow(x, 3) - 12 * x, precision);
  ASSERT_NEAR(b.hermite(7, x), 128 * pow(x, 7) - 1344 * pow(x, 5) + 3360 * pow(x, 3) - 1680 * x, precision);
}

//==============================================================================

TEST(Basis, laguerre)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-14;

  double x;
  INT m;
  m = 2;
  x = -2.238;
  ASSERT_NEAR(b.laguerre(m, 0, x), 1, precision);
  ASSERT_NEAR(b.laguerre(m, 1, x), 1 + m - x, precision);
  ASSERT_NEAR(b.laguerre(m, 2, x), 0.5 * pow(x, 2) - (m + 2) * x + 0.5 * (m + 1) * (m + 2), precision);
  m = -3;
  x = 5.567;
  ASSERT_NEAR(b.laguerre(m, 0, x), 1, precision);
  ASSERT_NEAR(b.laguerre(m, 1, x), 1 + m - x, precision);
  ASSERT_NEAR(b.laguerre(m, 2, x), 0.5 * pow(x, 2) - (m + 2) * x + 0.5 * (m + 1) * (m + 2), precision);
}

//============================================================================

TEST(Basis, legendre)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-10;

  double x;
  x = -2.238;
  ASSERT_NEAR(b.legendre(0, x), 1, precision);
  ASSERT_NEAR(b.legendre(1, x), x, precision);
  ASSERT_NEAR(b.legendre(2, x), (3 * pow(x, 2) - 1) / 2, precision);
  ASSERT_NEAR(b.legendre(3, x), (5 * pow(x, 3) - 3 * x) / 2, precision);
  ASSERT_NEAR(b.legendre(7, x), (429 * pow(x, 7) - 693 * pow(x, 5) + 315 * pow(x, 3) - 35 * x) / 16, precision);
  x = 4.123;
  ASSERT_NEAR(b.legendre(0, x), 1, precision);
  ASSERT_NEAR(b.legendre(1, x), x, precision);
  ASSERT_NEAR(b.legendre(2, x), (3 * pow(x, 2) - 1) / 2, precision);
  ASSERT_NEAR(b.legendre(3, x), (5 * pow(x, 3) - 3 * x) / 2, precision);
  ASSERT_NEAR(b.legendre(7, x), (429 * pow(x, 7) - 693 * pow(x, 5) + 315 * pow(x, 3) - 35 * x) / 16, precision);
}

//==============================================================================

TEST(Basis, hermite_)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  arma::vec x = {-2.238, 4.123, 2.978, 0.234};

  arma::vec res0 = b.hermite(0, x);
  arma::vec res1 = b.hermite(1, x);
  arma::vec res2 = b.hermite(2, x);
  arma::vec res7 = b.hermite(7, x);

  for (UINT i = 0; i < x.n_elem; i++)
  {
    ASSERT_NEAR(res0(i), b.hermite(0, x(i)), precision);
    ASSERT_NEAR(res1(i), b.hermite(1, x(i)), precision);
    ASSERT_NEAR(res2(i), b.hermite(2, x(i)), precision);
    ASSERT_NEAR(res7(i), b.hermite(7, x(i)), precision);
  }
}

//==============================================================================

TEST(Basis, laguerre_)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  arma::vec x = {-2.238, 4.123, 2.978, 0.234};

  arma::vec res0 = b.laguerre(3, 0, x);
  arma::vec res1 = b.laguerre(3, 1, x);
  arma::vec res2 = b.laguerre(3, 2, x);
  arma::vec res7 = b.laguerre(3, 7, x);

  for (UINT i = 0; i < x.n_elem; i++)
  {
    ASSERT_NEAR(res0(i), b.laguerre(3, 0, x(i)), precision);
    ASSERT_NEAR(res1(i), b.laguerre(3, 1, x(i)), precision);
    ASSERT_NEAR(res2(i), b.laguerre(3, 2, x(i)), precision);
    ASSERT_NEAR(res7(i), b.laguerre(3, 7, x(i)), precision);
  }
}

//==============================================================================

TEST(Basis, legendre_)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  arma::vec x = {-2.238, 4.123, 2.978, 0.234};

  arma::vec res0 = b.legendre(0, x);
  arma::vec res1 = b.legendre(1, x);
  arma::vec res2 = b.legendre(2, x);
  arma::vec res7 = b.legendre(7, x);

  for (UINT i = 0; i < x.n_elem; i++)
  {
    ASSERT_NEAR(res0(i), b.legendre(0, x(i)), precision);
    ASSERT_NEAR(res1(i), b.legendre(1, x(i)), precision);
    ASSERT_NEAR(res2(i), b.legendre(2, x(i)), precision);
    ASSERT_NEAR(res7(i), b.legendre(7, x(i)), precision);
  }
}

//==============================================================================

TEST(Basis, zPart)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  double z = -2.238;

  ASSERT_NEAR(b.zPartScalar(z,  0, 0), BASIS_ZPART_TARGET00, precision);
  ASSERT_NEAR(b.zPartScalar(z,  1, 1), BASIS_ZPART_TARGET01, precision);
  ASSERT_NEAR(b.zPartScalar(z,  2, 0), BASIS_ZPART_TARGET02, precision);
  ASSERT_NEAR(b.zPartScalar(z,  3, 1), BASIS_ZPART_TARGET03, precision);
  ASSERT_NEAR(b.zPartScalar(z,  4, 0), BASIS_ZPART_TARGET04, precision);
  ASSERT_NEAR(b.zPartScalar(z,  5, 1), BASIS_ZPART_TARGET05, precision);
  ASSERT_NEAR(b.zPartScalar(z,  6, 0), BASIS_ZPART_TARGET06, precision);
  ASSERT_NEAR(b.zPartScalar(z,  7, 1), BASIS_ZPART_TARGET07, precision);
  ASSERT_NEAR(b.zPartScalar(z,  8, 0), BASIS_ZPART_TARGET08, precision);
  ASSERT_NEAR(b.zPartScalar(z,  9, 1), BASIS_ZPART_TARGET09, precision);
  ASSERT_NEAR(b.zPartScalar(z, 10, 0), BASIS_ZPART_TARGET10, precision);
  ASSERT_NEAR(b.zPartScalar(z, 11, 1), BASIS_ZPART_TARGET11, precision);
  ASSERT_NEAR(b.zPartScalar(z, 12, 0), BASIS_ZPART_TARGET12, precision);
  ASSERT_NEAR(b.zPartScalar(z, 13, 1), BASIS_ZPART_TARGET13, precision);
  ASSERT_NEAR(b.zPartScalar(z, 14, 0), BASIS_ZPART_TARGET14, precision);
  ASSERT_NEAR(b.zPartScalar(z, 15, 1), BASIS_ZPART_TARGET15, precision);
  ASSERT_NEAR(b.zPartScalar(z, 16, 0), BASIS_ZPART_TARGET16, precision);
  ASSERT_NEAR(b.zPartScalar(z, 17, 1), BASIS_ZPART_TARGET17, precision);
  ASSERT_NEAR(b.zPartScalar(z, 18, 0), BASIS_ZPART_TARGET18, precision);
  ASSERT_NEAR(b.zPartScalar(z, 19, 1), BASIS_ZPART_TARGET19, precision);
  ASSERT_NEAR(b.zPartScalar(z, 20, 0), BASIS_ZPART_TARGET20, precision);
}

//==============================================================================

TEST(Basis, zPart_big)
{
  Basis b(12.3, 1.3, 5.2, 24, 24, 1.1);

  double precision = 1e-15;

  double z = -2.238;

  ASSERT_NEAR(b.zPartScalar(z,  0, 0), BASIS_ZPART_BIG_TARGET00, precision);
  ASSERT_NEAR(b.zPartScalar(z,  1, 1), BASIS_ZPART_BIG_TARGET01, precision);
  ASSERT_NEAR(b.zPartScalar(z,  2, 0), BASIS_ZPART_BIG_TARGET02, precision);
  ASSERT_NEAR(b.zPartScalar(z,  3, 1), BASIS_ZPART_BIG_TARGET03, precision);
  ASSERT_NEAR(b.zPartScalar(z,  4, 0), BASIS_ZPART_BIG_TARGET04, precision);
  ASSERT_NEAR(b.zPartScalar(z,  5, 1), BASIS_ZPART_BIG_TARGET05, precision);
  ASSERT_NEAR(b.zPartScalar(z,  6, 0), BASIS_ZPART_BIG_TARGET06, precision);
  ASSERT_NEAR(b.zPartScalar(z,  7, 1), BASIS_ZPART_BIG_TARGET07, precision);
  ASSERT_NEAR(b.zPartScalar(z,  8, 0), BASIS_ZPART_BIG_TARGET08, precision);
  ASSERT_NEAR(b.zPartScalar(z,  9, 1), BASIS_ZPART_BIG_TARGET09, precision);
  ASSERT_NEAR(b.zPartScalar(z, 10, 0), BASIS_ZPART_BIG_TARGET10, precision);
  ASSERT_NEAR(b.zPartScalar(z, 11, 1), BASIS_ZPART_BIG_TARGET11, precision);
  ASSERT_NEAR(b.zPartScalar(z, 12, 0), BASIS_ZPART_BIG_TARGET12, precision);
  ASSERT_NEAR(b.zPartScalar(z, 13, 1), BASIS_ZPART_BIG_TARGET13, precision);
  ASSERT_NEAR(b.zPartScalar(z, 14, 0), BASIS_ZPART_BIG_TARGET14, precision);
  ASSERT_NEAR(b.zPartScalar(z, 15, 1), BASIS_ZPART_BIG_TARGET15, precision);
  ASSERT_NEAR(b.zPartScalar(z, 16, 0), BASIS_ZPART_BIG_TARGET16, precision);
  ASSERT_NEAR(b.zPartScalar(z, 17, 1), BASIS_ZPART_BIG_TARGET17, precision);
  ASSERT_NEAR(b.zPartScalar(z, 18, 0), BASIS_ZPART_BIG_TARGET18, precision);
  ASSERT_NEAR(b.zPartScalar(z, 19, 1), BASIS_ZPART_BIG_TARGET19, precision);
  ASSERT_NEAR(b.zPartScalar(z, 20, 0), BASIS_ZPART_BIG_TARGET20, precision);
}

//==============================================================================

TEST(Basis, rPart)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  double r = -2.238;

  ASSERT_NEAR(b.rPartScalar(r,  0,  2), BASIS_RPART_TARGET00, precision);
  ASSERT_NEAR(b.rPartScalar(r,  1,  1), BASIS_RPART_TARGET01, precision);
  ASSERT_NEAR(b.rPartScalar(r,  2,  3), BASIS_RPART_TARGET02, precision);
  ASSERT_NEAR(b.rPartScalar(r,  3,  2), BASIS_RPART_TARGET03, precision);
  ASSERT_NEAR(b.rPartScalar(r,  4,  1), BASIS_RPART_TARGET04, precision);
  ASSERT_NEAR(b.rPartScalar(r,  5,  4), BASIS_RPART_TARGET05, precision);
  ASSERT_NEAR(b.rPartScalar(r,  6,  3), BASIS_RPART_TARGET06, precision);
  ASSERT_NEAR(b.rPartScalar(r,  7,  4), BASIS_RPART_TARGET07, precision);
  ASSERT_NEAR(b.rPartScalar(r,  8,  2), BASIS_RPART_TARGET08, precision);
  ASSERT_NEAR(b.rPartScalar(r,  9,  5), BASIS_RPART_TARGET09, precision);
  ASSERT_NEAR(b.rPartScalar(r, 10,  2), BASIS_RPART_TARGET10, precision);
  ASSERT_NEAR(b.rPartScalar(r, 11, 11), BASIS_RPART_TARGET11, precision);
  ASSERT_NEAR(b.rPartScalar(r, 12,  3), BASIS_RPART_TARGET12, precision);
  ASSERT_NEAR(b.rPartScalar(r, 13, 12), BASIS_RPART_TARGET13, precision);
  ASSERT_NEAR(b.rPartScalar(r, 14,  1), BASIS_RPART_TARGET14, precision);
  ASSERT_NEAR(b.rPartScalar(r, 15, 14), BASIS_RPART_TARGET15, precision);
  ASSERT_NEAR(b.rPartScalar(r, 16,  3), BASIS_RPART_TARGET16, precision);
  ASSERT_NEAR(b.rPartScalar(r, 17, 14), BASIS_RPART_TARGET17, precision);
  ASSERT_NEAR(b.rPartScalar(r, 18,  2), BASIS_RPART_TARGET18, precision);
  ASSERT_NEAR(b.rPartScalar(r, 19, 15), BASIS_RPART_TARGET19, precision);
  ASSERT_NEAR(b.rPartScalar(r, 20,  5), BASIS_RPART_TARGET20, precision);
}

//==============================================================================

TEST(Basis, rPart_big)
{
  Basis b(12.3, 5.3, 1.2, 24, 24, 1.1);

  double precision = 1e-15;

  double r = -2.238;

  ASSERT_NEAR(b.rPartScalar(r,  0,  2), BASIS_RPART_BIG_TARGET00, precision);
  ASSERT_NEAR(b.rPartScalar(r,  1,  1), BASIS_RPART_BIG_TARGET01, precision);
  ASSERT_NEAR(b.rPartScalar(r,  2,  3), BASIS_RPART_BIG_TARGET02, precision);
  ASSERT_NEAR(b.rPartScalar(r,  3,  2), BASIS_RPART_BIG_TARGET03, precision);
  ASSERT_NEAR(b.rPartScalar(r,  4,  1), BASIS_RPART_BIG_TARGET04, precision);
  ASSERT_NEAR(b.rPartScalar(r,  5,  4), BASIS_RPART_BIG_TARGET05, precision);
  ASSERT_NEAR(b.rPartScalar(r,  6,  3), BASIS_RPART_BIG_TARGET06, precision);
  ASSERT_NEAR(b.rPartScalar(r,  7,  4), BASIS_RPART_BIG_TARGET07, precision);
  ASSERT_NEAR(b.rPartScalar(r,  8,  2), BASIS_RPART_BIG_TARGET08, precision);
  ASSERT_NEAR(b.rPartScalar(r,  9,  5), BASIS_RPART_BIG_TARGET09, precision);
  ASSERT_NEAR(b.rPartScalar(r, 10,  2), BASIS_RPART_BIG_TARGET10, precision);
  ASSERT_NEAR(b.rPartScalar(r, 11, 11), BASIS_RPART_BIG_TARGET11, precision);
  ASSERT_NEAR(b.rPartScalar(r, 12,  3), BASIS_RPART_BIG_TARGET12, precision);
  ASSERT_NEAR(b.rPartScalar(r, 13, 12), BASIS_RPART_BIG_TARGET13, precision);
  ASSERT_NEAR(b.rPartScalar(r, 14,  1), BASIS_RPART_BIG_TARGET14, precision);
  ASSERT_NEAR(b.rPartScalar(r, 15, 14), BASIS_RPART_BIG_TARGET15, precision);
  ASSERT_NEAR(b.rPartScalar(r, 16,  3), BASIS_RPART_BIG_TARGET16, precision);
  ASSERT_NEAR(b.rPartScalar(r, 17, 14), BASIS_RPART_BIG_TARGET17, precision);
  ASSERT_NEAR(b.rPartScalar(r, 18,  2), BASIS_RPART_BIG_TARGET18, precision);
  ASSERT_NEAR(b.rPartScalar(r, 19, 15), BASIS_RPART_BIG_TARGET19, precision);
  ASSERT_NEAR(b.rPartScalar(r, 20,  5), BASIS_RPART_BIG_TARGET20, precision);
}

//==============================================================================

TEST(Basis, zPart_)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-15;

  arma::vec z = {-2.238, 4.123, 2.978, 0.234};

  arma::vec res0 = b.zPart(z, 0, 0);
  arma::vec res1 = b.zPart(z, 1, 1);
  arma::vec res2 = b.zPart(z, 2, 0);
  arma::vec res3 = b.zPart(z, 3, 1);
  arma::vec res4 = b.zPart(z, 4, 0);
  arma::vec res5 = b.zPart(z, 5, 1);
  arma::vec res6 = b.zPart(z, 6, 0);
  arma::vec res7 = b.zPart(z, 7, 1);
  arma::vec res8 = b.zPart(z, 8, 0);
  arma::vec res9 = b.zPart(z, 9, 1);

  for (UINT i = 0; i < z.n_elem; i++)
  {
    ASSERT_NEAR(res0(i), b.zPartScalar(z(i), 0, 0), precision);
    ASSERT_NEAR(res1(i), b.zPartScalar(z(i), 1, 1), precision);
    ASSERT_NEAR(res2(i), b.zPartScalar(z(i), 2, 0), precision);
    ASSERT_NEAR(res3(i), b.zPartScalar(z(i), 3, 1), precision);
    ASSERT_NEAR(res4(i), b.zPartScalar(z(i), 4, 0), precision);
    ASSERT_NEAR(res5(i), b.zPartScalar(z(i), 5, 1), precision);
    ASSERT_NEAR(res6(i), b.zPartScalar(z(i), 6, 0), precision);
    ASSERT_NEAR(res7(i), b.zPartScalar(z(i), 7, 1), precision);
    ASSERT_NEAR(res8(i), b.zPartScalar(z(i), 8, 0), precision);
    ASSERT_NEAR(res9(i), b.zPartScalar(z(i), 9, 1), precision);
  }
}

//==============================================================================

TEST(Basis, rPart_)
{
  Basis b(12.3, 1.3, 1.2, 12, 20, 1.1);

  double precision = 1e-13;

  arma::vec r = {2.238, 4.123, 2.978, 0.234};

  arma::vec res0 = b.rPart(r, 0, 2);
  arma::vec res1 = b.rPart(r, 1, 1);
  arma::vec res2 = b.rPart(r, 2, 3);
  arma::vec res3 = b.rPart(r, 3, 2);
  arma::vec res4 = b.rPart(r, 4, 1);
  arma::vec res5 = b.rPart(r, 5, 4);
  arma::vec res6 = b.rPart(r, 6, 3);
  arma::vec res7 = b.rPart(r, 7, 4);
  arma::vec res8 = b.rPart(r, 8, 2);
  arma::vec res9 = b.rPart(r, 9, 5);

  for (UINT i = 0; i < r.n_elem; i++)
  {
    ASSERT_NEAR(res0(i), b.rPartScalar(r(i), 0, 2), precision);
    ASSERT_NEAR(res1(i), b.rPartScalar(r(i), 1, 1), precision);
    ASSERT_NEAR(res2(i), b.rPartScalar(r(i), 2, 3), precision);
    ASSERT_NEAR(res3(i), b.rPartScalar(r(i), 3, 2), precision);
    ASSERT_NEAR(res4(i), b.rPartScalar(r(i), 4, 1), precision);
    ASSERT_NEAR(res5(i), b.rPartScalar(r(i), 5, 4), precision);
    ASSERT_NEAR(res6(i), b.rPartScalar(r(i), 6, 3), precision);
    ASSERT_NEAR(res7(i), b.rPartScalar(r(i), 7, 4), precision);
    ASSERT_NEAR(res8(i), b.rPartScalar(r(i), 8, 2), precision);
    ASSERT_NEAR(res9(i), b.rPartScalar(r(i), 9, 5), precision);
  }
}

//==============================================================================

TEST(Basis, calcWDN)
{
  State state("examples/42Ca_deformed_2x9.msg.gz");
  Basis &basis = state.basis;
  basis.calcWDN();

  for (INT m = 0; m < basis.mMax; m++)
  {
    INT n_zMax = basis.n_zMax(m, 0);
    arma::mat WDN = basis.W(m) * basis.N(m);
    arma::mat R = arma::eye(n_zMax * 2, n_zMax * 2);
    R.submat(0, n_zMax, n_zMax - 1, n_zMax * 2 - 1) = basis.Tab.submat(0, 0, n_zMax - 1, n_zMax - 1);
    R.submat(n_zMax, 0, n_zMax * 2 - 1, n_zMax - 1) = basis.Tab.submat(0, 0, n_zMax - 1, n_zMax - 1).t();
    arma::mat result = WDN.t() * R * WDN - arma::eye(WDN.n_cols, WDN.n_cols);
    // Error for each m block
    ASSERT_NEAR(arma::norm(result, "inf"), 0.0, 1e-10);
  }

  arma::mat M  = basis.HOtoOR.t();
  arma::mat Mi = basis.ORtoHO.t();

  // Test 001
  arma::mat result1 = basis.S - basis.S.t();
  ASSERT_NEAR(arma::norm(result1, "inf"), 0.0, 1e-20);
  // Test 002
  arma::mat result2 = Mi * M - arma::eye(basis.ORqn.nb, basis.ORqn.nb);
  ASSERT_NEAR(arma::norm(result2, "inf"), 0.0, 1e-10);
  // Test 003
  arma::mat result3 = M.t() * basis.S * M - arma::eye(basis.ORqn.nb, basis.ORqn.nb);
  ASSERT_NEAR(arma::norm(result3, "inf"), 0.0, 1e-10);
  // Test 004
  arma::mat result4 = basis.S * M - Mi.t();
  ASSERT_NEAR(arma::norm(result4, "inf"), 0.0, 1e-10);

  /*
  // Test 005
  arma::mat result5 = basis.S * M * M.t() - arma::eye(arma::size(basis.S));
  ASSERT_NEAR(arma::norm(result5, "inf"), 0.0, 1e-10);
  // Test 006
  arma::mat result6 = M * M.t() * basis.S - arma::eye(arma::size(basis.S));
  ASSERT_NEAR(arma::norm(result6, "inf"), 0.0, 1e-10);
  // Test 007
  arma::mat result7 = M * Mi - arma::eye(basis.HOqn.nb, basis.HOqn.nb);
  ASSERT_NEAR(arma::norm(result7, "inf"), 0.0, 1e-10);
  // Test 008
  arma::mat result8 = Mi - M.t() * (M * M.t()).i();
  ASSERT_NEAR(arma::norm(result8, "inf"), 0.0, 1e-10);
  */
  // Test 009
  arma::mat result9 = Mi - (M.t() * M).i() * M.t();
  ASSERT_NEAR(arma::norm(result9, "inf"), 0.0, 1e-10);
  // Neutron number
  arma::mat rn = 2 * state.rho(NEUTRON); // coeff 2 : negative AND positive time-reversal states
  result2 = basis.S * rn;
  ASSERT_NEAR(arma::trace(result2), 22.0, 1e-4);
  // Proton number
  arma::mat rp = 2 * state.rho(PROTON );
  result2 = basis.S * rp;
  ASSERT_NEAR(arma::trace(result2), 20.0, 1e-4);
  // Total
  arma::mat rt = rn + rp;
  result2 = basis.S * rt;
  ASSERT_NEAR(arma::trace(result2), 42.0, 1e-4);


  // Test the values of rho and kappa in the biorthogonal basis.
  // basis.calcWDN();
  // arma::mat rho(NEUTRON)HO = basis.S * state.rho(NEUTRON) * basis.S.t();
  // arma::mat rho(NEUTRON)HF = state.HOtoHFn * rho(NEUTRON)HO * state.HOtoHFn.t();
  // arma::vec ui2b = arma::sort(state.vecOcn, "descend");
  // arma::vec ui2 = arma::sort(arma::diagvec(rho(NEUTRON)HF), "descend");
  // ASSERT_NEAR(arma::norm(ui2 - ui2b, "inf"), 0., 1.0e-07);
  // arma::mat kappanHO = basis.S * state.kappan * basis.S.t();
  // arma::mat kappanHF = state.HOtoHFn * kappanHO * state.HOtoHFn.t();
  // arma::vec vi2b = arma::sort(arma::ones<arma::vec>(state.vecOcn.n_elem) - state.vecOcn, "ascend");
  // arma::vec uivib = arma::sort(arma::sqrt(ui2b % vi2b), "ascend");
  // arma::vec uivi = arma::sort(arma::diagvec(kappanHF), "ascend");
  // ASSERT_NEAR(arma::norm(uivi - uivib, "inf"), 0., 1.0e-09);
}

//==============================================================================

TEST(Basis, calcHharmo2ct)
{
  // Test the calculation of the "harmonic 2-center potential" matrix elements for a cylindrical 2-center basis.
  Basis basis("examples/42Ca_deformed_2x9.msg.gz");
  basis.calcHharmo2ct();

  ASSERT_NEAR(basis.Recouv( 0 +  9 * 2, 0 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET00, 1e-14);
  ASSERT_NEAR(basis.Recouv( 0 +  9 * 2, 1 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET01, 1e-14);
  ASSERT_NEAR(basis.Recouv( 1 +  9 * 2, 0 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET02, 1e-14);
  ASSERT_NEAR(basis.Recouv( 1 +  9 * 2, 1 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET03, 1e-14);
  ASSERT_NEAR(basis.Recouvz(0 + 11 * 2, 0 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET04, 1e-14);
  ASSERT_NEAR(basis.Recouvz(0 + 11 * 2, 1 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET05, 1e-13);
  ASSERT_NEAR(basis.Recouvz(1 + 11 * 2, 0 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET06, 1e-13);
  ASSERT_NEAR(basis.Recouvz(1 + 11 * 2, 1 + 10 * 2), BASIS_CALCHHARMO2CT_TARGET07, 1e-14);

  // (2, 0, 5, 0, 0) x (2, 0, 3, 0, 0)
  ASSERT_NEAR(basis.H(2)(5 + 0 * basis.n_zMax(2), 3 + 0 * basis.n_zMax(2)), BASIS_CALCHHARMO2CT_TARGET08, 1e-14);
  // (2, 0, 5, 0, 0) x (2, 0, 3, 1, 0)
  ASSERT_NEAR(basis.H(2)(5 + 0 * basis.n_zMax(2), 3 + 1 * basis.n_zMax(2)), BASIS_CALCHHARMO2CT_TARGET09, 1e-14);
  // (2, 0, 5, 1, 0) x (2, 0, 3, 0, 0)
  ASSERT_NEAR(basis.H(2)(5 + 1 * basis.n_zMax(2), 3 + 0 * basis.n_zMax(2)), BASIS_CALCHHARMO2CT_TARGET10, 1e-14);
  // (2, 0, 5, 1, 0) x (2, 0, 3, 1, 0)
  ASSERT_NEAR(basis.H(2)(5 + 1 * basis.n_zMax(2), 3 + 1 * basis.n_zMax(2)), BASIS_CALCHHARMO2CT_TARGET11, 1e-14);
  // (3, 0, 3, 0, 0) x (3, 0, 1, 0, 0)
  ASSERT_NEAR(basis.H(3)(3 + 0 * basis.n_zMax(3), 1 + 0 * basis.n_zMax(3)), BASIS_CALCHHARMO2CT_TARGET12, 1e-14);
  // (3, 0, 3, 0, 0) x (3, 0, 1, 1, 0)
  ASSERT_NEAR(basis.H(3)(3 + 0 * basis.n_zMax(3), 1 + 1 * basis.n_zMax(3)), BASIS_CALCHHARMO2CT_TARGET13, 1e-13);
  // (3, 0, 3, 1, 0) x (3, 0, 1, 0, 0)
  ASSERT_NEAR(basis.H(3)(3 + 1 * basis.n_zMax(3), 1 + 0 * basis.n_zMax(3)), BASIS_CALCHHARMO2CT_TARGET14, 1e-13);
  // (3, 0, 3, 1, 0) x (3, 0, 1, 1, 0)
  ASSERT_NEAR(basis.H(3)(3 + 1 * basis.n_zMax(3), 1 + 1 * basis.n_zMax(3)), BASIS_CALCHHARMO2CT_TARGET15, 1e-14);
}

//==============================================================================

TEST(Basis, getTransformationMatrix)
{
  {
    // 1ct identical transformation
    Basis b(-3.0, 0.60, 0.70, 11, 24, 1.3);
    arma::mat t = b.getTransformationMatrix(b);
    ASSERT_NEAR(arma::norm(t - arma::eye(t.n_cols, t.n_rows), "inf"), 0.0, 1e-12);
  }
  {
    // 2ct identical transformation
    Basis b( 8.0, 0.60, 0.70, 11, 24, 1.3);
    arma::mat t = b.getTransformationMatrix(b);
    ASSERT_NEAR(arma::norm(t - arma::eye(t.n_cols, t.n_rows), "inf"), 0.0, 1e-12);

  }
}

//==============================================================================

TEST(Basis, getOverlapMatrixHOZ)
{

  Basis b1( 3.0, 0.60, 0.70, 11, 24, 1.3);
  Basis b2(11.0, 0.62, 0.75, 12, 24, 1.0);

  arma::mat overlapz = b1.getOverlapMatrixHOZ(b2);

  // Let's test everything !!!

  Discrete d(&b1, Mesh::gaussLaguerreHermite(200, 200));

  for (INT n_z1 = 0; n_z1 < b1.n_zGlobalMax; n_z1++)
  {
    for (INT d1 = 0; d1 < b1.dMax; d1++)
    {
      INT izd1 = b1.QNn_zMaxd.find({n_z1, d1});

      for (INT n_z2 = 0; n_z2 < b2.n_zGlobalMax; n_z2++)
      {
        for (INT d2 = 0; d2 < b2.dMax; d2++)
        {
          INT izd2 = b2.QNn_zMaxd.find({n_z2, d2});

          double numOver = 0.0;

          arma::vec zFunc = b1.zPart(d.mesh.az.p, n_z1, d1)
                            % b2.zPart(d.mesh.az.p, n_z2, d2)
                            % d.mesh.az.we;

          numOver = arma::accu(zFunc);

          ASSERT_NEAR(numOver, overlapz(izd1, izd2), 1e-11);
        }
      }
    }
  }

  b1 = Basis("examples/42Ca_deformed_2x9.msg.gz");
  b2 = Basis("examples/42Ca_deformed_1x11.msg.gz");
  ASSERT_NEAR(b1.getOverlapMatrixHOZ(b2)(b1.QNn_zMaxd.find({3, 0}), b2.QNn_zMaxd.find({0, 0})), BASIS_GETOVERLAPMATRIXHOZ_TARGET00, 1e-12);

  b1 = Basis("examples/42Ca_deformed_2x9.msg.gz");
  b2 = Basis("examples/42Ca_deformed_2x9.msg.gz");
  ASSERT_NEAR(b1.getOverlapMatrixHOZ(b2)(b1.QNn_zMaxd.find({3, 0}), b2.QNn_zMaxd.find({5, 1})), BASIS_GETOVERLAPMATRIXHOZ_TARGET01, 1e-12);

  b1 = Basis("examples/42Ca_deformed_1x11.msg.gz");
  b2 = Basis("examples/42Ca_deformed_1x11.msg.gz");
  ASSERT_NEAR(arma::norm(b1.getOverlapMatrixHOZ(b2) - arma::eye(b1.QNn_zMaxd.nb, b1.QNn_zMaxd.nb), "inf"), 0.0, 1e-12);
}

//==============================================================================

TEST(Basis, getOverlapMatrixHOR)
{
  Basis b1 = Basis(-1.0, 0.60, 0.70, 11, 24, 1.3);
  Basis b2 = Basis(11.0, 0.62, 0.75, 12, 24, 1.0);

  arma::mat overlapr = b1.getOverlapMatrixHOR(b2);

  // Let's test everything !!!

  Discrete d(&b1, Mesh::gaussLaguerreHermite(200, 200));

  for (INT m1 = 0; m1 < b1.mMax; m1++)
  {
    for (INT n1 = 0; n1 < b1.nMax(m1); n1++)
    {
      INT imn1 = b1.QNmn.find({m1, n1});

      for (INT m2 = 0; m2 < b2.mMax; m2++)
      {
        for (INT n2 = 0; n2 < b2.nMax(m2); n2++)
        {
          INT imn2 = b2.QNmn.find({m2, n2});

          double numOver = 0.0;

          if (m1 == m2)
          {
            arma::vec rFunc = b1.rPart(d.mesh.ax.p, m1, n1)
                              % b2.rPart(d.mesh.ax.p, m2, n2)
                              % d.mesh.ax.p
                              % d.mesh.ax.we;

            numOver = arma::accu(rFunc) * 2.0 * PI;

          }

          ASSERT_NEAR(numOver, overlapr(imn1, imn2), 1e-04); // TODO: why so large ?
        }
      }
    }
  }

  b1 = Basis("examples/42Ca_deformed_1x11.msg.gz");
  b2 = Basis("examples/42Ca_deformed_2x9.msg.gz");
  ASSERT_NEAR(b1.getOverlapMatrixHOR(b2)(b1.QNmn.find({2, 2}), b2.QNmn.find({2, 1})), BASIS_GETOVERLAPMATRIXHOR_TARGET00, 1e-14);
  ASSERT_NEAR(b1.getOverlapMatrixHOR(b2)(b1.QNmn.find({2, 1}), b2.QNmn.find({2, 2})), BASIS_GETOVERLAPMATRIXHOR_TARGET01, 1e-14);
}

//==============================================================================

TEST(Basis, calcTab)
{
  Basis basis("examples/42Ca_deformed_2x9.msg.gz");
  basis.calcTab();
  ASSERT_NEAR(basis.Tab(5, 3),           BASIS_CALCTAB_TARGET00, 1e-12);
  ASSERT_NEAR(basis.Tab(3, 5),           BASIS_CALCTAB_TARGET01, 1e-12);
  ASSERT_NEAR(basis.Tab(2, 7),           BASIS_CALCTAB_TARGET02, 1e-12);
  ASSERT_NEAR(basis.Tab(7, 2),           BASIS_CALCTAB_TARGET03, 1e-12);
  ASSERT_NEAR(basis.tabzd(10, 10, 0, 1), BASIS_CALCTAB_TARGET04, 1e-12);
  ASSERT_NEAR(basis.tabzd( 9, 10, 0, 1), BASIS_CALCTAB_TARGET05, 1e-12);
}

//==============================================================================

TEST(Basis, tabzd)
{
  {
    // 1ct state
    Basis basis("examples/42Ca_deformed_1x11.msg.gz");
    basis.calcTab();
    // Exact quadrature
    Discrete d(&basis, Mesh::gaussLaguerreHermite(20, 30));

    for (UINT ia = 0; ia < basis.HOqn.nb; ia++)
    {
      for (UINT ib = ia; ib < basis.HOqn.nb; ib++)
      {
        INT ma   = basis.HOqn(ia)(0);
        INT na   = basis.HOqn(ia)(1);
        INT n_za = basis.HOqn(ia)(2);
        INT da   = basis.HOqn(ia)(3);
        INT sa   = basis.HOqn(ia)(4);
        INT mb   = basis.HOqn(ib)(0);
        INT nb   = basis.HOqn(ib)(1);
        INT n_zb = basis.HOqn(ib)(2);
        INT db   = basis.HOqn(ib)(3);
        INT sb   = basis.HOqn(ib)(4);

        if (ma != mb) continue; // different m

        if (sa != sb) continue; // different spin

        arma::vec rFunc = basis.rPartNormReduced(d.mesh.ax.p, ma, na)
                          % basis.rPartNormReduced(d.mesh.ax.p, mb, nb)
                          % d.mesh.ax.w;
        arma::vec zFunc = basis.zPartNormReduced(d.mesh.az.p, n_za)
                          % basis.zPartNormReduced(d.mesh.az.p, n_zb)
                          % d.mesh.az.w;
        double expectedResult = (na == nb) ? basis.tabzd(n_za, n_zb, da, db) : 0.0;
        double result = arma::accu(rFunc) * arma::accu(zFunc);
        ASSERT_NEAR(result, expectedResult, 1e-10);
      }
    }
  }

  {
    // 2ct state
    Basis basis("examples/42Ca_deformed_2x9.msg.gz");
    basis.calcTab();
    basis.calcTalmanz();
    // Exact quadrature
    Discrete d(&basis, Mesh::gaussLaguerreHermite(20, 30));

    for (UINT ia = 0; ia < basis.HOqn.nb; ia++)
    {
      for (UINT ib = ia; ib < basis.HOqn.nb; ib++)
      {
        INT ma   = basis.HOqn(ia)(0);
        INT na   = basis.HOqn(ia)(1);
        INT n_za = basis.HOqn(ia)(2);
        INT da   = basis.HOqn(ia)(3);
        INT sa   = basis.HOqn(ia)(4);
        INT mb   = basis.HOqn(ib)(0);
        INT nb   = basis.HOqn(ib)(1);
        INT n_zb = basis.HOqn(ib)(2);
        INT db   = basis.HOqn(ib)(3);
        INT sb   = basis.HOqn(ib)(4);

        if (ma != mb) continue; // different m

        if (sa != sb) continue; // different spin

        if (da == db)
        {
          arma::vec rFunc = basis.rPartNormReduced(d.mesh.ax.p, ma, na)
                            % basis.rPartNormReduced(d.mesh.ax.p, mb, nb)
                            % d.mesh.ax.w;
          arma::vec zFunc = basis.zPartNormReduced(d.mesh.az.p, n_za)
                            % basis.zPartNormReduced(d.mesh.az.p, n_zb)
                            % d.mesh.az.w;
          double expectedResult = (na == nb) ? basis.tabzd(n_za, n_zb, da, db) : 0.0;
          double result = arma::accu(rFunc) * arma::accu(zFunc);
          ASSERT_NEAR(result, expectedResult, 1e-10);
        }
        else
        {
          arma::vec rFunc = basis.rPartNormReduced(d.mesh.ax.p, ma, na)
                            % basis.rPartNormReduced(d.mesh.ax.p, mb, nb)
                            % d.mesh.ax.w;
          double zResult = 0.0;

          for (INT n_zc = 0; n_zc < n_za + n_zb + 1; n_zc++)
          {
            zResult += basis.talmanz(n_za, da, n_zb, db)(n_zc)
                       * arma::accu(  basis.zPartNormReduced(d.mesh.az.p, n_zc)
                                      % basis.zPartNormReduced(d.mesh.az.p, 0)
                                      % d.mesh.az.w  );
          }

          double expectedResult = (na == nb) ? basis.tabzd(n_za, n_zb, da, db) : 0.0;
          double result = arma::accu(rFunc) * zResult;
          ASSERT_NEAR(result, expectedResult, 1e-10);
        }
      }
    }
  }
}

//==============================================================================


TEST(Basis, getDataTree)
{
  // without prefix
  std::string filename = "misc/data/test_2ct_256Fm.hfb3";
  DataTree dataTree(filename);
  Basis basis(dataTree);

  DataTree d = basis.getDataTree();
  DataTree result;
  result.set("basis/nOscil", basis.nOscil);
  result.set("basis/d_0", basis.d_0);
  result.set("basis/b_r", basis.b_r);
  result.set("basis/b_z", basis.b_z);
  result.set("basis/g_q", basis.g_q);
  result.set("basis/n_zMax", basis.n_zMaxImposed);

  ASSERT_TRUE(result == d);

  result.set("basis/d_0", 20.0);
  ASSERT_FALSE(result == d);

  //with prefix
  DataTree prefixD = basis.getDataTree("state/");
  DataTree resultD;
  resultD.set("state/basis/nOscil", basis.nOscil);
  resultD.set("state/basis/d_0", basis.d_0);
  resultD.set("state/basis/b_r", basis.b_r);
  resultD.set("state/basis/b_z", basis.b_z);
  resultD.set("state/basis/g_q", basis.g_q);
  resultD.set("state/basis/n_zMax", basis.n_zMaxImposed);

  ASSERT_TRUE(resultD == prefixD);
}

//==============================================================================

/*
TEST(Basis, info)
{
}
*/

//==============================================================================

TEST(Basis, getBasisDistance)
{
  double precision = 1e-15;

  Basis b1(12.3, 1.3, 1.2, 12, 20, 1.1);
  Basis b2(12.3, 1.3, 1.2, 12, 20, 1.1);

  ASSERT_NEAR(b1.getBasisDistance(b2), 0.0, precision);
}

//==============================================================================

TEST(Basis, calcTalmanz)
{
  {
    // 1ct state
    Basis basis(0.0, 5.3, 5.2, 20, 24, 1.0);
    basis.calcTalmanz();

    double precision = 1e-9;

    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 0), BASIS_CALCTALMANZ_TARGET00, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 1), BASIS_CALCTALMANZ_TARGET01, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 2), BASIS_CALCTALMANZ_TARGET02, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 3), BASIS_CALCTALMANZ_TARGET03, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 4), BASIS_CALCTALMANZ_TARGET04, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 5), BASIS_CALCTALMANZ_TARGET05, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 6), BASIS_CALCTALMANZ_TARGET06, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 7), BASIS_CALCTALMANZ_TARGET07, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 8), BASIS_CALCTALMANZ_TARGET08, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)( 9), BASIS_CALCTALMANZ_TARGET09, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(10), BASIS_CALCTALMANZ_TARGET10, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(11), BASIS_CALCTALMANZ_TARGET11, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(12), BASIS_CALCTALMANZ_TARGET12, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(13), BASIS_CALCTALMANZ_TARGET13, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(14), BASIS_CALCTALMANZ_TARGET14, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(15), BASIS_CALCTALMANZ_TARGET15, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(16), BASIS_CALCTALMANZ_TARGET16, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(17), BASIS_CALCTALMANZ_TARGET17, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(18), BASIS_CALCTALMANZ_TARGET18, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(19), BASIS_CALCTALMANZ_TARGET19, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(20), BASIS_CALCTALMANZ_TARGET20, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(21), BASIS_CALCTALMANZ_TARGET21, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(22), BASIS_CALCTALMANZ_TARGET22, precision);
    ASSERT_NEAR(basis.talmanz(11, 0, 11, 0)(23), BASIS_CALCTALMANZ_TARGET23, precision);
    double z = 4.1;

    for (INT n_za = 0; n_za < basis.n_zGlobalMax + 1; n_za++)
    {
      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax + 1; n_zb++)
      {
        ASSERT_GE(basis.talmanz(n_za, 0, n_zb, 0).n_rows, n_za + n_zb + 1);

        double p = basis.zPartScalar(z, n_za) * basis.zPartScalar(z, n_zb);
        double sum = 0;

        INT n_zcMin = abs(n_za - n_zb);
        INT n_zcMax = n_za + n_zb + 1;

        for (INT n_zc = n_zcMin; n_zc < n_zcMax; n_zc++)
        {
          sum += basis.talmanz(n_za, 0, n_zb, 0)(n_zc)
                 * basis.zPartScalar(z, n_zc);
        }

        sum *= basis.zPartScalar(z, 0);

        ASSERT_NEAR(p, sum, 1e-7);
      }
    }
  }

  {
    // 2ct state
    Basis basis(12.0, 5.3, 5.2, 20, 24, 1.1);
    basis.calcTalmanz();

    double precision = 1e-10;

    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 0), BASIS_CALCTALMANZ_TARGET24, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 1), BASIS_CALCTALMANZ_TARGET25, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 2), BASIS_CALCTALMANZ_TARGET26, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 3), BASIS_CALCTALMANZ_TARGET27, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 4), BASIS_CALCTALMANZ_TARGET28, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 5), BASIS_CALCTALMANZ_TARGET29, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 6), BASIS_CALCTALMANZ_TARGET30, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 7), BASIS_CALCTALMANZ_TARGET31, precision);

    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 3), BASIS_CALCTALMANZ_TARGET32, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 1)( 3), BASIS_CALCTALMANZ_TARGET33, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 0)( 3), BASIS_CALCTALMANZ_TARGET34, precision);
    ASSERT_NEAR(basis.talmanz( 9, 0, 8, 1)( 3), BASIS_CALCTALMANZ_TARGET35, precision);
    ASSERT_NEAR(basis.talmanz( 9, 1, 8, 0)( 3), BASIS_CALCTALMANZ_TARGET36, precision);
    ASSERT_NEAR(basis.talmanz( 9, 1, 8, 1)( 3), BASIS_CALCTALMANZ_TARGET37, precision);
    ASSERT_NEAR(basis.talmanz( 9, 1, 8, 0)( 3), BASIS_CALCTALMANZ_TARGET38, precision);
    ASSERT_NEAR(basis.talmanz( 9, 1, 8, 1)( 3), BASIS_CALCTALMANZ_TARGET39, precision);

    double z = 4.1;

    for (INT n_za = 0; n_za < basis.n_zGlobalMax + 1; n_za++)
    {
      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax + 1; n_zb++)
      {
        for (INT da = 0; da < basis.dMax; da++)
        {
          for (INT db = 0; db < basis.dMax; db++)
          {
            ASSERT_GE(basis.talmanz(n_za, da, n_zb, db).n_rows, n_za + n_zb + 1);

            double dab = (1.0 - da - db) * basis.d_0 / 2.0;
            double p = basis.zPartScalar(z, n_za, da) * basis.zPartScalar(z, n_zb, db);
            double sum = 0;

            INT n_zcMin = 0;
            INT n_zcMax = n_za + n_zb + 1;

            for (INT n_zc = n_zcMin; n_zc < n_zcMax; n_zc++)
            {
              sum += basis.talmanz(n_za, da, n_zb, db)(n_zc)
                     * basis.zPartScalar(z - dab, n_zc);
            }

            sum *= basis.zPartScalar(z - dab, 0);
            ASSERT_NEAR(p, sum, 1e-7);
          }
        }
      }
    }
  }
}

//==============================================================================

TEST(Basis, getFullMoshinskyz)
{
  {
    // 1ct basis
    Basis basis(0.0, 5.3, 5.2, 20, 24, 1.1);
    Multi<arma::vec> fullMoshinskyz = basis.getFullMoshinskyz();

    double z1 = 4.1;
    double z2 = 1.2;
    double Z = (z1 + z2) / sqrt(2.0);
    double z = (z1 - z2) / sqrt(2.0);

    for (INT n_za = 0; n_za < basis.n_zGlobalMax * 2 + 1; n_za++)
    {
      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax * 2 + 1; n_zb++)
      {
        double p = basis.zPartScalar(z1, n_za, 0) * basis.zPartScalar(z2, n_zb, 0);
        double sum = 0;
        bool truncated = false;

        for (INT n_zc = 0; n_zc < n_za + n_zb + 1; n_zc++)
        {
          INT n_zd = n_za + n_zb - n_zc;

          if (n_zc > 40) // (n_zGlobalMax * 4 - 3) is too big, and only moshinskyz()(0) will be used after.
          {
            truncated = true;
            break;
          }

          sum += fullMoshinskyz(n_za, n_zb)(n_zc)
                 * basis.zPartScalar(Z, n_zc)
                 * basis.zPartScalar(z, n_zd);
        }

        if (!truncated)
        {
          ASSERT_NEAR(p, sum, 1e-10);
        }
      }
    }
  }

  {
    // 2ct basis
    Basis basis(12.0, 5.3, 5.2, 20, 24, 1.1);
    Multi<arma::vec>fullMoshinskyz_2ct = basis.getFullMoshinskyz();

    double z1 = 4.1;
    double z2 = 1.2;
    double Z = (z1 + z2) / sqrt(2.0);
    double z = (z1 - z2) / sqrt(2.0);

    for (INT n_za = 0; n_za < basis.n_zGlobalMax * 2 + 1; n_za++)
    {
      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax * 2 + 1; n_zb++)
      {
        for (INT da = 0; da < basis.dMax; da++)
        {
          for (INT db = 0; db < basis.dMax; db++)
          {
            double Dab = (1.0 - da - db) * basis.d_0 / sqrt(2.0);
            double dab = (db - da) * basis.d_0 / sqrt(2.0);
            double p = basis.zPartScalar(z1, n_za, da) * basis.zPartScalar(z2, n_zb, db);
            double sum = 0;
            bool truncated = false;

            for (INT n_zc = 0; n_zc < n_za + n_zb + 1; n_zc++)
            {
              INT n_zd = n_za + n_zb - n_zc;

              if (n_zc > 40) // (n_zGlobalMax * 4 - 3) is too big, and only moshinskyz()(0) will be used after.
              {
                truncated = true;
                break;
              }

              sum += fullMoshinskyz_2ct(n_za, n_zb)(n_zc)
                     * basis.zPartScalar(Z - Dab, n_zc)
                     * basis.zPartScalar(z - dab, n_zd);
            }

            if (!truncated)
            {
              ASSERT_NEAR(p, sum, 1e-10);
            }
          }
        }
      }
    }
  }
}

//==============================================================================

TEST(Basis, getFullMoshinskyr)
{
  // 2ct basis
  Basis basis(0.0, 5.3, 5.2, 20, 24, 1.1);

  double r1 = 1.2;
  double r2 = 1.1;
  double R = (r1 + r2) / sqrt(2.0);
  double r = (r1 - r2) / sqrt(2.0);
  INT nbok0 = 0;
  INT nbok1 = 0;
  INT nbtot = 0;

  for (INT ma = 0; ma < basis.mMax; ma++)
  {
    for (INT na = 0; na < basis.nGlobalMax; na++)
    {
      INT xa = na * 2 + abs(ma);

      for (INT mc = 0; mc < basis.mMax; mc++)
      {
        for (INT nc = 0; nc < basis.nGlobalMax; nc++)
        {
          if (Tools::irand(1000) != 0) continue;

          INT xc = nc * 2 + abs(mc);
          double p = basis.rPartScalar(r1, ma, na) * basis.rPartScalar(r2, mc, nc);
          double sum0 = 0;
          double sum1 = 0;

          for (INT nmu = 0; nmu < na + nc + 1; nmu++)
          {
            for (INT mmu = -na - nc + nmu; mmu < na + nc - nmu + mc + ma + 1; mmu++)
            {
              INT xmu = nmu * 2 + abs(mmu);
              INT ml = ma + mc - mmu;
              INT nl = (xa + xc - xmu - abs(ml)) / 2;
              double mosh0 = basis.getFullMoshinskyr(ma, na, mc, nc, mmu, nmu, ml, nl);
              double term0 = mosh0 * basis.rPartScalar(R, mmu, nmu) * basis.rPartScalar(r, ml, nl);
              sum0 += term0;
            }
          }

          for (INT nmu = 0; nmu < na + nc + MIN(mc, ma) + 1; nmu++)
          {
            for (INT mmu = -na - nc + nmu - ma; mmu < na + nc - nmu + mc + 1; mmu++)
            {
              INT xmu = nmu * 2 + abs(mmu);
              INT ml = -ma + mc - mmu;
              INT nl = (xa + xc - xmu - abs(ml)) / 2;
              double mosh0 = basis.getFullMoshinskyr(-ma, na, mc, nc, mmu, nmu, ml, nl);
              double term0 = mosh0 * basis.rPartScalar(R, mmu, nmu) * basis.rPartScalar(r, ml, nl);
              sum1 += term0;
            }
          }

          if (fabs(p - sum0) < 1e-11) nbok0++;

          if (fabs(p - sum1) < 1e-11) nbok1++;

          nbtot++;
        }
      }
    }
  }

  ASSERT_EQ(nbok0, nbtot);
  ASSERT_EQ(nbok1, nbtot);
}

//==============================================================================

/*
TEST(Basis, niceStr)
{
}
*/

//==============================================================================


TEST(Basis, operator_equality)
{
  Basis a(10.0, 1.8, 1.8, 11, 24, 1.3);
  Basis b(10.0, 1.8, 1.8, 11, 24, 1.3);

  ASSERT_TRUE(a == b);
  ASSERT_EQ(a, b);

  Basis c(10.0, 1.8, 1.8, 22, 24, 1.3);
  ASSERT_FALSE(a == c);

  Basis d(10.0 + 1e-8, 1.8, 1.8, 11, 24, 1.3);
  Basis e(10.0 + 1e-9, 1.8, 1.8, 11, 24, 1.3);
  ASSERT_FALSE(a == d);
  ASSERT_TRUE(a == e);

  b.b_r = 1.3;
  ASSERT_FALSE(a == b);

  Basis copy = a;
  ASSERT_TRUE(a == copy);
}

//==============================================================================


TEST(Basis, operator_inequality)
{
  Basis a(10.0, 1.8, 1.8, 11, 24, 1.3);
  Basis b(10.0, 1.8, 1.8, 11, 24, 1.3);
  ASSERT_FALSE(a != b);

  Basis c(10.0, 1.8, 1.8, 22, 24, 1.3);
  ASSERT_TRUE(a != c);

  Basis d(10.0 + 1e-8, 1.8, 1.8, 11, 24, 1.3);
  Basis e(10.0 + 1e-9, 1.8, 1.8, 11, 24, 1.3);
  ASSERT_TRUE(a != d);
  ASSERT_FALSE(a != e);

  b.b_r = 1.3;
  ASSERT_TRUE(a != b);

  Basis copy = a;
  ASSERT_FALSE(a != copy);
}

//==============================================================================

TEST(Basis, calcTalmanr)
{
  {
    // 1ct state
    Basis basis(0.0, 5.3, 5.2, 12, 28, 1.1);

    basis.calcTalmanr();

    double r = 4.1;
    double e = basis.rPartScalar(r, 0, 0);

    for (INT m_a = 0; m_a < basis.mMax + 1; m_a++)
    {
      INT Ma = MAX(m_a - 1,0);
      for (INT n_a = 0; n_a < basis.nMax(Ma) + 1; n_a++)
      {
       	for (INT m_b = 0; m_b < basis.mMax + 1; m_b++)
        {
          INT Mb = MAX(m_b - 1,0);
          INT m_c = abs(m_b - m_a);

          for (INT n_b = 0; n_b < basis.nMax(Mb) + 1; n_b++)
          {
            INT n_c_min = (abs(2 * n_a + m_a - 2 * n_b - m_b) - abs(m_a - m_b)) / 2;
            INT n_c_max = (abs(2 * n_a + m_a + 2 * n_b + m_b) - abs(m_a - m_b)) / 2 + 1;

            if (n_c_min < 0) n_c_min = 0;

            double p = basis.rPartScalar(r, m_a, n_a) * basis.rPartScalar(r, m_b, n_b);
            double sum = 0;

            for (INT n_c = n_c_min; n_c < n_c_max; n_c++)
            {
              sum += basis.talmanr(m_a, n_a, m_b, n_b)(n_c) * basis.rPartScalar(r, m_c, n_c);
            }

            sum *= e;
            ASSERT_NEAR(p, sum, 1e-12);
          }
        }
      }
    }
  }

  {
    // 1ct state - Negative part
    Basis basis(0.0, 5.3, 5.2, 12, 28, 1.1);

    basis.calcTalmanr();

    double r = 4.1;
    double e = basis.rPartScalar(r, 0, 0);

    for (INT m_a = 0; m_a < basis.mMax; m_a++)
    {
      for (INT n_a = 0; n_a < basis.nMax(m_a); n_a++)
      {
        for (INT m_b = 0; m_b < basis.mMax; m_b++)
        {
          INT m_c = abs(m_b + m_a);

          for (INT n_b = 0; n_b < basis.nMax(m_b); n_b++)
          {
            INT n_c_min = (abs(2 * n_a + m_a - 2 * n_b - m_b) - abs(-m_a - m_b)) / 2;
            INT n_c_max = (abs(2 * n_a + m_a + 2 * n_b + m_b) - abs(-m_a - m_b)) / 2 + 1;

            if (n_c_min < 0) n_c_min = 0;

            double p = basis.rPartScalar(r, m_a, n_a) * basis.rPartScalar(r, m_b, n_b);
            double sum = 0;

            for (INT n_c = n_c_min; n_c < n_c_max; n_c++)
            {
              sum += basis.talmanr(m_a, n_a, -m_b, n_b)(n_c) * basis.rPartScalar(r, m_c, n_c);
            }

            sum *= e;

            ASSERT_NEAR(p, sum, 1e-13);
          }
        }
      }
    }
  }
}

//==============================================================================

TEST(Basis, calcMoshinskyr)
{

  Basis basis(0.0, 2.3, 2.2, 14, 24, 1.1);

  basis.calcMoshinskyr();

  INT M = basis.mMax;
  INT maxn = basis.Nmaxr;

  for (INT mmu = -2 * M + 2; mmu < 2 * M - 1; mmu++)
  {
    for (INT nmu = 0; nmu < maxn; nmu++)
    {
      for (INT nnu = 0; nnu < maxn; nnu++)
      {
        ASSERT_NEAR(       basis.moshinskyr(mmu, nmu)(      nnu),
                    basis.getFullMoshinskyr(mmu, nmu, -mmu, nnu, 0, 0, 0, nmu + nnu + abs(mmu)), 1e-14);
      }
    }
  }
}

//==============================================================================


TEST(Basis, calcMoshinskyz)
{
  {
    Basis basis(0.0, 5.3, 5.2, 20, 24, 1.1);
    basis.calcMoshinskyz();

    Multi<arma::vec> fullMoshinskyz = basis.getFullMoshinskyz();

    for (INT n_za = 0; n_za < basis.n_zGlobalMax * 2 + 1; n_za++)
    {
      for (INT n_zb = 0; n_zb < basis.n_zGlobalMax * 2 + 1; n_zb++)
      {
        ASSERT_NEAR(basis.moshinskyz(n_za, n_zb), fullMoshinskyz(n_za,n_zb)(0), 1e-15);
      }
    }
  }
}

