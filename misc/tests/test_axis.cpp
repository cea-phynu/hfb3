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
 *  \brief Test suite for the Axis class.
 */

#include "gtest/gtest.h"
#include "axis.h"
#include "tools.h"

//==============================================================================

TEST(Axis, Axis)
{
  Axis a(Axis::REGULAR, 20, -1.0, 1.0);
  ASSERT_TRUE(a.nb == 20);
  ASSERT_TRUE(a.p.n_rows == 20);
  a = Axis(Axis::GAUSS_LEGENDRE, 10, -1.0, 1.0);

  ASSERT_NEAR(a.p(3), -4.3339539412924719E-001, 1e-15);
  ASSERT_NEAR(a.w(3),  2.6926671930999635E-001, 1e-15);
  a = Axis(Axis::GAUSS_LAGUERRE, 20);
  ASSERT_NEAR(a.p(3),  1.7073065310283438E+000, 1e-15);
  ASSERT_NEAR(a.w(6),  6.2025508445722370E-003, 1e-15);
  a = Axis(Axis::GAUSS_HERMITE, 30);
  ASSERT_NEAR(a.p(26),  4.9889189685899440E+000, 1e-15);
  ASSERT_NEAR(a.w(22),  2.9387252289229878E-005, 1e-15);

  arma::vec func;

  // Calculate an integral using a Gauss-Legendre quadrature.
  a = Axis(Axis::GAUSS_LEGENDRE, 60, -1.0, 1.0);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.w;
  ASSERT_NEAR(arma::accu(func),  1.2376439265798188E000, 1e-10);

  // Calculate an integral using a Gauss-Laguerre quadrature.
  a = Axis(Axis::GAUSS_LAGUERRE, 60);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.w % arma::exp(a.p);
  ASSERT_NEAR(arma::accu(func),  6.7164671087916214E-001, 1e-15);

  // Calculate an integral using a Gauss-Hermite quadrature.
  a = Axis(Axis::GAUSS_HERMITE, 60);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.w % arma::exp(a.p % a.p);
  ASSERT_NEAR(arma::accu(func),  1.3432934170262217E+000, 1e-14);

  // Calculate an integral using a shifted Gauss-Legendre quadrature.
  a = Axis(Axis::GAUSS_LEGENDRE, 60, -3.0, -1.0);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.w;
  ASSERT_NEAR(arma::accu(func),  5.2822948618813680E-002, 1e-10);

  // Calculate an integral using a Gauss-Laguerre quadrature.
  a = Axis(Axis::GAUSS_LAGUERRE, 60);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.we;
  ASSERT_NEAR(arma::accu(func),  6.7164671087916214E-001, 1e-15);

  // Calculate an integral using a Gauss-Hermite quadrature.
  a = Axis(Axis::GAUSS_HERMITE, 60);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.we;
  ASSERT_NEAR(arma::accu(func),  1.3432934170262217E+000, 1e-14);

  // Calculate an integral using a shifted Gauss-Legendre quadrature.
  a = Axis(Axis::GAUSS_LEGENDRE, 60, -3.0, -1.0);
  func = arma::exp(a.p % a.p * -1) / (1.0 + a.p % a.p) % a.we;
  ASSERT_NEAR(arma::accu(func),  5.2822948618813680E-002, 1e-10);
}

//==============================================================================

TEST(Axis, operator_equality)
{
  Axis a(Axis::REGULAR, 20, -1.0, 1.0);
  Axis b(Axis::REGULAR, 20, -1.0, 1.0);
  ASSERT_TRUE(a == b);
  ASSERT_EQ(a, b);

  Axis c(Axis::REGULAR, 20, -2.0, 1.0);
  ASSERT_FALSE(a == c);

  Axis d(Axis::GAUSS_LEGENDRE, 30);
  Axis e(Axis::GAUSS_LEGENDRE, 30);
  ASSERT_TRUE(d == e);

  e.w[3] = 1.7;
  ASSERT_FALSE(d == e);

  e = Axis(Axis::GAUSS_LEGENDRE, 30);
  e.we[3] = 1.7;
  ASSERT_FALSE(d == e);

  Axis f(Axis::GAUSS_LEGENDRE, 40);
  ASSERT_FALSE(d == f);
}

//==============================================================================

TEST(Axis, operator_inequality)
{
  Axis a(Axis::REGULAR, 20, -1.0, 1.0);
  Axis b(Axis::REGULAR, 20, -1.0, 1.0);
  ASSERT_FALSE(a != b);

  Axis c(Axis::REGULAR, 20, -2.0, 1.0);
  ASSERT_TRUE(a != c);

  Axis d(Axis::GAUSS_LEGENDRE, 30);
  Axis e(Axis::GAUSS_LEGENDRE, 30);
  ASSERT_FALSE(d != e);

  e.w[3] = 1.7;
  ASSERT_TRUE(d != e);

  e = Axis(Axis::GAUSS_LEGENDRE, 30);
  e.we[3] = 1.7;
  ASSERT_TRUE(d != e);


  Axis f(Axis::GAUSS_LEGENDRE, 40);
  ASSERT_TRUE(d != f);
}
