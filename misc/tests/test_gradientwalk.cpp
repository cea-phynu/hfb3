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
 *  \brief Test suite for the GaussianProcess class.
 */

#include "gtest/gtest.h"
#include "gradientwalk.h"
#include "tools.h"

//==============================================================================

double metaFunc(const arma::vec &coords, const arma::vec &shift)
{
  double result = 0.0;

  for (int i = 0; i < coords.n_rows; i++)
  {
    result -= pow((coords(i) - shift(i)), 2);
  }

  return result;
}

//==============================================================================
//
arma::vec loop(GradientWalk &gw, arma::vec minimumPos)
{

  arma::vec next;

  INT iter = 0;
  double convergence = 999;
  while(convergence > 1e-50)
  {
    arma::vec x = gw.getEval();
    if (!x.empty()) {
      double val = metaFunc(x, minimumPos);
      gw.addEval(x, val);
      Tools::info(PF("[%04d] ", iter) + Tools::vecToStr(x) + PF(" %10.6f", val));
    }

    convergence = gw.getConvergence();

    iter++;
  }

  return gw.getOptimum();
}

//==============================================================================

TEST(GradientWalk, optim2d)
{
  // msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO};
  // useColors = true;

  arma::vec minimumPos = {2.06,  2.50};
  arma::vec startPos   = {1.70,  1.60};

  GradientWalk gw;

  gw.xMin           = { 1.5, 1.5};
  gw.xMax           = { 3.0, 3.0};
  gw.xStep          = { 0.01, 0.01};
  gw.initialCoords = startPos;
  gw.init();

  arma::vec finalPos = loop(gw, minimumPos);

  Tools::info("Final pos: " + Tools::vecToStr(finalPos));

  ASSERT_NEAR(arma::norm(finalPos - minimumPos, "inf"), 0.0, 1e-01);

  double val = metaFunc(finalPos, minimumPos);
  ASSERT_NEAR(val, 0.0, 2e-01);
}

//==============================================================================

TEST(GradientWalk, optim3d)
{
  // msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO};
  // useColors = true;

  arma::vec minimumPos = {0.1,  0.2, 0.3};
  arma::vec startPos   = {0.5, -0.4, 4.3};

  GradientWalk gw;

  gw.xMin           = { -1.0, -1.0,  1.0};
  gw.xMax           = {  1.0,  1.0,  5.0};
  gw.xStep          = {  0.1,  0.1,  0.1};
  gw.initialCoords  = startPos;
  gw.init();

  arma::vec finalPos = loop(gw, minimumPos);

  Tools::mesg("UNITST", "Final pos: " + Tools::vecToStr(finalPos));

  ASSERT_NEAR(arma::norm(finalPos - minimumPos, "inf"), 0.0, 1e-04);

  double val = metaFunc(finalPos, minimumPos);
  ASSERT_NEAR(val, 0.0, 2e-01);
}


//==============================================================================

// TEST(GaussianProcess, optim3d_realistic)
// {
//   // msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO};
//   // useColors = true;
//
//   arma::vec minimumPos = {2.06,  2.50, 8.60};
//   arma::vec startPos   = {1.9, 2.1, 8.7};
//
//   GaussianProcess gp;
//
//   gp.xMin           = startPos * 0.9;
//   gp.xMax           = startPos * 1.1;
//   gp.xStep          = { 0.01,  0.01,  0.01};
//   gw.initialCoords  = startPos;
//   gw.init();
//
//   arma::mat init = gp.getStartingPoints(startPos);
//   INT iter = 0;
//
//   for (uint i = 0; i < init.n_cols; i++)
//   {
//     arma::vec x = init.col(i);
//     double val = metaFunc(x, minimumPos);
//
//     gp.addEval(x, val);
//
//     // if (i % 2 == 0) gp.addEval(x, val);
//     // else            gp.addBadEval(x);
//
//     Tools::info(PF("[%04d] ", iter) + Tools::vecToStr(x) + PF(" %10.6f", val));
//     iter++;
//   }
//
//   double det = gp.calcInference(GaussianProcess::CHOL);
//   // gp.checkSpace();
//
//   while(det > 1e-50)
//   {
//
//     arma::mat toCalc;
//
//     toCalc = gp.getPointsToEval(1, "min");
//
//     arma::vec x = toCalc.col(0);
//     double val = metaFunc(x, minimumPos);
//
//     if (iter % 5 == 0) gp.addBadEval(x);
//     else               gp.addEval(x, val);
//
//     det = gp.calcInference(GaussianProcess::CHOL);
//     gp.checkSpace();
//     Tools::mesg("UNITST", PF("EVAL [%04d] ", iter) + Tools::vecToStr(x) + PF(" %10.6f %e", val, det));
//
//     iter++;
//   }
//   arma::vec finalPos = gp.getOptimum();
//
//   Tools::mesg("UNITST", "Final pos: " + Tools::vecToStr(finalPos));
//
//   ASSERT_NEAR(arma::norm(finalPos - minimumPos, "inf"), 0.0, 1e-01);
//
//   double val = metaFunc(finalPos, minimumPos);
//   ASSERT_NEAR(val, 0.0, 2e-01);
// }

