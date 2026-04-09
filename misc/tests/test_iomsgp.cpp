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
 *  \brief Test suite for the IOmsgp class.
 */

#include "gtest/gtest.h"
#include "io_msgp.h"
#include "datatree.h"

//============================================================================

/*
TEST(IOmsgp, IOmsgp)
{
}
*/

//============================================================================

TEST(IOmsgp, saveDataTree)
{
  forceValidDataTree = true;

  DataTree d0;

  // Creation of the DataTree

  d0.set("bool", false);
  d0.set("int", 4);
  d0.set("double", 5.8);
  d0.set("string", "text message");

  VEC V;
  V.randu(20);
  d0.set("vec", V);
  IVEC iv1 = arma::randi<IVEC>(5);
  d0.set("ivec1", iv1);
  d0.set("ivec2", IVEC());
  UVEC uv1 = arma::randi<UVEC>(5);
  d0.set("uvec1", uv1);
  d0.set("uvec2", UVEC());

  MAT M;
  M.randu(20, 13);
  d0.set("mat", M);
  IMAT im1 = arma::randi<IMAT>(5, 5);
  d0.set("imat1", iv1);
  d0.set("imat2", IMAT());
  UMAT um1 = arma::randi<UMAT>(5, 5);
  d0.set("umat1", uv1);
  d0.set("umat2", UMAT());

  CUBE C;
  C.randu(2, 13, 5);
  d0.set("cube", C);
  ICUBE ic1 = arma::randi<ICUBE>(5, 5, 5);
  d0.set("icube1", ic1);
  d0.set("icube2", ICUBE());
  UCUBE uc1 = arma::randi<UCUBE>(5, 5, 5);
  d0.set("ucube1", uc1);
  d0.set("ucube2", UCUBE());

  Multi<VEC  > multiV;
  multiV(1) = VEC();
  d0.set("multiVec", multiV);
  Multi<IVEC > multiIV;
  multiIV(1) = IVEC();
  d0.set("multiIV", multiIV);
  Multi<UVEC > multiUV;
  multiUV(1) = UVEC();
  d0.set("multiUV", multiUV);

  Multi<MAT  > multiM;
  multiM(1) = MAT();
  d0.set("multiM", multiM);
  Multi<IMAT > multiIM;
  multiIM(1) = IMAT();
  d0.set("multiIM", multiIM);
  Multi<UMAT > multiUM;
  multiUM(1) = UMAT();
  d0.set("multiUM", multiUM);

  Multi<CUBE  > multiC;
  multiC(1) = CUBE();
  d0.set("multiC", multiC);
  Multi<ICUBE  > multiIC;
  multiIC(1) = ICUBE();
  d0.set("multiIC", multiIC);
  Multi<UCUBE  > multiUC;
  multiUC(1) = UCUBE();
  d0.set("multiUC", multiUC);

  Multi<double > multiDouble;
  multiDouble(1, 4) = 1.5;
  d0.set("multDouble", multiDouble);

  std::string fileName = "/tmp/test_hfb3.msg.gz";
  IOmsgp().saveDataTree(d0, fileName);

  DataTree d1 = DataTree(fileName);

  ASSERT_EQ(d0, d1);

  d0.set("mat2", M);

  ASSERT_NE(d0, d1);
}

//============================================================================

TEST(IOmsgp, fromContent)
{
  forceValidDataTree = true;

  DataTree d0;

  // Creation of the DataTree

  d0.set("bool", false);
  d0.set("int", 4);
  d0.set("double", 5.8);
  d0.set("string", "text message");

  VEC V;
  V.randu(20);
  d0.set("vec", V);
  IVEC iv1 = arma::randi<IVEC>(5);
  d0.set("ivec1", iv1);
  d0.set("ivec2", IVEC());
  UVEC uv1 = arma::randi<UVEC>(5);
  d0.set("uvec1", uv1);
  d0.set("uvec2", UVEC());

  MAT M;
  M.randu(20, 13);
  d0.set("mat", M);
  IMAT im1 = arma::randi<IMAT>(5, 5);
  d0.set("imat1", iv1);
  d0.set("imat2", IMAT());
  UMAT um1 = arma::randi<UMAT>(5, 5);
  d0.set("umat1", uv1);
  d0.set("umat2", UMAT());

  CUBE C;
  C.randu(2, 13, 5);
  d0.set("cube", C);
  ICUBE ic1 = arma::randi<ICUBE>(5, 5, 5);
  d0.set("icube1", ic1);
  d0.set("icube2", ICUBE());
  UCUBE uc1 = arma::randi<UCUBE>(5, 5, 5);
  d0.set("ucube1", uc1);
  d0.set("ucube2", UCUBE());

  Multi<VEC  > multiV;
  multiV(1) = VEC();
  d0.set("multiVec", multiV);
  Multi<IVEC > multiIV;
  multiIV(1) = IVEC();
  d0.set("multiIV", multiIV);
  Multi<UVEC > multiUV;
  multiUV(1) = UVEC();
  d0.set("multiUV", multiUV);

  Multi<MAT  > multiM;
  multiM(1) = MAT();
  d0.set("multiM", multiM);
  Multi<IMAT > multiIM;
  multiIM(1) = IMAT();
  d0.set("multiIM", multiIM);
  Multi<UMAT > multiUM;
  multiUM(1) = UMAT();
  d0.set("multiUM", multiUM);

  Multi<CUBE  > multiC;
  multiC(1) = CUBE();
  d0.set("multiC", multiC);
  Multi<ICUBE  > multiIC;
  multiIC(1) = ICUBE();
  d0.set("multiIC", multiIC);
  Multi<UCUBE  > multiUC;
  multiUC(1) = UCUBE();
  d0.set("multiUC", multiUC);

  Multi<double > multiDouble;
  multiDouble(1, 4) = 1.5;
  d0.set("multDouble", multiDouble);

  // Data saving
  std::string content = IOmsgp().serializeDataTree(d0);

  // Data recovery
  DataTree result = IOmsgp().fromContent(content);

  ASSERT_EQ(d0, result);
}

//============================================================================

TEST(IOmsgp, serializeDataTree)
{
  forceValidDataTree = true;

  DataTree d0;

  d0.set("int", 1);
  d0.set("dbl"  , 2.0);
  d0.set("str"  , "3.00");

  arma::vec  tempVec  = arma::ones(3) + 0.1;
  arma::mat  tempMat  = arma::ones(3, 3);
  arma::cube tempCube = arma::ones(3, 3, 3);
  IVEC  tempIVec  = arma::ones<IVEC >(3);
  IMAT  tempIMat  = arma::ones<IMAT >(3, 3);
  ICUBE tempICube = arma::ones<ICUBE >(3, 3, 3);
  Multi<IVEC > tempMIVec;
  tempMIVec(1) = arma::ones<IVEC >(3);
  Multi<IMAT > tempMIMat;
  tempMIMat(2,8,0) = arma::ones<IMAT >(3, 3);
  Multi<ICUBE > tempMICube;
  tempMICube(3,5,1) = arma::ones<ICUBE >(3, 3, 3);
  Multi<arma::vec> tempMVec;
  tempMVec(1,1,2) = arma::ones(3);
  Multi<arma::mat> tempMMat;
  tempMMat(2) = arma::ones(3, 3);
  Multi<arma::cube> tempMCube;
  tempMCube(3,6,7) = arma::ones(3, 3, 3);

  d0.set("vec"  , tempVec);
  d0.set("mat"  , tempMat);
  d0.set("cube" , tempCube);
  d0.set("ivec" , tempIVec);
  d0.set("imat" , tempIMat);
  d0.set("icube", tempICube);
  d0.set("multiIVec", tempMIVec);
  d0.set("multiIMat", tempMIMat);
  d0.set("multiICube", tempMICube);
  d0.set("multiVec", tempMVec);
  d0.set("multiMat", tempMMat);
  d0.set("multiCube", tempMCube);

  std::string bytes = IOmsgp::serializeDataTree(d0);

  DataTree d1 = DataTree::fromContent(bytes);

  ASSERT_EQ(d0, d1);
}
