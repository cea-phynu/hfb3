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
  DataTree d0;

  arma::vec V;
  V.randu(20);

  arma::mat M;
  M.randu(20, 13);

  Multi<IVEC > multiV;
  multiV(1) = arma::randi<IVEC>(3);

  Multi<IMAT > imultiM;
  imultiM(2,7,8) = arma::randi<IMAT>(3, 3);
  imultiM(2,4) = IMAT();

  Multi<arma::mat> multiM;
  multiM(2) = arma::randu(3, 3);
  multiM(4,9) = arma::randu(3, 4);

  d0.set("string", "text message");
  d0.set("string", "text message override");
  d0.set("int", 4);
  d0.set("double", 5.8);
  d0.set("vec", V);
  d0.set("ivec", IVEC());
  d0.set("mat", M);
  d0.set("multiVec", multiV);
  d0.set("multiIMat", imultiM);
  d0.set("multiMat", multiM);

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
  DataTree d0;

  // Creation of the DataTree
  arma::vec V;
  V.randu(20);

  IVEC iv1 = arma::randi<IVEC>(5);

  arma::mat M;
  M.randu(20, 13);

  Multi<IVEC > multiIV;
  multiIV(1) = IVEC();

  Multi<IMAT > multiIM;
  multiIM(2, 4) = arma::randi<IMAT >(3, 3);

  Multi<arma::mat> multiM;
  multiM(2, 4, 5) = arma::randu(3, 3);

  d0.set("string", "text message");
  d0.set("int", 4);
  d0.set("double", 5.8);
  d0.set("vec", V);
  d0.set("ivec1", iv1);
  d0.set("ivec2", IVEC());
  d0.set("mat", M);
  d0.set("multIVec", multiIV);
  d0.set("multiIMat", multiIM);
  d0.set("multiMat", multiM);

  // Data saving
  std::string content = IOmsgp().serializeDataTree(d0);

  // Data recovery
  DataTree result = IOmsgp().fromContent(content);

  ASSERT_EQ(d0, result);
}

//============================================================================

TEST(IOmsgp, serializeDataTree)
{
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
