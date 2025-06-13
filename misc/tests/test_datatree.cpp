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
 *  \brief Test suite for the Datatree class.
 */

#include "gtest/gtest.h"
#include "datatree.h"

//==============================================================================

TEST(DataTree, DataTree)
{
  DataTree d;
  std::string result;
  // Example std::strings
  d.set("songs/first", std::string("La fille du Bédouin..."));
  d.set("songs/second", std::string("C'est un fameux trois-mâts..."));
  // Example cylindrical basis
  d.set("basis0/b_r", 0.118656);
  d.set("basis0/b_z", 0.253325);
  d.set("basis0/d_0", 11.5285);
  d.set("basis0/nOscil", 11);
  d.set("basis0/n_zMaxImposed", 24);
  d.set("basis0/gq", 1.6);
  d.set("basis0/type", std::string("Cylindrical"));

  // Example triaxial basis
  d.set("basis1/beta", 0.8);
  d.set("basis1/gamma", 112.0 / 180 * PI);
  d.set("basis1/d_0", 10);
  d.set("basis1/nOscil", 16);
  d.set("basis1/type", std::string("Triaxial"));

  // Example arma::vec
  arma::vec testVec = arma::randu(12);
  d.set("testVec", testVec);
  // Example arma::mat
  arma::mat testMat = arma::randu(12, 8);
  d.set("testMat", testMat);
  // Example arma::cube
  arma::cube testCube = arma::randu(12, 2, 7);
  d.set("testCube", testCube);
  // Example IVEC
  IVEC testIVec = arma::randi<IVEC >(12);
  d.set("testIVec", testIVec);
  // Example IMAT
  IMAT testIMat = arma::randi<IMAT >(12, 8);
  d.set("testIMat", testIMat);
  // Example ICUBE>
  ICUBE testICube = arma::randi<ICUBE >(12, 2, 7);
  d.set("testICube", testICube);
  // Example Multi<IVEC >
  Multi<IVEC > testIMultiVec;
  testIMultiVec(1) = arma::randi<IVEC >(12);
  testIMultiVec(2) = arma::randi<IVEC >(12);
  d.set("testIMultiVec", testIMultiVec);
  // Example Multi<IMAT >
  Multi<IMAT > testIMultiMat;
  testIMultiMat(96) = arma::randi<IMAT >(12, 8);
  testIMultiMat(5) = arma::randi<IMAT >(12, 8);
  d.set("testIMultiMat", testIMultiMat);
  // Example Multi<ICUBE>
  Multi<ICUBE > testIMultiCube;
  testIMultiCube(7) = arma::randi<ICUBE >(12, 2, 7);
  testIMultiCube(16) = arma::randi<ICUBE >(12, 2, 7);
  d.set("testIMultiCube", testIMultiCube);
  // Example Multi<arma::vec>
  Multi<arma::vec> testMultiVec;
  testMultiVec(1) = arma::randu(12);
  testMultiVec(2) = arma::randu(12);
  d.set("testMultiVec", testMultiVec);
  // Example Multi<arma::mat>
  Multi<arma::mat> testMultiMat;
  testMultiMat(96) = arma::randu(12, 8);
  testMultiMat(5) = arma::randu(12, 8);
  d.set("testMultiMat", testMultiMat);
  // Example Multi<arma::cube>
  Multi<arma::cube> testMultiCube;
  testMultiCube(7) = arma::randu(12, 2, 7);
  testMultiCube(16) = arma::randu(12, 2, 7);
  d.set("testMultiCube", testMultiCube);

  d.get(result, "songs/second");
  ASSERT_EQ(result, std::string("C'est un fameux trois-mâts..."));
  ASSERT_EQ(d.getS("songs/second"), std::string("C'est un fameux trois-mâts..."));

  d.get(result, "songs/first");
  ASSERT_EQ(result, std::string("La fille du Bédouin..."));
  ASSERT_EQ(d.getS("songs/first") , std::string("La fille du Bédouin..."));

  INT resultInt = -1;
  d.get(resultInt, "basis0/nOscil");
  ASSERT_EQ(resultInt, 11);
  ASSERT_EQ(d.getI("basis0/nOscil"), 11);

  double resultDouble = 0.0;
  d.get(resultDouble, "basis0/gq");
  ASSERT_NEAR(resultDouble, 1.6, 1e-16);
  ASSERT_NEAR(d.getD("basis0/gq"), 1.6, 1e-16);

  arma::vec testVec2;
  d.get(testVec2, "testVec");
  ASSERT_NEAR(arma::abs(testVec2 - testVec).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getV("testVec") - testVec).max(), 0.0, 1e-13);

  arma::mat testMat2;
  d.get(testMat2, "testMat");
  ASSERT_NEAR(arma::abs(testMat2 - testMat).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getM("testMat") - testMat).max(), 0.0, 1e-13);

  arma::cube testCube2;
  d.get(testCube2, "testCube");
  ASSERT_NEAR(arma::abs(testCube2 - testCube).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getC("testCube") - testCube).max(), 0.0, 1e-13);

  IVEC testIVec2;
  d.get(testIVec2, "testIVec");
  ASSERT_NEAR(arma::abs(testIVec2 - testIVec).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getIV("testIVec") - testIVec).max(), 0.0, 1e-13);

  IMAT testIMat2;
  d.get(testIMat2, "testIMat");
  ASSERT_NEAR(arma::abs(testIMat2 - testIMat).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getIM("testIMat") - testIMat).max(), 0.0, 1e-13);

  ICUBE testICube2;
  d.get(testICube2, "testICube");
  ASSERT_NEAR(arma::abs(testICube2 - testICube).max(), 0.0, 1e-13);
  ASSERT_NEAR(arma::abs(d.getIC("testICube") - testICube).max(), 0.0, 1e-13);

  Multi<IVEC > testIMultiVec2;
  d.get(testIMultiVec2, "testIMultiVec");
  for(auto &key : testIMultiVec.getKeys())
  {
    ASSERT_NEAR(arma::abs(testIMultiVec2(key) - testIMultiVec(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMIV("testIMultiVec")(key) - testIMultiVec(key)).max(), 0.0, 1e-13);
  }

  Multi<IMAT > testIMultiMat2;
  d.get(testIMultiMat2, "testIMultiMat");
  for(auto &key : testIMultiMat.getKeys())
  {
    ASSERT_NEAR(arma::abs(testIMultiMat2(key) - testIMultiMat(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMIM("testIMultiMat")(key) - testIMultiMat(key)).max(), 0.0, 1e-13);
  }

  Multi<ICUBE > testIMultiCube2;
  d.get(testIMultiCube2, "testIMultiCube");
  for(auto &key : testIMultiCube.getKeys())
  {
    ASSERT_NEAR(arma::abs(testIMultiCube2(key) - testIMultiCube(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMIC("testIMultiCube")(key) - testIMultiCube(key)).max(), 0.0, 1e-13);
  }

  Multi<arma::vec> testMultiVec2;
  d.get(testMultiVec2, "testMultiVec");
  for(auto &key : testMultiVec.getKeys())
  {
    ASSERT_NEAR(arma::abs(testMultiVec2(key) - testMultiVec(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMV("testMultiVec")(key) - testMultiVec(key)).max(), 0.0, 1e-13);
  }

  Multi<arma::mat> testMultiMat2;
  d.get(testMultiMat2, "testMultiMat");
  for(auto &key : testMultiMat.getKeys())
  {
    ASSERT_NEAR(arma::abs(testMultiMat2(key) - testMultiMat(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMM("testMultiMat")(key) - testMultiMat(key)).max(), 0.0, 1e-13);
  }

  Multi<arma::cube> testMultiCube2;
  d.get(testMultiCube2, "testMultiCube");
  for(auto &key : testMultiCube.getKeys())
  {
    ASSERT_NEAR(arma::abs(testMultiCube2(key) - testMultiCube(key)).max(), 0.0, 1e-13);
    ASSERT_NEAR(arma::abs(d.getMC("testMultiCube")(key) - testMultiCube(key)).max(), 0.0, 1e-13);
  }
}

//==============================================================================


TEST(DataTree, isEmpty)
{
  DataTree d;
  ASSERT_TRUE(d.isEmpty());

  arma::vec testVec = arma::vec();
  d.set("testVec", testVec);
  ASSERT_FALSE(d.isEmpty());
}

//==============================================================================


TEST(DataTree, size)
{
  DataTree d;
  ASSERT_EQ(d.size(), 0);

  arma::vec testVec = arma::randu(12);
  d.set("testVec", testVec);
  ASSERT_EQ(d.size(), 1);

  Multi<IMAT > testMultiMat;
  testMultiMat(1) = arma::randi<IMAT >(12, 8);
  testMultiMat(2) = arma::randi<IMAT >(12, 8);
  d.set("testMultiMat", testMultiMat);

  Multi<ICUBE > testMultiCube;
  testMultiCube(1) = arma::randi<ICUBE >(12, 2, 7);
  testMultiCube(7) = arma::randi<ICUBE >(12, 2, 7);
  d.set("testMultiCube", testMultiCube);

  d.set("testMultiIvec", Multi<IVEC>());

  ASSERT_EQ(d.size(), 4);
}


//==============================================================================

TEST(DataTree, clear)
{
  DataTree d;
  arma::vec testVec = arma::randu(12);
  d.set("testVec", testVec);
  ASSERT_FALSE(d.isEmpty());

  d.clear();
  ASSERT_TRUE(d.isEmpty());
}


//==============================================================================

TEST(DataTree, clean)
{
  DataTree d1;
  arma::vec testVec = arma::randu(12);
  arma::vec testVec2 = arma::randu(12);
  d1.set("vec1", testVec);
  d1.set("vec2", arma::vec());
  d1.set("vec3", testVec2);
  d1.set("mat", arma::mat());
  DataTree d2 = d1;

  ASSERT_EQ(d1, d2);

  d1.clean();
  ASSERT_NE(d1, d2);
}


//==============================================================================

TEST(DataTree, merge)
{
  arma::vec v1 = arma::randu(12);
  arma::vec v2 = arma::randu(12);
  std::string s = "plop";

  // DataTree 1
  DataTree d1;
  d1.set("v1", v1);
  d1.set("myString", s);

  // DataTree 2
  DataTree d2;
  d2.set("v2", v2);

  // DataTree Result
  DataTree resu;
  resu.set("v1", v1);
  resu.set("v2", v2);
  resu.set("myString", s);


  d1.merge(d2);
  ASSERT_EQ(d1, resu);

  // overwritting test (same key, different value)
  DataTree d22;
  d22.set("v2", v1);
  d2.merge(d22);
  ASSERT_EQ(d2, d22);
}


//==============================================================================

TEST(DataTree, operatorEquals)
{
  DataTree d;
  d.set("INT", 9);
  d.set("double", 9.8);
  d.set("string", "test");
  d.set("vec"  , arma::randu(70     ).eval());
  d.set("mat"  , arma::randu(70,7   ).eval());
  d.set("cube" , arma::randu(70,9,24).eval());
  d.set("ivec" , arma::randi<IVEC >(70     ));
  d.set("imat" , arma::randi<IMAT >(70,7   ));
  d.set("icube", arma::randi<ICUBE>(70,9,24));

  // TODO: add tests of uvec, umat, ucube

  Multi<arma::vec> multiV;
  multiV(3,8,9) = arma::randu(5);
  multiV(8,9,3) = arma::randu(3);

  Multi<IVEC> multiIV;
  multiIV(8) = arma::randi<IVEC>(3);

  Multi<arma::mat> multiM;
  multiM(8,1) = arma::randu(3,8);

  Multi<IMAT> multiIM;
  multiIM(0,7) = arma::randi<IMAT>(3,8);

  multiIM(13,-5,4) = IMAT();
  Multi<arma::cube> multiC;

  multiC(4,7,2) = arma::cube();
  Multi<ICUBE> multiIC;

  multiIC(4,4,4) = arma::randi<ICUBE>(4,4,4);

  d.set("multiVec"  , multiV );
  d.set("multiIVec" , multiIV);
  d.set("multiMat"  , multiM );
  d.set("multiIMat" , multiIM);
  d.set("multiCube" , multiC );
  d.set("multiICube", multiIC);

  DataTree copy = d;
  ASSERT_EQ(d, copy);

  DataTree d2;
  d2.set("multiVec"  , multiV );
  d2.set("multiIVec" , multiIV);
  d2.set("multiMat"  , multiM );
  d2.set("multiIMat" , multiIM);
  d2.set("multiCube" , multiC );
  d2.set("multiICube", multiIC);

  DataTree d2Copy;
  Multi<arma::vec> multiV2;
  multiV2(8,9,3) = multiV(8,9,3);
  multiV2(3,8,9) = multiV(3,8,9);

  Multi<IMAT> multiIM2;
  multiIM2(13,-5,4) = multiIM(13,-5,4);
  multiIM2(0,7) = multiIM(0,7);

  d2Copy.set("multiVec"  , multiV2 );
  d2Copy.set("multiIVec" , multiIV );
  d2Copy.set("multiMat"  , multiM  );
  d2Copy.set("multiIMat" , multiIM2);
  d2Copy.set("multiCube" , multiC  );
  d2Copy.set("multiICube", multiIC );

  ASSERT_EQ(d2, d2Copy);

  d2Copy.set("emptyMultiVec", Multi<arma::vec>());
  ASSERT_NE(d2, d2Copy);

  DataTree d3;
  arma::vec v1 = arma::randu(7);
  arma::vec v2 = arma::randu(76);
  arma::vec v3 = arma::randu(8);
  d3.set("v1", v1);
  d3.set("v2", v2);
  d3.set("v3", v3);

  DataTree d3Copy;
  d3Copy.set("v3", v3);
  d3Copy.set("v1", v1);
  d3Copy.set("v2", v2);

  ASSERT_EQ(d3, d3Copy);
}

//==============================================================================

TEST(DataTree, getSection)
{
  DataTree dt;
  dt.set("test0/val0", 5);
  dt.set("test0/val1", "plop");
  dt.set("test0/val2/b", 6);
  dt.set("test0/val2/a", 7);
  dt.set("test0/val3", 8);
  dt.set("test0/val4", arma::mat());

  DataTree dt2 = dt.getSection("test0/");
  ASSERT_EQ(dt2.getS("val1"), "plop");

  DataTree dt3 = dt.getSection("test0/val2/");
  ASSERT_EQ(dt3.getI("a"), 7);
}

//==============================================================================

TEST(DataTree, getDefault)
{
  auto dt0 = DataTree::getDefault();
  auto dt1 = DataTree();
  dt1.setDefault();

  ASSERT_EQ(dt0, dt1);
}

