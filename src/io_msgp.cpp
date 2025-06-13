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

#include "io_msgp.h"
#include "tools.h"
#include <sstream>
#include "msgpack.hpp"
#include "datatree.h"
#include "gzstream.h"

/** \file
 *  \brief Methods of the IOmsgp class.
 *
 *  The way we use the msgpack format to store our data may be changed in the
 *  near future. Please wait before stealing this class for your own project if
 *  you want to preserve interoperability.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 *  Create an empty DataTree.
 */

IOmsgp::IOmsgp(void)
{
  DBG_ENTER;
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree object from MSGP formatted content.
 */

DataTree IOmsgp::fromContent(const std::string &content)
{
  DBG_ENTER;

  DataTree result;

  msgpack::object_handle oh = msgpack::unpack(content.data(), content.size());
  msgpack::object obj = oh.get();
  std::map<std::string, msgpack::object> map = obj.convert();

  for (auto &k : map)
  {
    try
    {
      auto val = k.second.as<INT>();
      result.set(k.first, val);
    }
    catch (...) { }

    try
    {
      if (!result.contains(k.first))
      {
        auto val = k.second.as<double>();
        result.set(k.first, val);
      }
    }
    catch (...) { }

    try
    {
      if (!result.contains(k.first))
      {
        auto val = k.second.as<std::string>();
        result.set(k.first, val);
      }
    }
    catch (...) { }

    try
    {
      if (!result.contains(k.first))
      {
        auto val = k.second.as<std::tuple<std::string, std::string> >();
        std::string dataType = std::get<0>(val);

        if (dataType == "IV")
        {
          result.set(k.first, array2ivec(std::get<1>(val)));
        }
        if (dataType == "IM")
        {
          result.set(k.first, array2imat(std::get<1>(val)));
        }
        if (dataType == "IC")
        {
          result.set(k.first, array2icube(std::get<1>(val)));
        }
        if (dataType == "UV")
        {
          result.set(k.first, array2uvec(std::get<1>(val)));
        }
        if (dataType == "UM")
        {
          result.set(k.first, array2umat(std::get<1>(val)));
        }
        if (dataType == "UC")
        {
          result.set(k.first, array2ucube(std::get<1>(val)));
        }
        if (dataType == "V")
        {
          result.set(k.first, array2vec(std::get<1>(val)));
        }
        if (dataType == "M")
        {
          result.set(k.first, array2mat(std::get<1>(val)));
        }
        if (dataType == "C")
        {
          result.set(k.first, array2cube(std::get<1>(val)));
        }
      }
    }
    catch (...) { }

    try
    {
      if (!result.contains(k.first))
      {
        auto val = k.second.as<std::tuple<std::string, std::map<std::vector<INT>, std::string> > >();
        std::string dataType = std::get<0>(val);

        if (dataType == "UV")
        {
          result.set(k.first, array2multiUVec(std::get<1>(val)));
        }
        if (dataType == "UM")
        {
          result.set(k.first, array2multiUMat(std::get<1>(val)));
        }
        if (dataType == "UC")
        {
          result.set(k.first, array2multiUCube(std::get<1>(val)));
        }
        if (dataType == "IV")
        {
          result.set(k.first, array2multiIVec(std::get<1>(val)));
        }
        if (dataType == "IM")
        {
          result.set(k.first, array2multiIMat(std::get<1>(val)));
        }
        if (dataType == "IC")
        {
          result.set(k.first, array2multiICube(std::get<1>(val)));
        }
        if (dataType == "V")
        {
          result.set(k.first, array2multiVec(std::get<1>(val)));
        }
        if (dataType == "M")
        {
          result.set(k.first, array2multiMat(std::get<1>(val)));
        }
        if (dataType == "C")
        {
          result.set(k.first, array2multiCube(std::get<1>(val)));
        }
      }
    }
    catch (...) { }

  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Write a DataTree instance.
 */

void IOmsgp::saveDataTree(const DataTree &d, const std::string &filename)
{
  DBG_ENTER;

  //  std::ofstream fp(filename.c_str());
  ogzstream fp(filename.c_str());
  std::string str = serializeDataTree(d);
  fp.write(str.c_str(), str.size());
  fp.close();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Serialize a DataTree instance.
 */

std::string IOmsgp::serializeDataTree(const DataTree &d)
{
  DBG_ENTER;

  msgpack::sbuffer buffer;
  msgpack::packer<msgpack::sbuffer> pk(&buffer);

  //d.clean();
  INT size = (INT)(d.size() - d.emptyKeys.size());

  pk.pack_map(size);

  for (auto &it : d.intMap)
  {
    pk.pack(it.first);
    pk.pack(it.second);
  }

  for (auto &it : d.doubleMap)
  {
    pk.pack(it.first);
    pk.pack(it.second);
  }

  for (auto &it : d.stringMap)
  {
    pk.pack(it.first);
    pk.pack(it.second);
  }

  for (auto &it : d.vecMap)
  {
    pk.pack(it.first);
    pk.pack(vec2array(it.second));
  }

  for (auto &it : d.matMap)
  {
    pk.pack(it.first);
    pk.pack(mat2array(it.second));
  }

  for (auto &it : d.cubeMap)
  {
    pk.pack(it.first);
    pk.pack(cube2array(it.second));
  }

  for (auto &it : d.uvecMap)
  {
    pk.pack(it.first);
    pk.pack(uvec2array(it.second));
  }

  for (auto &it : d.umatMap)
  {
    pk.pack(it.first);
    pk.pack(umat2array(it.second));
  }

  for (auto &it : d.ucubeMap)
  {
    pk.pack(it.first);
    pk.pack(ucube2array(it.second));
  }

  for (auto &it : d.multiUVecMap)
  {
    pk.pack(it.first);
    pk.pack(multiUVec2array(it.second));
  }

  for (auto &it : d.multiUMatMap)
  {
    pk.pack(it.first);
    pk.pack(multiUMat2array(it.second)); 
  }

  for (auto &it : d.multiUCubeMap)
  {
    pk.pack(it.first);
    pk.pack(multiUCube2array(it.second));
  }

  for (auto &it : d.ivecMap)
  {
    pk.pack(it.first);
    pk.pack(ivec2array(it.second));
  }

  for (auto &it : d.imatMap)
  {
    pk.pack(it.first);
    pk.pack(imat2array(it.second));
  }

  for (auto &it : d.icubeMap)
  {
    pk.pack(it.first);
    pk.pack(icube2array(it.second));
  }

  for (auto &it : d.multiIVecMap)
  {
    pk.pack(it.first);
    pk.pack(multiIVec2array(it.second));
  }

  for (auto &it : d.multiIMatMap)
  {
    pk.pack(it.first);
    pk.pack(multiIMat2array(it.second)); 
  }

  for (auto &it : d.multiICubeMap)
  {
    pk.pack(it.first);
    pk.pack(multiICube2array(it.second));
  }

  for (auto &it : d.multiVecMap)
  {
    pk.pack(it.first);
    pk.pack(multiVec2array(it.second));
  }

  for (auto &it : d.multiMatMap)
  {
    pk.pack(it.first);
    pk.pack(multiMat2array(it.second)); 
  }

  for (auto &it : d.multiCubeMap)
  {
    pk.pack(it.first);
    pk.pack(multiCube2array(it.second));
  }

  DBG_RETURN(std::string(buffer.data(), buffer.size()));
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an arma::vec to an arma::vec.
 */

arma::vec IOmsgp::array2vec(const std::string &v)
{
  arma::vec m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an arma::mat to an arma::mat.
 */

arma::mat IOmsgp::array2mat(const std::string &v)
{
  arma::mat m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an arma::cube to an arma::cube.
 */

arma::cube IOmsgp::array2cube(const std::string &v)
{
  arma::cube m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an arma::vec to a string.
 */

std::tuple<std::string, std::string> IOmsgp::vec2array(const arma::vec &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("V", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an arma::mat to a string.
 */

std::tuple<std::string, std::string> IOmsgp::mat2array(const arma::mat &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("M", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an arma::cube to a string.
 */

std::tuple<std::string, std::string> IOmsgp::cube2array(const arma::cube &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("C", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an IVEC to an IVEC.
 */

IVEC IOmsgp::array2ivec(const std::string &v)
{
  IVEC m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an IMAT to an IMAT.
 */

IMAT IOmsgp::array2imat(const std::string &v)
{
  IMAT m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an ICUBE to an ICUBE.
 */

ICUBE IOmsgp::array2icube(const std::string &v)
{
  ICUBE m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an UVEC to an UVEC.
 */

UVEC IOmsgp::array2uvec(const std::string &v)
{
  UVEC m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an UMAT to an UMAT.
 */

UMAT IOmsgp::array2umat(const std::string &v)
{
  UMAT m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert a string contening the data of an UCUBE to an UCUBE.
 */

UCUBE IOmsgp::array2ucube(const std::string &v)
{
  UCUBE m;

  std::istringstream iss(v);
  m.load(iss, arma::arma_binary);

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an IVEC to a string.
 */

std::tuple<std::string, std::string> IOmsgp::ivec2array(const IVEC &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("IV", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an IMAT to a string.
 */

std::tuple<std::string, std::string> IOmsgp::imat2array(const IMAT &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("IM", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an ICUBE to a string.
 */

std::tuple<std::string, std::string> IOmsgp::icube2array(const ICUBE &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("IC", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an UVEC to a string.
 */

std::tuple<std::string, std::string> IOmsgp::uvec2array(const UVEC &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("UV", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an UMAT to a string.
 */

std::tuple<std::string, std::string> IOmsgp::umat2array(const UMAT &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("UM", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an UCUBE to a string.
 */

std::tuple<std::string, std::string> IOmsgp::ucube2array(const UCUBE &v)
{
  std::ostringstream oss;
  v.save(oss, arma::arma_binary);

  return std::make_tuple("UC", oss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<IVEC>.
 */

Multi<IVEC > IOmsgp::array2multiIVec(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<IVEC > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    IVEC vec;
    vec.load(iss, arma::arma_binary);
    m(key) = vec;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<IMAT>.
 */

Multi<IMAT > IOmsgp::array2multiIMat(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<IMAT > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    IMAT mat;
    mat.load(iss, arma::arma_binary);
    m(key) = mat;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<ICUBE>.
 */

Multi<ICUBE > IOmsgp::array2multiICube(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<ICUBE > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    ICUBE cube;
    cube.load(iss, arma::arma_binary);
    m(key) = cube;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<UVEC>.
 */

Multi<UVEC > IOmsgp::array2multiUVec(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<UVEC > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    UVEC vec;
    vec.load(iss, arma::arma_binary);
    m(key) = vec;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<UMAT>.
 */

Multi<UMAT > IOmsgp::array2multiUMat(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<UMAT > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    UMAT mat;
    mat.load(iss, arma::arma_binary);
    m(key) = mat;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<UCUBE>.
 */

Multi<UCUBE > IOmsgp::array2multiUCube(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<UCUBE > m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    UCUBE cube;
    cube.load(iss, arma::arma_binary);
    m(key) = cube;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<IVEC> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiIVec2array(const Multi<IVEC > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("IV", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<IMAT> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiIMat2array(const Multi<IMAT > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("IM", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<ICUBE> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiICube2array(const Multi<ICUBE > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("IC", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<UVEC> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiUVec2array(const Multi<UVEC > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("UV", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<UMAT> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiUMat2array(const Multi<UMAT > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("UM", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<ICUBE> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiUCube2array(const Multi<UCUBE > &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("UC", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<arma::vec>.
 */

Multi<arma::vec> IOmsgp::array2multiVec(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<arma::vec> m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    arma::vec vec;
    vec.load(iss, arma::arma_binary);
    m(key) = vec;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<arma::mat>.
 */

Multi<arma::mat> IOmsgp::array2multiMat(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<arma::mat> m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    arma::mat mat;
    mat.load(iss, arma::arma_binary);
    m(key) = mat;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an array to an Multi<arma::cube>.
 */

Multi<arma::cube> IOmsgp::array2multiCube(const std::map<std::vector<INT>, std::string> &v)
{
  Multi<arma::cube> m;

  for(auto it = v.begin(); it != v.end(); ++it) 
  {
    // key
    std::vector<INT> key(it->first.begin(), it->first.end());

    // value
    std::istringstream iss(it->second);
    arma::cube cube;
    cube.load(iss, arma::arma_binary);
    m(key) = cube;
  }

  return m;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<arma::vec> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiVec2array(const Multi<arma::vec> &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("V", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<arma::mat> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiMat2array(const Multi<arma::mat> &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("M", result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert an Multi<arma::cube> to an array.
 */

std::tuple<std::string, std::map<std::vector<INT>, std::string> > IOmsgp::multiCube2array(const Multi<arma::cube> &v)
{
  std::map<std::vector<INT>, std::string> result;

  for (auto &key : v.getKeys())
  {
    std::ostringstream oss;
    v(key).save(oss, arma::arma_binary);
    std::string binary_data = oss.str();

    std::vector<INT> myVector(key.begin(), key.end());
    result.emplace(myVector, binary_data);
  }

  return std::make_tuple("C", result);
}
