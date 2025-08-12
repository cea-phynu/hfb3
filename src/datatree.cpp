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

#include <regex>
#include "datatree.h"
#include "tools.h"
#include "io_amedee.h"
#include "io_berger.h"
#include "io_hfb3.h"
#include "io_msgp.h"
#include "io_json.h"
#include "general.h"
#include "io_hfb3.h"

/** \file
 *  \brief Methods of the DataTree class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/// Helper macro to insert keys matching a pattern to an std::map.
#define FIND_AND_INSERT(K, M) for (auto &K : M) if (K.first.find(pattern) != std::string::npos) result.insert(K.first)

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor with or without a filename to load.
 */

DataTree::DataTree(const std::string &filename)
{
  general.setGlobalValidKeys();

  if (filename != "")
  {
    (*this) = fromContent(Tools::readFile(filename));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Save a DataTree to a .msg.gz file.

void DataTree::save(const std::string &filename, bool verbose) const
{
  IOmsgp().saveDataTree(*this, filename);
  if (verbose) Tools::mesg("Saving", "DataTree saved in " + PF_GREEN(filename));
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Return the number of keys.

UINT DataTree::size(void) const
{
  UINT size = 0;
  size += boolMap.size();
  size += intMap.size();
  size += doubleMap.size();
  size += stringMap.size();
  size += vecMap.size();
  size += matMap.size();
  size += cubeMap.size();
  size += uvecMap.size();
  size += umatMap.size();
  size += ucubeMap.size();
  size += ivecMap.size();
  size += imatMap.size();
  size += icubeMap.size();
  size += multiIVecMap.size();
  size += multiIMatMap.size();
  size += multiICubeMap.size();
  size += multiUVecMap.size();
  size += multiUMatMap.size();
  size += multiUCubeMap.size();
  size += multiVecMap.size();
  size += multiMatMap.size();
  size += multiCubeMap.size();
  size += emptyKeys.size();
  return size;
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Return a DataTree from a file content.

DataTree DataTree::fromContent(const std::string &content)
{
  DBG_ENTER;

  std::string deflated = content;

  // Decompress if needed
  try
  {
    deflated = Tools::decompressString(content);
  }
  catch (...)
  {
  }

  DataTree d = IOberger().fromContent(deflated);

  if (d.isEmpty()) d = IOjson().fromContent(deflated);

  if (d.isEmpty()) d = IOamedee().fromContent(deflated);

  if (d.isEmpty()) d = IOhfb3().fromContent(deflated);

  if (d.isEmpty()) d = IOmsgp().fromContent(deflated);

  ASSERT(!d.isEmpty(), "content format error");

  DBG_RETURN(d);
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Check if the DataTree is empty.

bool DataTree::isEmpty(void) const
{
  return ( size() == 0 );
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Clear the tree.

void DataTree::clear(void)
{
  boolMap.clear();
  intMap.clear();
  doubleMap.clear();
  stringMap.clear();
  vecMap.clear();
  matMap.clear();
  cubeMap.clear();
  ivecMap.clear();
  imatMap.clear();
  icubeMap.clear();
  uvecMap.clear();
  umatMap.clear();
  ucubeMap.clear();
  multiIVecMap.clear();
  multiIMatMap.clear();
  multiICubeMap.clear();
  multiUVecMap.clear();
  multiUMatMap.clear();
  multiUCubeMap.clear();
  multiVecMap.clear();
  multiMatMap.clear();
  multiCubeMap.clear();
  emptyKeys.clear();
}

//==============================================================================
//==============================================================================
//==============================================================================

/// Clean the tree (remove all empty arma object).

void DataTree::clean(void)
{
  for (auto it = vecMap.begin()       ;it != vecMap.end()       ;) { if (it->second.empty()) it = vecMap.erase(it)       ;else (it++);}
  for (auto it = matMap.begin()       ;it != matMap.end()       ;) { if (it->second.empty()) it = matMap.erase(it)       ;else (it++);}
  for (auto it = cubeMap.begin()      ;it != cubeMap.end()      ;) { if (it->second.empty()) it = cubeMap.erase(it)      ;else (it++);}
  for (auto it = ivecMap.begin()      ;it != ivecMap.end()      ;) { if (it->second.empty()) it = ivecMap.erase(it)      ;else (it++);}
  for (auto it = imatMap.begin()      ;it != imatMap.end()      ;) { if (it->second.empty()) it = imatMap.erase(it)      ;else (it++);}
  for (auto it = icubeMap.begin()     ;it != icubeMap.end()     ;) { if (it->second.empty()) it = icubeMap.erase(it)     ;else (it++);}
  for (auto it = uvecMap.begin()      ;it != uvecMap.end()      ;) { if (it->second.empty()) it = uvecMap.erase(it)      ;else (it++);}
  for (auto it = umatMap.begin()      ;it != umatMap.end()      ;) { if (it->second.empty()) it = umatMap.erase(it)      ;else (it++);}
  for (auto it = ucubeMap.begin()     ;it != ucubeMap.end()     ;) { if (it->second.empty()) it = ucubeMap.erase(it)     ;else (it++);}

  for (auto it = multiVecMap.begin()  ;it != multiVecMap.end()  ;) { if (it->second.empty()) it = multiVecMap.erase(it)  ;else (it++);}
  for (auto it = multiMatMap.begin()  ;it != multiMatMap.end()  ;) { if (it->second.empty()) it = multiMatMap.erase(it)  ;else (it++);}
  for (auto it = multiCubeMap.begin() ;it != multiCubeMap.end() ;) { if (it->second.empty()) it = multiCubeMap.erase(it) ;else (it++);}
  for (auto it = multiUVecMap.begin() ;it != multiUVecMap.end() ;) { if (it->second.empty()) it = multiUVecMap.erase(it) ;else (it++);}
  for (auto it = multiUMatMap.begin() ;it != multiUMatMap.end() ;) { if (it->second.empty()) it = multiUMatMap.erase(it) ;else (it++);}
  for (auto it = multiUCubeMap.begin();it != multiUCubeMap.end();) { if (it->second.empty()) it = multiUCubeMap.erase(it);else (it++);}
  for (auto it = multiIVecMap.begin() ;it != multiIVecMap.end() ;) { if (it->second.empty()) it = multiIVecMap.erase(it) ;else (it++);}
  for (auto it = multiIMatMap.begin() ;it != multiIMatMap.end() ;) { if (it->second.empty()) it = multiIMatMap.erase(it) ;else (it++);}
  for (auto it = multiICubeMap.begin();it != multiICubeMap.end();) { if (it->second.empty()) it = multiICubeMap.erase(it);else (it++);} }

//==============================================================================
//==============================================================================
//==============================================================================

/** Delete a value from the data tree.
 *
 *  \param key The key.
 *  \return the number of deleted values.
 */

INT DataTree::del(const std::string &key)
{
  INT nbDeleted = 0;

  if (key[key.size() - 1] == '/')
  {
    for (auto &i : boolMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        boolMap.erase(i.first);
      }
    }

    for (auto &i : intMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        intMap.erase(i.first);
      }
    }

    for (auto &i : doubleMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        doubleMap.erase(i.first);
      }
    }

    for (auto &i : stringMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        stringMap.erase(i.first);
      }
    }

    for (auto &i : vecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        vecMap.erase(i.first);
      }
    }

    for (auto &i : matMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        matMap.erase(i.first);
      }
    }

    for (auto &i : cubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        cubeMap.erase(i.first);
      }
    }

    for (auto &i : uvecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        uvecMap.erase(i.first);
      }
    }

    for (auto &i : umatMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        umatMap.erase(i.first);
      }
    }

    for (auto &i : ucubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        ucubeMap.erase(i.first);
      }
    }

    for (auto &i : ivecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        ivecMap.erase(i.first);
      }
    }

    for (auto &i : imatMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        imatMap.erase(i.first);
      }
    }

    for (auto &i : icubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        icubeMap.erase(i.first);
      }
    }

    for (auto &i : multiUVecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiUVecMap.erase(i.first);
      }
    }

    for (auto &i : multiUMatMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiUMatMap.erase(i.first);
      }
    }

    for (auto &i : multiUCubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiUCubeMap.erase(i.first);
      }
    }

    for (auto &i : multiIVecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiIVecMap.erase(i.first);
      }
    }

    for (auto &i : multiIMatMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiIMatMap.erase(i.first);
      }
    }

    for (auto &i : multiICubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiICubeMap.erase(i.first);
      }
    }

    for (auto &i : multiVecMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiVecMap.erase(i.first);
      }
    }

    for (auto &i : multiMatMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiMatMap.erase(i.first);
      }
    }

    for (auto &i : multiCubeMap)
    {
      if (i.first.compare(0, key.size(), key) == 0)
      {
        nbDeleted++;
        multiCubeMap.erase(i.first);
      }
    }
  }
  else
  {
    if (boolMap.count(key) == 1)
    {
      nbDeleted++;
      boolMap.erase(key);
    }

    if (intMap.count(key) == 1)
    {
      nbDeleted++;
      intMap.erase(key);
    }

    if (doubleMap.count(key) == 1)
    {
      nbDeleted++;
      doubleMap.erase(key);
    }

    if (stringMap.count(key) == 1)
    {
      nbDeleted++;
      stringMap.erase(key);
    }

    if (vecMap.count(key) == 1)
    {
      nbDeleted++;
      vecMap.erase(key);
    }

    if (matMap.count(key) == 1)
    {
      nbDeleted++;
      matMap.erase(key);
    }

    if (cubeMap.count(key) == 1)
    {
      nbDeleted++;
      cubeMap.erase(key);
    }

    if (uvecMap.count(key) == 1)
    {
      nbDeleted++;
      uvecMap.erase(key);
    }

    if (umatMap.count(key) == 1)
    {
      nbDeleted++;
      umatMap.erase(key);
    }

    if (ucubeMap.count(key) == 1)
    {
      nbDeleted++;
      ucubeMap.erase(key);
    }

    if (multiUVecMap.count(key) == 1)
    {
      nbDeleted++;
      multiUVecMap.erase(key);
    }

    if (multiUMatMap.count(key) == 1)
    {
      nbDeleted++;
      multiUMatMap.erase(key);
    }

    if (multiUCubeMap.count(key) == 1)
    {
      nbDeleted++;
      multiUCubeMap.erase(key);
    }

    if (ivecMap.count(key) == 1)
    {
      nbDeleted++;
      ivecMap.erase(key);
    }

    if (imatMap.count(key) == 1)
    {
      nbDeleted++;
      imatMap.erase(key);
    }

    if (icubeMap.count(key) == 1)
    {
      nbDeleted++;
      icubeMap.erase(key);
    }

    if (multiIVecMap.count(key) == 1)
    {
      nbDeleted++;
      multiIVecMap.erase(key);
    }

    if (multiIMatMap.count(key) == 1)
    {
      nbDeleted++;
      multiIMatMap.erase(key);
    }

    if (multiICubeMap.count(key) == 1)
    {
      nbDeleted++;
      multiICubeMap.erase(key);
    }

    if (multiVecMap.count(key) == 1)
    {
      nbDeleted++;
      multiVecMap.erase(key);
    }

    if (multiMatMap.count(key) == 1)
    {
      nbDeleted++;
      multiMatMap.erase(key);
    }

    if (multiCubeMap.count(key) == 1)
    {
      nbDeleted++;
      multiCubeMap.erase(key);
    }
  }

  return nbDeleted;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, std::string) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::set(const std::string &key, const char *val)
{
  set(key, std::string(val));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, std::string) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::set(const std::string &key, const std::string &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "S"), "Invalid DataTree: wrong type for key: '" + key + "'");

  stringMap[key] = std::regex_replace(val, std::regex("\n"), " ");
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Delete a key or create an empty key.
 *
 *  \param key The key.
 */

void DataTree::set(const std::string &key)
{
  del(key);
  emptyKeys.insert(key);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, bool) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::set(const std::string &key, const bool &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "B"), "Invalid DataTree: wrong type for key: '" + key + "'");

  boolMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, INT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::set(const std::string &key, INT val)
{
  if (strict_mode)
    ASSERT(isValid(key, "I"), "Invalid DataTree: wrong type for key: '" + key + "'");

  intMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, double) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::set(const std::string &key, double val)
{
  if (strict_mode)
    ASSERT(isValid(key, "D"), "Invalid DataTree: wrong type for key: '" + key + "'");

  doubleMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::vec) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::vec object.
 */

void DataTree::set(const std::string &key, const arma::vec &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "V"), "Invalid DataTree: wrong type for key: '" + key + "'");

  vecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::mat) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::mat object.
 */

void DataTree::set(const std::string &key, const arma::mat &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "M"), "Invalid DataTree: wrong type for key: '" + key + "'");

  matMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::cube) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::cube object.
 */

void DataTree::set(const std::string &key, const arma::cube &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "C"), "Invalid DataTree: wrong type for key: '" + key + "'");

  cubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, IVEC) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The IVEC object.
 */

void DataTree::set(const std::string &key, const IVEC &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "IV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  ivecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, IMAT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The IMAT object.
 */

void DataTree::set(const std::string &key, const IMAT &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "IM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  imatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, ICUBE) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The ICUBE object.
 */

void DataTree::set(const std::string &key, const ICUBE &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "IC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  icubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<IVEC >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<IVEC > object.
 */

void DataTree::set(const std::string &key, const Multi<IVEC> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MIV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiIVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<IMAT >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<IMAT > object.
 */

void DataTree::set(const std::string &key, const Multi<IMAT> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MIM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiIMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<ICUBE>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<ICUBE> object.
 */

void DataTree::set(const std::string &key, const Multi<ICUBE> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MIC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiICubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UVEC) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UVEC object.
 */

void DataTree::set(const std::string &key, const UVEC &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "UV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  uvecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UMAT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UMAT object.
 */

void DataTree::set(const std::string &key, const UMAT &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "UM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  umatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UCUBE) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UCUBE object.
 */

void DataTree::set(const std::string &key, const UCUBE &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "UC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  ucubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UVEC >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UVEC > object.
 */

void DataTree::set(const std::string &key, const Multi<UVEC> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MUV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiUVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UMAT >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UMAT > object.
 */

void DataTree::set(const std::string &key, const Multi<UMAT> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MUM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiUMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UCUBE>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UCUBE> object.
 */

void DataTree::set(const std::string &key, const Multi<UCUBE> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MUC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiUCubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::vec>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::vec> object.
 */

void DataTree::set(const std::string &key, const Multi<arma::vec> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::mat>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::m  t> object.
 */

void DataTree::set(const std::string &key, const Multi<arma::mat> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::cube>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::cube> object.
 */

void DataTree::set(const std::string &key, const Multi<arma::cube> &val)
{
  if (strict_mode)
    ASSERT(isValid(key, "MC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  multiCubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an INT value from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(INT &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "I"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (intMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no int key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = intMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get a bool value from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(bool &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "B"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (boolMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no bool key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = boolMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get a double value from the tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(double &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "D"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (doubleMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no double key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = doubleMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an std::string from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(std::string &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "S"), "Invalid DataTree: wrong type for key: '" + key + "'");


  if (stringMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no std::string key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = stringMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::vec from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(arma::vec &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "V"), "Invalid DataTree: wrong type for key: '" + key + "'");


  if (vecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no arma::vec key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = vecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::mat from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(arma::mat &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "M"), "Invalid DataTree: wrong type for key: '" + key + "'");


  if (matMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no arma::mat key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = matMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::cube from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(arma::cube &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "C"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (cubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no arma::cube key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = cubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UVEC from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(UVEC &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "UV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (uvecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no UVEC key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = uvecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UMAT from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(UMAT &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "UM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (umatMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no UMAT key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = umatMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UCUBE from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(UCUBE &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "UC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (ucubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no UCUBE key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = ucubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UVEC > from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<UVEC > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MUV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiUVecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<UVEC> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiUVecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UMAT > from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<UMAT > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MUM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiUMatMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<UMAT> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiUMatMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UCUBE> from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<UCUBE > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MUC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiUCubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<UCUBE> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiUCubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an IVEC from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(IVEC &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "IV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (ivecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no IVEC key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = ivecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an IMAT from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(IMAT &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "IM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (imatMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no IMAT key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = imatMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an ICUBE from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(ICUBE &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "IC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (icubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no ICUBE key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = icubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<IVEC > from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<IVEC > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MIV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiIVecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<IVEC> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiIVecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<IMAT > from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<IMAT > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MIM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiIMatMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<IMAT> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiIMatMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<ICUBE> from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<ICUBE > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MIC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiICubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<ICUBE> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiICubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::vec> from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<arma::vec> &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MV"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiVecMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<arma::vec> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiVecMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::mat> from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<arma::mat> &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MM"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiMatMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<arma::mat> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiMatMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::cube> from the data tree.
 *
 *  \param target The target variable.
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 *  \return Returns true if the target variable has been set to the new value.
 */

bool DataTree::get(Multi<arma::cube > &target, const std::string &key, bool ignore) const
{
  DBG_ENTER;

  if (strict_mode)
    ASSERT(isValid(key, "MC"), "Invalid DataTree: wrong type for key: '" + key + "'");

  if (multiCubeMap.count(key) == 0)
  {
    if (!ignore)
    {
      ERROR("no Multi<arma::cube> key '" + key + "'");
    }

    DBG_RETURN(false);
  }

  target = multiCubeMap.find(key)->second;
  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get a bool from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

bool DataTree::getB(const std::string &key, bool ignore, bool *isFound) const
{
  if (strict_mode)
    ASSERT(isValid(key, "B"), "Invalid DataTree: wrong type for key: '" + key + "'");

  bool result = false;
  bool found = false;
  found = get(result, key, ignore);
  if (isFound != NULL) *isFound = found;
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an INT from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

INT DataTree::getI(const std::string &key, bool ignore) const
{
  INT result = 0;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get a double from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

double DataTree::getD(const std::string &key, bool ignore) const
{
  double result = 0.0;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an std::string from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

std::string DataTree::getS(const std::string &key, bool ignore) const
{
  std::string result = "";
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::vec from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

arma::vec DataTree::getV(const std::string &key, bool ignore) const
{
  arma::vec result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::mat from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

arma::mat DataTree::getM(const std::string &key, bool ignore) const
{
  arma::mat result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an arma::cube from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

arma::cube DataTree::getC(const std::string &key, bool ignore) const
{
  arma::cube result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UVEC from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

UVEC DataTree::getUV(const std::string &key, bool ignore) const
{
  UVEC result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UMAT from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

UMAT DataTree::getUM(const std::string &key, bool ignore) const
{
  UMAT result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an UCUBE from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

UCUBE DataTree::getUC(const std::string &key, bool ignore) const
{
  UCUBE result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UVEC > from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<UVEC > DataTree::getMUV(const std::string &key, bool ignore) const
{
  Multi<UVEC > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UMAT > from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<UMAT > DataTree::getMUM(const std::string &key, bool ignore) const
{
  Multi<UMAT > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<UCUBE> from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<UCUBE > DataTree::getMUC(const std::string &key, bool ignore) const
{
  Multi<UCUBE > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an IVEC from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

IVEC DataTree::getIV(const std::string &key, bool ignore) const
{
  IVEC result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an IMAT from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

IMAT DataTree::getIM(const std::string &key, bool ignore) const
{
  IMAT result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an ICUBE from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

ICUBE DataTree::getIC(const std::string &key, bool ignore) const
{
  ICUBE result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<IVEC > from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<IVEC > DataTree::getMIV(const std::string &key, bool ignore) const
{
  Multi<IVEC > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<IMAT > from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<IMAT > DataTree::getMIM(const std::string &key, bool ignore) const
{
  Multi<IMAT > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<ICUBE> from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<ICUBE > DataTree::getMIC(const std::string &key, bool ignore) const
{
  Multi<ICUBE > result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::vec> from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<arma::vec> DataTree::getMV(const std::string &key, bool ignore) const
{
  Multi<arma::vec> result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::mat> from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<arma::mat> DataTree::getMM(const std::string &key, bool ignore) const
{
  Multi<arma::mat> result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get an Multi<arma::cube> from the data tree.
 *
 *  \param key The key.
 *  \param ignore If true, print a warning if the key does not exist.
 */

Multi<arma::cube> DataTree::getMC(const std::string &key, bool ignore) const
{
  Multi<arma::cube> result;
  get(result, key, ignore);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string DataTree::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"int       ", Tools::infoStr((INT)intMap.size())},
      {"double    ", Tools::infoStr((INT)doubleMap.size())},
      {"string    ", Tools::infoStr((INT)stringMap.size())},
      {"vec       ", Tools::infoStr((INT)vecMap.size())},
      {"mat       ", Tools::infoStr((INT)matMap.size())},
      {"cube      ", Tools::infoStr((INT)cubeMap.size())},
      {"uvec      ", Tools::infoStr((INT)uvecMap.size())},
      {"umat      ", Tools::infoStr((INT)umatMap.size())},
      {"ucube     ", Tools::infoStr((INT)ucubeMap.size())},
      {"multiUVec ", Tools::infoStr((INT)multiUVecMap.size())},
      {"multiUMat ", Tools::infoStr((INT)multiUMatMap.size())},
      {"multiUCube", Tools::infoStr((INT)multiUCubeMap.size())},
      {"ivec      ", Tools::infoStr((INT)ivecMap.size())},
      {"imat      ", Tools::infoStr((INT)imatMap.size())},
      {"icube     ", Tools::infoStr((INT)icubeMap.size())},
      {"multiIVec ", Tools::infoStr((INT)multiIVecMap.size())},
      {"multiIMat ", Tools::infoStr((INT)multiIMatMap.size())},
      {"multiICube", Tools::infoStr((INT)multiICubeMap.size())},
      {"multiVec  ", Tools::infoStr((INT)multiVecMap.size())},
      {"multiMat  ", Tools::infoStr((INT)multiMatMap.size())},
      {"multiCube ", Tools::infoStr((INT)multiCubeMap.size())},
      {"empty ", PF("%d", emptyKeys.size())},
    }, isShort);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"DataTree", ""},
      {"int       ", mapToStr("INT   ", intMap)},
      {"double    ", mapToStr("double", doubleMap)},
      {"string    ", mapToStr("string", stringMap)},
      {"vec       ", mapToStr("vec   ", vecMap)},
      {"mat       ", mapToStr("mat   ", matMap)},
      {"cube      ", mapToStr("cube  ", cubeMap)},
      {"uvec      ", mapToStr("uvec  ", uvecMap)},
      {"umat      ", mapToStr("umat  ", umatMap)},
      {"ucube     ", mapToStr("ucube ", ucubeMap)},
      {"multiUVec ", multiMapToStr("multiUVec  ", multiUVecMap)},
      {"multiUMat ", multiMapToStr("multiUMat  ", multiUMatMap)},
      {"multiUCube", multiMapToStr("multiUCube ", multiUCubeMap)},
      {"ivec      ", mapToStr("ivec  ", ivecMap)},
      {"imat      ", mapToStr("imat  ", imatMap)},
      {"icube     ", mapToStr("icube ", icubeMap)},
      {"multiIVec ", multiMapToStr("multiIVec  ", multiIVecMap)},
      {"multiIMat ", multiMapToStr("multiIMat  ", multiIMatMap)},
      {"multiICube", multiMapToStr("multiICube ", multiICubeMap)},
      {"multiVec  ", multiMapToStr("multiVec  ", multiVecMap)},
      {"multiMat  ", multiMapToStr("multiMat  ", multiMatMap)},
      {"multiCube ", multiMapToStr("multiCube ", multiCubeMap)},
      {"empty ", PF("%d", emptyKeys.size())},
    }, isShort);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return an abstract of the dataTree.
 */

std::string DataTree::abstract(void) const
{
  std::stringstream result;

  if (intMap.size() != 0)
  {
    result << Tools::color("red") << "=== integer ===" << Tools::color() << std::endl;

    for (std::map<std::string, INT>::const_iterator it = intMap.begin(); it != intMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << it->second << std::endl;
  }

  if (doubleMap.size() != 0)
  {
    result << Tools::color("red") << "=== double ===" << Tools::color() << std::endl;

    for (std::map<std::string, double>::const_iterator it = doubleMap.begin(); it != doubleMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << it->second << std::endl;
  }

  if (stringMap.size() != 0)
  {
    result << Tools::color("red") << "=== string ===" << Tools::color() << std::endl;

    for (std::map<std::string, std::string>::const_iterator it = stringMap.begin(); it != stringMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << it->second << std::endl;
  }

  if (vecMap.size() != 0)
  {
    result << Tools::color("red") << "=== arma::vec ===" << Tools::color() << std::endl;

    for (std::map<std::string, arma::vec>::const_iterator it = vecMap.begin(); it != vecMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << ")" << std::endl;
  }

  if (matMap.size() != 0)
  {
    result << Tools::color("red") << "=== arma::mat ===" << Tools::color() << std::endl;

    for (std::map<std::string, arma::mat>::const_iterator it = matMap.begin(); it != matMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << ")" << std::endl;
  }

  if (cubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== arma::cube ===" << Tools::color() << std::endl;

    for (std::map<std::string, arma::cube>::const_iterator it = cubeMap.begin(); it != cubeMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << "x" << it->second.n_slices << ")" << std::endl;
  }

  if (uvecMap.size() != 0)
  {
    result << Tools::color("red") << "=== UVEC ===" << Tools::color() << std::endl;

    for (std::map<std::string, UVEC >::const_iterator it = uvecMap.begin(); it != uvecMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << ")" << std::endl;
  }

  if (umatMap.size() != 0)
  {
    result << Tools::color("red") << "=== UMAT ===" << Tools::color() << std::endl;

    for (std::map<std::string, UMAT >::const_iterator it = umatMap.begin(); it != umatMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << ")" << std::endl;
  }

  if (ucubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== UCUBE ===" << Tools::color() << std::endl;

    for (std::map<std::string, UCUBE >::const_iterator it = ucubeMap.begin(); it != ucubeMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << "x" << it->second.n_slices << ")" << std::endl;
  }

  if (multiUVecMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<UVEC > ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<UVEC > >::const_iterator it = multiUVecMap.begin(); it != multiUVecMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << ")" << std::endl;
      }
    }
  }

  if (multiUMatMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<UMAT > ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<UMAT > >::const_iterator it = multiUMatMap.begin(); it != multiUMatMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols << ")" << std::endl;
      }
    }
  }

  if (multiUCubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<UCUBE> ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<UCUBE > >::const_iterator it = multiUCubeMap.begin(); it != multiUCubeMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols  << "x" << it->second(key).n_slices << ")" << std::endl;
      }
    }
  }

  if (ivecMap.size() != 0)
  {
    result << Tools::color("red") << "=== IVEC ===" << Tools::color() << std::endl;

    for (std::map<std::string, IVEC >::const_iterator it = ivecMap.begin(); it != ivecMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << ")" << std::endl;
  }

  if (imatMap.size() != 0)
  {
    result << Tools::color("red") << "=== IMAT ===" << Tools::color() << std::endl;

    for (std::map<std::string, IMAT >::const_iterator it = imatMap.begin(); it != imatMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << ")" << std::endl;
  }

  if (icubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== ICUBE ===" << Tools::color() << std::endl;

    for (std::map<std::string, ICUBE >::const_iterator it = icubeMap.begin(); it != icubeMap.end(); ++it)
      result << Tools::color("blue") << it->first << ": " << Tools::color() << "(" << it->second.n_rows << "x" << it->second.n_cols << "x" << it->second.n_slices << ")" << std::endl;
  }

  if (multiIVecMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<IVEC > ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<IVEC > >::const_iterator it = multiIVecMap.begin(); it != multiIVecMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << ")" << std::endl;
      }
    }
  }

  if (multiIMatMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<IMAT > ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<IMAT > >::const_iterator it = multiIMatMap.begin(); it != multiIMatMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols << ")" << std::endl;
      }
    }
  }

  if (multiICubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<ICUBE> ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<ICUBE > >::const_iterator it = multiICubeMap.begin(); it != multiICubeMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols  << "x" << it->second(key).n_slices << ")" << std::endl;
      }
    }
  }

  if (multiVecMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<arma::vec> ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<arma::vec> >::const_iterator it = multiVecMap.begin(); it != multiVecMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << ")" << std::endl;
      }
    }
  }

  if (multiMatMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<arma::mat> ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<arma::mat> >::const_iterator it = multiMatMap.begin(); it != multiMatMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols << ")" << std::endl;
      }
    }
  }

  if (multiCubeMap.size() != 0)
  {
    result << Tools::color("red") << "=== Multi<arma::cube> ===" << Tools::color() << std::endl;

    for (std::map<std::string, Multi<arma::cube> >::const_iterator it = multiCubeMap.begin(); it != multiCubeMap.end(); ++it)
    {
      result << Tools::color("blue") << it->first << ": " << Tools::color() << std::endl;

      for (auto &key : it->second.getKeys())
      {
        result << it->second.keyToStr(key) << ": " << "(" << it->second(key).n_rows << "x" << it->second(key).n_cols  << "x" << it->second(key).n_slices << ")" << std::endl;
      }
    }
  }

  std::string tempStr = result.str();

  if (!tempStr.empty()) tempStr.pop_back(); // remove the last std::endl;

  return tempStr;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a list of keys containing a pattern.
 */

std::set<std::string> DataTree::getKeys(const std::string &pattern) const
{
  std::set<std::string> result;


  FIND_AND_INSERT(key, intMap       );
  FIND_AND_INSERT(key, doubleMap    );
  FIND_AND_INSERT(key, stringMap    );
  FIND_AND_INSERT(key, vecMap       );
  FIND_AND_INSERT(key, matMap       );
  FIND_AND_INSERT(key, cubeMap      );
  FIND_AND_INSERT(key, ivecMap      );
  FIND_AND_INSERT(key, imatMap      );
  FIND_AND_INSERT(key, icubeMap     );
  FIND_AND_INSERT(key, uvecMap      );
  FIND_AND_INSERT(key, umatMap      );
  FIND_AND_INSERT(key, ucubeMap     );
  FIND_AND_INSERT(key, multiIVecMap );
  FIND_AND_INSERT(key, multiIMatMap );
  FIND_AND_INSERT(key, multiICubeMap);
  FIND_AND_INSERT(key, multiUVecMap );
  FIND_AND_INSERT(key, multiUMatMap );
  FIND_AND_INSERT(key, multiUCubeMap);
  FIND_AND_INSERT(key, multiVecMap  );
  FIND_AND_INSERT(key, multiMatMap  );
  FIND_AND_INSERT(key, multiCubeMap );

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

DataTree DataTree::operator+(const DataTree &other) const
{
  // warning: A = B + C may overwrite B's entries with C's entries
  DataTree result(*this);
  result.merge(other);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Merge the current instance with an other DataTree instance (overwriting identical keys).
 */

void DataTree::merge(const DataTree &other)
{
  for (auto &key : other.boolMap        ) set(key.first, key.second);
  for (auto &key : other.intMap         ) set(key.first, key.second);
  for (auto &key : other.doubleMap      ) set(key.first, key.second);
  for (auto &key : other.stringMap      ) set(key.first, key.second);
  for (auto &key : other.uvecMap        ) set(key.first, key.second);
  for (auto &key : other.umatMap        ) set(key.first, key.second);
  for (auto &key : other.ucubeMap       ) set(key.first, key.second);
  for (auto &key : other.ivecMap        ) set(key.first, key.second);
  for (auto &key : other.imatMap        ) set(key.first, key.second);
  for (auto &key : other.icubeMap       ) set(key.first, key.second);
  for (auto &key : other.vecMap         ) set(key.first, key.second);
  for (auto &key : other.matMap         ) set(key.first, key.second);
  for (auto &key : other.cubeMap        ) set(key.first, key.second);
  for (auto &key : other.multiUVecMap   ) set(key.first, key.second);
  for (auto &key : other.multiUMatMap   ) set(key.first, key.second);
  for (auto &key : other.multiUCubeMap  ) set(key.first, key.second);
  for (auto &key : other.multiIVecMap   ) set(key.first, key.second);
  for (auto &key : other.multiIMatMap   ) set(key.first, key.second);
  for (auto &key : other.multiICubeMap  ) set(key.first, key.second);
  for (auto &key : other.multiVecMap    ) set(key.first, key.second);
  for (auto &key : other.multiMatMap    ) set(key.first, key.second);
  for (auto &key : other.multiCubeMap   ) set(key.first, key.second);

  for (auto &key : other.emptyKeys) del(key);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Validate the keys.
 */

void DataTree::validate(void) const
{
  // Look for duplicates

  std::set<std::string> keys;
  for (auto &key : boolMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : intMap       ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : doubleMap    ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : uvecMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : umatMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : ucubeMap     ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : ivecMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : imatMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : icubeMap     ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : vecMap       ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : matMap       ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : cubeMap      ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiUVecMap ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiUMatMap ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiUCubeMap) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiIVecMap ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiIMatMap ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiICubeMap) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiVecMap  ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiMatMap  ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : multiCubeMap ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };
  for (auto &key : stringMap    ) { if (!(keys.insert(key.first).second)) ERROR("Invalid DataTree: duplicate key found: '" + key.first + "'"); };

  // Ckeck key types
  for (auto &key : boolMap      ) { if (!(isValid(key.first, "B"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : intMap       ) { if (!(isValid(key.first, "I"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : doubleMap    ) { if (!(isValid(key.first, "D"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : uvecMap      ) { if (!(isValid(key.first, "UV" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : umatMap      ) { if (!(isValid(key.first, "UM" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : ucubeMap     ) { if (!(isValid(key.first, "UC" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : ivecMap      ) { if (!(isValid(key.first, "IV" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : imatMap      ) { if (!(isValid(key.first, "IM" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : icubeMap     ) { if (!(isValid(key.first, "IC" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : vecMap       ) { if (!(isValid(key.first, "V"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : matMap       ) { if (!(isValid(key.first, "M"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : cubeMap      ) { if (!(isValid(key.first, "C"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiUVecMap ) { if (!(isValid(key.first, "MUV"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiUMatMap ) { if (!(isValid(key.first, "MUM"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiUCubeMap) { if (!(isValid(key.first, "MUC"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiIVecMap ) { if (!(isValid(key.first, "MIV"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiIMatMap ) { if (!(isValid(key.first, "MIM"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiICubeMap) { if (!(isValid(key.first, "MIC"))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiVecMap  ) { if (!(isValid(key.first, "MV" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiMatMap  ) { if (!(isValid(key.first, "MM" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : multiCubeMap ) { if (!(isValid(key.first, "MC" ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
  for (auto &key : stringMap    ) { if (!(isValid(key.first, "S"  ))) ERROR("Invalid DataTree: wrong type or invalid key: '" + key.first + "'"); };
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Test if a key exists.
 */

bool DataTree::contains(const std::string &key, const std::string &type) const
{

  if ((type == "int"   ) || (type == "")) for (auto &k : intMap        ) if (k.first == key) return true;
  if ((type == "double") || (type == "")) for (auto &k : doubleMap     ) if (k.first == key) return true;
  if ((type == "string") || (type == "")) for (auto &k : stringMap     ) if (k.first == key) return true;
  if ((type == "uvec"  ) || (type == "")) for (auto &k : uvecMap       ) if (k.first == key) return true;
  if ((type == "umat"  ) || (type == "")) for (auto &k : umatMap       ) if (k.first == key) return true;
  if ((type == "ucube" ) || (type == "")) for (auto &k : ucubeMap      ) if (k.first == key) return true;
  if ((type == "ivec"  ) || (type == "")) for (auto &k : ivecMap       ) if (k.first == key) return true;
  if ((type == "imat"  ) || (type == "")) for (auto &k : imatMap       ) if (k.first == key) return true;
  if ((type == "icube" ) || (type == "")) for (auto &k : icubeMap      ) if (k.first == key) return true;
  if ((type == "vec"   ) || (type == "")) for (auto &k : vecMap        ) if (k.first == key) return true;
  if ((type == "mat"   ) || (type == "")) for (auto &k : matMap        ) if (k.first == key) return true;
  if ((type == "cube"  ) || (type == "")) for (auto &k : cubeMap       ) if (k.first == key) return true;
  if ((type == "Mumat" ) || (type == "")) for (auto &k : multiUVecMap  ) if (k.first == key) return true;
  if ((type == "Mimat" ) || (type == "")) for (auto &k : multiUMatMap  ) if (k.first == key) return true;
  if ((type == "Micube") || (type == "")) for (auto &k : multiUCubeMap ) if (k.first == key) return true;
  if ((type == "Mivec" ) || (type == "")) for (auto &k : multiIVecMap  ) if (k.first == key) return true;
  if ((type == "Mimat" ) || (type == "")) for (auto &k : multiIMatMap  ) if (k.first == key) return true;
  if ((type == "Micube") || (type == "")) for (auto &k : multiICubeMap ) if (k.first == key) return true;
  if ((type == "Mvec"  ) || (type == "")) for (auto &k : multiVecMap   ) if (k.first == key) return true;
  if ((type == "Mmat"  ) || (type == "")) for (auto &k : multiMatMap   ) if (k.first == key) return true;
  if ((type == "Mcube" ) || (type == "")) for (auto &k : multiCubeMap  ) if (k.first == key) return true;

  return false;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Test if a key containing a pattern exists.
 */

bool DataTree::containsPattern(const std::string &pattern) const
{
  for (auto &key : boolMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : intMap       ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : doubleMap    ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : stringMap    ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : uvecMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : umatMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : ucubeMap     ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : ivecMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : imatMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : icubeMap     ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : vecMap       ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : matMap       ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : cubeMap      ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiUVecMap ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiUMatMap ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiUCubeMap) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiIVecMap ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiIMatMap ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiICubeMap) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiVecMap  ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiMatMap  ) if (key.first.find(pattern) != std::string::npos) return true;
  for (auto &key : multiCubeMap ) if (key.first.find(pattern) != std::string::npos) return true;

  return false;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the type of a key.
 */

const std::string DataTree::getType(const std::string &key) const
{
  for (auto &o : general.globalValidKeys)
  {
    if (o.key != key ) continue;

    return o.type;
  }

  return "";
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Test if a couple (key, type) is valid.
 */

bool DataTree::isValid(const std::string &key, const std::string &type) const
{
  for (auto &o : general.globalValidKeys)
  {
    if (o.key != key ) continue;
    if ((type != "") && (o.type != type)) continue;

    return true;
  }

  DEBUG("%ld keys", general.globalValidKeys.size());
  INFO("invalid key: " + key + " type: " + type);

  return false;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get a DataTree instance filled with default values for valid keys.
  */

DataTree DataTree::getDefault(void)
{
  DBG_ENTER;

  DataTree result;
  result.setDefault();

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Fill with default values.
  */

void DataTree::setDefault(void)
{
  DBG_ENTER;

  for (auto &o : general.globalValidKeys)
  {
    if (!o.defaultValue.empty())
      IOhfb3::updateDataTree(*this, o.key, o.defaultValue, o.type);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check equality between DataTree instances.
 */

bool DataTree::operator==(const DataTree &other) const
{
  if (intMap != other.intMap) return false;
  if (doubleMap != other.doubleMap) return false;
  if (stringMap != other.stringMap) return false;
  if (!mapStringArmaEquals(uvecMap, other.uvecMap )) return false;
  if (!mapStringArmaEquals(umatMap, other.umatMap )) return false;
  if (!mapStringArmaEquals(ucubeMap, other.ucubeMap)) return false;
  if (!mapStringArmaEquals(ivecMap, other.ivecMap )) return false;
  if (!mapStringArmaEquals(imatMap, other.imatMap )) return false;
  if (!mapStringArmaEquals(icubeMap, other.icubeMap)) return false;
  if (!mapStringArmaEquals(vecMap, other.vecMap  )) return false;
  if (!mapStringArmaEquals(matMap, other.matMap  )) return false;
  if (!mapStringArmaEquals(cubeMap, other.cubeMap )) return false;
  if (!mapStringArmaEquals(multiIVecMap, other.multiIVecMap  )) return false;
  if (!mapStringArmaEquals(multiIMatMap, other.multiIMatMap  )) return false;
  if (!mapStringArmaEquals(multiICubeMap, other.multiICubeMap )) return false;
  if (!mapStringArmaEquals(multiUVecMap, other.multiUVecMap  )) return false;
  if (!mapStringArmaEquals(multiUMatMap, other.multiUMatMap  )) return false;
  if (!mapStringArmaEquals(multiUCubeMap, other.multiUCubeMap )) return false;
  if (!mapStringArmaEquals(multiVecMap, other.multiVecMap  )) return false;
  if (!mapStringArmaEquals(multiMatMap, other.multiMatMap  )) return false;
  if (!mapStringArmaEquals(multiCubeMap, other.multiCubeMap )) return false;
  if (emptyKeys != other.emptyKeys) return false;

  return true;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, std::string) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::setS(const std::string &key, const std::string &val)
{
  stringMap[key] = std::regex_replace(val, std::regex("\n"), " ");
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, std::string) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::setSC(const std::string &key, const char *val)
{
  set(key, std::string(val));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, bool) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::setB(const std::string &key, const bool &val)
{
  set(key, val);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or set an empty key.
 *
 *  \param key The key.
 */

void DataTree::setE(const std::string &key)
{
  del(key);
  emptyKeys.insert(key);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, INT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::setI(const std::string &key, INT val)
{
  if (doubleMap.find(key) != doubleMap.end()) // set as double if doubleMap key already exists
    doubleMap[key] = double(val);
  else
    intMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, double) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The value.
 */

void DataTree::setD(const std::string &key, double val)
{
  intMap.erase(key);
  doubleMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::vec) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::vec object.
 */

void DataTree::setV(const std::string &key, const arma::vec &val)
{
  vecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::mat) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::mat object.
 */

void DataTree::setM(const std::string &key, const arma::mat &val)
{
  matMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, arma::cube) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The arma::cube object.
 */

void DataTree::setC(const std::string &key, const arma::cube &val)
{
  cubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, IVEC) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The IVEC object.
 */

void DataTree::setIV(const std::string &key, const IVEC &val)
{
  ivecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, IMAT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The IMAT object.
 */

void DataTree::setIM(const std::string &key, const IMAT &val)
{
  imatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, ICUBE) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The ICUBE object.
 */

void DataTree::setIC(const std::string &key, const ICUBE &val)
{
  icubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<IVEC >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<IVEC > object.
 */

void DataTree::setMIV(const std::string &key, const Multi<IVEC> &val)
{
  multiIVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<IMAT >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<IMAT > object.
 */

void DataTree::setMIM(const std::string &key, const Multi<IMAT> &val)
{
  multiIMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<ICUBE>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<ICUBE> object.
 */

void DataTree::setMIC(const std::string &key, const Multi<ICUBE> &val)
{
  multiICubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UVEC) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UVEC object.
 */

void DataTree::setUV(const std::string &key, const UVEC &val)
{
  uvecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UMAT) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UMAT object.
 */

void DataTree::setUM(const std::string &key, const UMAT &val)
{
  umatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, UCUBE) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The UCUBE object.
 */

void DataTree::setUC(const std::string &key, const UCUBE &val)
{
  ucubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UVEC >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UVEC > object.
 */

void DataTree::setMUV(const std::string &key, const Multi<UVEC> &val)
{
  multiUVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UMAT >) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UMAT > object.
 */

void DataTree::setMUM(const std::string &key, const Multi<UMAT> &val)
{
  multiUMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<UCUBE>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<UCUBE> object.
 */

void DataTree::setMUC(const std::string &key, const Multi<UCUBE> &val)
{
  multiUCubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::vec>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::vec> object.
 */

void DataTree::setMV(const std::string &key, const Multi<arma::vec> &val)
{
  multiVecMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::mat>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::m  t> object.
 */

void DataTree::setMM(const std::string &key, const Multi<arma::mat> &val)
{
  multiMatMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create or update a (key, Multi<arma::cube>) couple in the data tree.
 *
 *  \param key The key.
 *  \param val The Multi<arma::cube> object.
 */

void DataTree::setMC(const std::string &key, const Multi<arma::cube> &val)
{
  multiCubeMap[key] = val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree from a section of an existing DataTree.
 *
 *  \param sectionName The section name.
 */

DataTree DataTree::getSection(const std::string &sectionName) const
{
  DBG_ENTER;

  DataTree result;

  if (intMap.size() != 0)
  {
    for (std::map<std::string, INT>::const_iterator it = intMap.begin(); it != intMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (doubleMap.size() != 0)
  {
    for (std::map<std::string, double>::const_iterator it = doubleMap.begin(); it != doubleMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (stringMap.size() != 0)
  {
    for (std::map<std::string, std::string>::const_iterator it = stringMap.begin(); it != stringMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (vecMap.size() != 0)
  {
    for (std::map<std::string, arma::vec>::const_iterator it = vecMap.begin(); it != vecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (matMap.size() != 0)
  {
    for (std::map<std::string, arma::mat>::const_iterator it = matMap.begin(); it != matMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (cubeMap.size() != 0)
  {
    for (std::map<std::string, arma::cube>::const_iterator it = cubeMap.begin(); it != cubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (uvecMap.size() != 0)
  {
    for (std::map<std::string, UVEC >::const_iterator it = uvecMap.begin(); it != uvecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (umatMap.size() != 0)
  {
    for (std::map<std::string, UMAT >::const_iterator it = umatMap.begin(); it != umatMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (ucubeMap.size() != 0)
  {
    for (std::map<std::string, UCUBE >::const_iterator it = ucubeMap.begin(); it != ucubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiUVecMap.size() != 0)
  {
    for (std::map<std::string, Multi<UVEC > >::const_iterator it = multiUVecMap.begin(); it != multiUVecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiUMatMap.size() != 0)
  {
    for (std::map<std::string, Multi<UMAT > >::const_iterator it = multiUMatMap.begin(); it != multiUMatMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiUCubeMap.size() != 0)
  {
    for (std::map<std::string, Multi<UCUBE > >::const_iterator it = multiUCubeMap.begin(); it != multiUCubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (ivecMap.size() != 0)
  {
    for (std::map<std::string, IVEC >::const_iterator it = ivecMap.begin(); it != ivecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (imatMap.size() != 0)
  {
    for (std::map<std::string, IMAT >::const_iterator it = imatMap.begin(); it != imatMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (icubeMap.size() != 0)
  {
    for (std::map<std::string, ICUBE >::const_iterator it = icubeMap.begin(); it != icubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiIVecMap.size() != 0)
  {
    for (std::map<std::string, Multi<IVEC > >::const_iterator it = multiIVecMap.begin(); it != multiIVecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiIMatMap.size() != 0)
  {
    for (std::map<std::string, Multi<IMAT > >::const_iterator it = multiIMatMap.begin(); it != multiIMatMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiICubeMap.size() != 0)
  {
    for (std::map<std::string, Multi<ICUBE > >::const_iterator it = multiICubeMap.begin(); it != multiICubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiVecMap.size() != 0)
  {
    for (std::map<std::string, Multi<arma::vec> >::const_iterator it = multiVecMap.begin(); it != multiVecMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiMatMap.size() != 0)
  {
    for (std::map<std::string, Multi<arma::mat> >::const_iterator it = multiMatMap.begin(); it != multiMatMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  if (multiCubeMap.size() != 0)
  {
    for (std::map<std::string, Multi<arma::cube> >::const_iterator it = multiCubeMap.begin(); it != multiCubeMap.end(); ++it)
      if (it->first.rfind(sectionName, 0) == 0) result.set(it->first.substr(sectionName.size(), std::string::npos), it->second);
  }

  DBG_RETURN(result);
}



