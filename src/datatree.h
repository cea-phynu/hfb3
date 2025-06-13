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

#ifndef DATATREE_H
#define DATATREE_H

/** \file
 *  \brief Headers for the DataTree class.
 */

#include "global.h"
#include "generic.h"
#include "multi.h"

/** \brief A data tree used as an interface between the solver and the external world.
 *
 *  This class represents the main way for the solver to load and save data.
 *
 *  Example:
 *
 * \code{.cpp}
    DataTree d;

    d.set("songs/first", std::string("La fille du Bédouin..."));
    d.set("songs/second", std::string("C'est un fameux trois-mâts..."));
    d.save("../data/save.msg.gz");

    d.clear();

    d.load(Tools::readFile("../data/save.msg.gz"));
    std::cout << d.getS("songs/second") << std::endl; // Print "C'est un fameux trois-mâts...".
    std::cout << d.getS("songs/first") << std::endl; // Print "La fille du Bédouin...".
    \endcode
 *
 */

class DataTree : public Generic
{
public:


  // TODO: add tests !

  DataTree(const std::string &filename = "");                          // #TEST#

  void set(const std::string &key, const std::string &val);
  void set(const std::string &key, const char *val);
  void set(const std::string &key, INT val);
  void set(const std::string &key, double val);
  void set(const std::string &key, const arma::vec &val);
  void set(const std::string &key, const arma::mat &val);
  void set(const std::string &key, const arma::cube &val);
  void set(const std::string &key, const UVEC &val);
  void set(const std::string &key, const UMAT &val);
  void set(const std::string &key, const UCUBE &val);
  void set(const std::string &key, const IVEC &val);
  void set(const std::string &key, const IMAT &val);
  void set(const std::string &key, const ICUBE &val);
  void set(const std::string &key, const Multi<UVEC > &val);
  void set(const std::string &key, const Multi<UMAT > &val);
  void set(const std::string &key, const Multi<UCUBE > &val);
  void set(const std::string &key, const Multi<IVEC > &val);
  void set(const std::string &key, const Multi<IMAT > &val);
  void set(const std::string &key, const Multi<ICUBE > &val);
  void set(const std::string &key, const Multi<arma::vec> &val);
  void set(const std::string &key, const Multi<arma::mat> &val);
  void set(const std::string &key, const Multi<arma::cube> &val);
  void set(const std::string &key, const bool &val);
  void set(const std::string &key);

  bool isEmpty(void) const;
  UINT size(void) const;

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  std::string abstract(void) const;

  std::set<std::string> getKeys(const std::string &pattern) const;
  bool contains(const std::string &key, const std::string &type = "") const;
  bool containsPattern(const std::string &pattern) const;

  bool get(INT &, const std::string &key, bool ignore = false) const;
  bool get(bool &, const std::string &key, bool ignore = false) const;
  bool get(double &, const std::string &key, bool ignore = false) const;
  bool get(std::string &, const std::string &key, bool ignore = false) const;

  bool get(arma::vec &, const std::string &key, bool ignore = false) const;
  bool get(arma::mat &, const std::string &key, bool ignore = false) const;

  bool get(arma::cube &, const std::string &key, bool ignore = false) const;
  bool get(IVEC &, const std::string &key, bool ignore = false) const;
  bool get(IMAT &, const std::string &key, bool ignore = false) const;
  bool get(ICUBE &, const std::string &key, bool ignore = false) const;
  bool get(UVEC &, const std::string &key, bool ignore = false) const;
  bool get(UMAT &, const std::string &key, bool ignore = false) const;
  bool get(UCUBE &, const std::string &key, bool ignore = false) const;
  bool get(Multi<IVEC > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<IMAT > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<ICUBE > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<UVEC > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<UMAT > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<UCUBE > &, const std::string &key, bool ignore = false) const;
  bool get(Multi<arma::vec> &, const std::string &key, bool ignore = false) const;
  bool get(Multi<arma::mat> &, const std::string &key, bool ignore = false) const;
  bool get(Multi<arma::cube> &, const std::string &key, bool ignore = false) const;

  DataTree getSection(const std::string &sectionName) const;

  bool         getB( const std::string &key, bool ignore = false, bool *isFound = NULL) const;

  // For python
  INT          getI( const std::string &key, bool ignore = false) const;
  double       getD( const std::string &key, bool ignore = false) const;
  std::string  getS( const std::string &key, bool ignore = false) const;
  arma::vec    getV( const std::string &key, bool ignore = false) const;
  arma::mat    getM( const std::string &key, bool ignore = false) const;
  arma::cube   getC( const std::string &key, bool ignore = false) const;
  IVEC   getIV(const std::string &key, bool ignore = false) const;
  IMAT   getIM(const std::string &key, bool ignore = false) const;
  ICUBE  getIC(const std::string &key, bool ignore = false) const;
  UVEC   getUV(const std::string &key, bool ignore = false) const;
  UMAT   getUM(const std::string &key, bool ignore = false) const;
  UCUBE  getUC(const std::string &key, bool ignore = false) const;
  Multi<IVEC >   getMIV(const std::string &key, bool ignore = false) const;
  Multi<IMAT >   getMIM(const std::string &key, bool ignore = false) const;
  Multi<ICUBE >  getMIC(const std::string &key, bool ignore = false) const;
  Multi<UVEC >   getMUV(const std::string &key, bool ignore = false) const;
  Multi<UMAT >   getMUM(const std::string &key, bool ignore = false) const;
  Multi<UCUBE >  getMUC(const std::string &key, bool ignore = false) const;
  Multi<arma::vec>   getMV(const std::string &key, bool ignore = false) const;
  Multi<arma::mat>   getMM(const std::string &key, bool ignore = false) const;
  Multi<arma::cube>  getMC(const std::string &key, bool ignore = false) const;

  void setS(const std::string &key, const std::string &val);
  void setSC(const std::string &key, const char *val);
  void setI(const std::string &key, INT val);
  void setD(const std::string &key, double val);
  void setV(const std::string &key, const arma::vec &val);
  void setM(const std::string &key, const arma::mat &val);
  void setC(const std::string &key, const arma::cube &val);
  void setUV(const std::string &key, const UVEC &val);
  void setUM(const std::string &key, const UMAT &val);
  void setUC(const std::string &key, const UCUBE &val);
  void setIV(const std::string &key, const IVEC &val);
  void setIM(const std::string &key, const IMAT &val);
  void setIC(const std::string &key, const ICUBE &val);
  void setMUV(const std::string &key, const Multi<UVEC > &val);
  void setMUM(const std::string &key, const Multi<UMAT > &val);
  void setMUC(const std::string &key, const Multi<UCUBE > &val);
  void setMIV(const std::string &key, const Multi<IVEC > &val);
  void setMIM(const std::string &key, const Multi<IMAT > &val);
  void setMIC(const std::string &key, const Multi<ICUBE > &val);
  void setMV(const std::string &key, const Multi<arma::vec> &val);
  void setMM(const std::string &key, const Multi<arma::mat> &val);
  void setMC(const std::string &key, const Multi<arma::cube> &val);
  void setB(const std::string &key, const bool &val);
  void setE(const std::string &key);
  void setDefault(void);
  static DataTree getDefault(void);

  std::list<KeyStruct > listOfKeys;

  void clear(void);
  void clean(void);
  void merge(const DataTree &other);
  DataTree operator+(const DataTree &other) const;

  INT del(const std::string &key);

  void validate(void) const;
  bool isValid(const std::string &key, const std::string &type) const;

  void save(const std::string &filename, bool verbose = true) const;

  static DataTree fromContent(const std::string &content);

  bool operator==(const DataTree &) const;
  bool operator!=(const DataTree &other) const
  {
    return !this->operator==(other);
  };

  //============================================================================
  //============================================================================
  //============================================================================

  /// Equality operator between std::map<std::string, arma::xxx> objects
  template<typename T> static bool mapStringArmaEquals(const std::map<std::string, T> &a, const std::map<std::string, T> &b)
  {
    std::vector<std::string> keysa;
    std::vector<std::string> keysb;
    std::transform(std::begin(a), std::end(a), std::back_inserter(keysa), [](std::pair<std::string, T> const & pair)
    {
      return pair.first;
    });
    std::transform(std::begin(b), std::end(b), std::back_inserter(keysb), [](std::pair<std::string, T> const & pair)
    {
      return pair.first;
    });

    std::sort(keysa.begin(), keysa.end());
    std::sort(keysb.begin(), keysb.end());
    if (keysa != keysb) return false;

    for (auto &key : keysa)
    {
      if (!(a.at(key) == b.at(key))) return false;
    }

    return true;
  }

  //============================================================================
  //============================================================================
  //============================================================================

  /// Return a std::string representation of an std::map<...>.
  template<typename V> static const std::string mapToStr(const std::string &title, const std::map<std::string, V> &map)
  {
    std::string result = "";

    if (map.size() == 0) return "none";

    std::vector<std::pair<std::string, std::string> > list;

    list.push_back({title, ""});

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      list.push_back({it->first, Tools::infoStr(it->second)});
    }

    result += Tools::treeStr(list, false);

    /*
    for (auto it = map.begin(); it != map.end(); ++it)
    {
      if (it != map.begin()) result += ", ";
      result += "\"" + it->first + "\"";
    }
    */

    return result;
  }

  /// Return a std::string representation of an std::map<...>.
  template<typename V> static const std::string multiMapToStr(const std::string &title, const std::map<std::string, V> &map)
  {
    std::string result = "";

    if (map.size() == 0) return "none";

    std::vector<std::pair<std::string, std::string> > list;

    list.push_back({title, ""});

    for (auto it = map.begin(); it != map.end(); ++it)
    {
      list.push_back({it->first, it->second.info(true)});
    }

    result += Tools::treeStr(list, false);

    return result;
  }

  //============================================================================
  //============================================================================
  //============================================================================

  /// Data tree for ints.
  std::map<std::string, INT> intMap;

  /// Data tree for doubles.
  std::map<std::string, double> doubleMap;

  /// Data tree for std::strings.
  std::map<std::string, std::string> stringMap;

  /// Data tree for IVEC.
  std::map<std::string, IVEC > ivecMap;

  /// Data tree for UVEC.
  std::map<std::string, UVEC > uvecMap;

  /// Data tree for IMAT.
  std::map<std::string, IMAT > imatMap;

  /// Data tree for UMAT.
  std::map<std::string, UMAT > umatMap;

  /// Data tree for ICUBE.
  std::map<std::string, ICUBE > icubeMap;

  /// Data tree for UCUBE.
  std::map<std::string, UCUBE > ucubeMap;

  /// Data tree for arma::vec.
  std::map<std::string, arma::vec> vecMap;

  /// Data tree for arma::mat.
  std::map<std::string, arma::mat > matMap;

  /// Data tree for arma::cube.
  std::map<std::string, arma::cube > cubeMap;

  /// Data tree for Multi<IVEC >.
  std::map<std::string, Multi<IVEC > > multiIVecMap;

  /// Data tree for Multi<UVEC >.
  std::map<std::string, Multi<UVEC > > multiUVecMap;

  /// Data tree for Multi<IMAT >.
  std::map<std::string, Multi<IMAT > > multiIMatMap;

  /// Data tree for Multi<UMAT >.
  std::map<std::string, Multi<UMAT > > multiUMatMap;

  /// Data tree for Multi<ICUBE>.
  std::map<std::string, Multi<ICUBE > > multiICubeMap;

  /// Data tree for Multi<UCUBE>.
  std::map<std::string, Multi<UCUBE > > multiUCubeMap;

  /// Data tree for Multi<arma::vec>.
  std::map<std::string, Multi<arma::vec> > multiVecMap;

  /// Data tree for Multi<arma::mat>.
  std::map<std::string, Multi<arma::mat> > multiMatMap;

  /// Data tree for Multi<arma::cube>.
  std::map<std::string, Multi<arma::cube> > multiCubeMap;

  /// Empty keys : keys that are automatically suppress if a merge appened.
  std::set<std::string> emptyKeys;
};


#endif // DATATREE_H
