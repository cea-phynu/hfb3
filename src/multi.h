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

#ifndef MULTI_H
#define MULTI_H

/** \file
 *  \brief Headers for the Multi class.
 */

#include "global.h"
#include <unordered_map>
#include <array>
#include "typo.h"
#include "tools.h"

//==============================================================================

/** \brief Specialized hasher for the std::vector<INT> keys.
 */

struct myHasher
{
  std::size_t operator()(const std::vector<INT> &a) const
  {
    std::size_t h = 0;

    for (auto &e : a)
    {
      h ^= std::hash<INT> {}(e) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }

    return h;
  }
};

//==============================================================================

/// The derived type for the Multi class.
#define BASE_TYPE std::unordered_map<std::vector<INT>, T, myHasher >
// #define BASE_TYPE tsl::robin_map<std::vector<INT>, T, myHasher >
// #define BASE_TYPE std::map<std::array<INT, N>, T>
// #define BASE_TYPE std::unordered_map<std::array<INT, N>, T, myHasher<N> >

/** \brief Generic x-D container.
 *
 *  This class implements methods to store objects using an std::unordered_map<std::vector<INT>, ...> container.
 */

template<class T>
class Multi : public BASE_TYPE
{
public:

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Access operator.
  template<class...T2>
  T & operator()(T2... args)
  {
    return this->operator[]({args...});
  }


  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Access operator (const).
  template<class...T2>
  const T & operator()(T2... args) const                                 // #TEST#
  {
    // if (this->find( {args...} ) == this->end())
    //   ERROR("key does not exist: " + this->keyToStr({args...}));

    return this->at({args...});
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Equality operators for Multi objects.
  bool operator==(const Multi<T> &other) const
  {
    std::vector<std::vector<INT> > keysa = this->getKeys();
    std::vector<std::vector<INT> > keysb = other.getKeys();

    if (keysa != keysb) return false;

    for (auto &key : keysa)
    {
      if(this->operator()(key) != other(key)) return false;
    }

    return true;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Test for emptiness.
  bool empty(void) const
  {
    return BASE_TYPE::empty();
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  // /// Remove a key.
  // template<class...T2>
  // void remove(T2... args)
  // {
  //   auto iter = this->find({args...});
  //
  //   if(iter != this->end())
  //     this->erase(iter);
  // }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Clear the object.
  void clear(void)
  {
    BASE_TYPE::clear();
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  std::vector<std::vector<INT> > getKeys(void) const
  {
    std::vector<std::vector<INT> > key;

    for (auto &it : *this)
    {
      key.push_back(it.first);
    }

    std::sort(key.begin(), key.end());

    return key;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  std::vector<T> getValues(void) const
  {
    std::vector<T> vals;

    for (auto &it : *this)
    {
      vals.push_back(it.second);
    }

    return vals;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Nice formatting of a key
  const std::string keyToStr(const std::vector<INT> &k) const
  {
    std::stringstream result;

    result << "[";

    INT N = (INT)(k.size());
    std::vector<INT>::const_iterator it = k.begin();
    for (INT i = 0; i < N; i++)
    {
      result << *it;

      if (i != N - 1) result << ",";

      ++it;
    }

    result << "]";

    return result.str();
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /** Return the number of elements.
   */

  UINT size(void) const
  {
    return getValues().size();
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /** Return some info.
   */

  const std::string info(bool isShort = false) const
  {
    std::string result = "";

    if (isShort)
    {
      result += Tools::treeStr(
      {
        {"type ", typo::type<T>},
        {"size ", Tools::infoStr(size())},
      }, true);
    }
    else
    {
      std::vector<std::pair<std::string, std::string> > content;
      content.push_back(std::make_pair("Content", ""));
      // INT i = 0;
      for (auto &key : this->getKeys())
      {
        content.push_back(std::make_pair(keyToStr(key), Tools::infoStr(this->operator()(key), isShort)));
        // i++;
      }

      result += Tools::treeStr(
      {
        {"Multi", ""},
        {"type  ", typo::type<T>},
        {"size  ", Tools::infoStr(size())},
        {"", Tools::treeStr(content)},
      }, false);
    }

    return result;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /** Check if a given element exists.
   */

  template<class...T2>
  bool contains(T2... args) const
  {
    return (this->find( {args...} ) != this->end());
  }
};

//================================================================================
//================================================================================
//================================================================================

/** \brief Type of the main container.
 */

#define MEM_TYPE arma::field<T>
// #define MEM_TYPE std::vector<T>

/** \brief Specific Fast x-D container.
 *
 *  This class implements a specific fast implementation of a generic x-D container.
 *  The key indices have to take contiguous values and span a regular hypercube.
 */

template<class T>
class FMulti
{
public:

  /// Dummy constructor.
  FMulti() = default;

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Constructor from list of spans (recommended)
  FMulti(std::vector<std::pair<INT, INT>> l)
  {
    IVEC _min(l.size());
    IVEC _max(l.size());

    for (int i=0; i < l.size(); i++) {
      _min(i) = l[i].first;
      _max(i) = l[i].second;
    }
    init(_min, _max);
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Access operator.
  template<typename...T2>
  T & operator()(T2... args)
  {
    UINT i = 0;
    UINT idx = 0;
    UINT mult = 1;
    // INFO("ASK");

#ifdef CHECKBOUNDS_FMULTI
    ASSERT(min.n_elem == sizeof...(args), "wrong number of indices in FMulti: %d indices used, should be %d", min.n_elem, sizeof...(args));
#endif

    for (auto arg : {args...})
    {

#ifdef CHECKBOUNDS_FMULTI
      ASSERT(arg >= min(i) && arg < max(i), PF("index %d error in FMulti: %d not in [%d, %d]", i, arg, min(i), max(i) - 1));
#endif

      idx += (arg-min(i)) * mult;
      mult *= span(i);
      i++;
    }
    // INFO("IDX: %d", idx);
    return mem[idx];
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Access operator (const).
  template<typename...T2>
  const T & operator()(T2... args) const                                 // #TEST#
  {
    UINT i = 0;
    UINT idx = 0;
    UINT mult = 1;
    for (auto arg : {args...})
    {

#ifdef CHECKBOUNDS_FMULTI
      ASSERT(arg >= min(i) && arg < max(i), PF("index %d error in FMulti: %d not in [%d, %d]", i, arg, min(i), max(i) - 1));
#endif

      idx += (arg-min(i)) * mult;
      mult *= span(i);
      i++;
    }

    return mem.at(idx);
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Test for emptiness.
  bool empty(void) const
  {
    //HACK: this is bad, but fast
    if (min.empty()) return true;

    return mem[0].empty();
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  UINT n_elem(void) const
  {
    UINT mult = 1;
    for (auto s: span) {
      mult *= UINT(s);
    }
    return mult;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// Sum size of all contained elements
  UINT size(void) const
  {
    UINT n = n_elem();
    UINT s = 0;
    for (UINT i = 0; i < n; i++) {
      s += mem[i].n_elem;
    }

    return s;
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

private:

  /// Initializer.
  void init(IVEC _min, IVEC _max)
  {
    ASSERT(_min.n_elem == _max.n_elem, "wrong sizes of min and max arguments in FMulti");

    min = _min;
    max = _max;
    span = max - min;

    UINT mult = 1;
    for (auto s: span) {
      mult *= UINT(s);
    }
    // INFO("ALLOC: %d", mult);
    mem = MEM_TYPE(mult);
  }

  //==============================================================================
  //==============================================================================
  //==============================================================================

  /// The container.
  MEM_TYPE mem;

  /// The minimal values of the keys.
  IVEC min;

  /// The maximal values of the keys.
  IVEC max;

  /// The span of the keys.
  IVEC span;
};

#endif // MULTI_H
