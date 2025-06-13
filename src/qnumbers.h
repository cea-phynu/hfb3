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

#ifndef QNUMBERS_H
#define QNUMBERS_H

/** \file
 *  \brief Headers for the Qnumbers class.
 */

#include "global.h"
#include "generic.h"
#include "tools.h"
#include <unordered_map>

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Function calculating the hash of an IVEC.
 *
 * This function calculates the hash of an IVEC. It allows IVEC to be used as a key in an std::unordered_map<>.
 */

struct VecHash
{
  std::size_t operator()(const IVEC &k) const
  {
    return Tools::ivec2hash(k);
  }
};

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Function testing the equality of two IVEC instances.
 *
 * This function tests the equality of two IVEC instances. It allows IVEC to be used as a key in an std::unordered_map<>.
 */

struct VecEqual
{
  bool operator()(const IVEC &a, const IVEC &b) const
  {
    return (a == b);
  }
};

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief The methods to store and manipulate quantum numbers sets.
 *
 *  This class implements methods to store and manipulate quantum numbers sets.
 */


class Qnumbers : public Generic
{
public:

  Qnumbers(UINT _n = 0);                                               // #TEST#
  Qnumbers(std::string, UINT);                                         // #TEST#
  void setNames(const std::vector<std::string> &l);                    // #TEST#
  INT find(const IVEC &) const;                                        // #TEST#
  void appendName(const std::string &s);                               // #TEST#
  UINT append(const IVEC &);                                           // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  void reorder(IMAT r);                                                // #TEST#
  INT operator()(UINT, UINT) const;                                    // #TEST#
  const Qnumbers extract(UINT) const;                                  // #TEST#
  IVEC operator()(UINT) const;                                         // #TEST#
  UINT calcBlocks(const UVEC &rows);                                   // #TEST#
  const Qnumbers sub(UVEC ids) const;                                  // #TEST#
  // Symmetry check
  void checkOmegaSym(const arma::mat &mat, const std::string &name);   // #TEST#
  void clear(void);                                                    // #TEST#
  const std::string niceStr(UINT id) const;                            // #TEST#

  /// The number of quantum numbers.
  UINT nbq;

  /// The number of values of quantum numbers.
  UINT nb;

  /// The quantum number names.
  std::vector<std::string> name;

  /// The block structure (if defined).
  std::vector<Qnumbers> blocks;

  /// The parent block filter (if defined).
  UVEC filter;

  /// The quantum numbers.
  IMAT mat;

private:
  /// The index map structure.
  std::unordered_map<IVEC, ulong, VecHash, VecEqual> index;


};

#endif
