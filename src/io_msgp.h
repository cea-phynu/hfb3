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

#ifndef IO_MSGP_H
#define IO_MSGP_H

/** \file
 *  \brief Headers for the IOmsgp class.
 */

#include "global.h"
#include "multi.h"

class DataTree;

/** \brief Interface to Msgpack files.
 *
 * This class provides methods to read / write .msg files (using msgpack).
 */

class IOmsgp
{
public:

  IOmsgp(void);                                                        // #TEST#
  static DataTree fromContent(const std::string &content);             // #TEST#
  static void saveDataTree(const DataTree &,
                           const std::string &filename);
  static std::string serializeDataTree(const DataTree &);

  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  static arma::vec  array2vec( const std::string &v);
  static arma::mat  array2mat( const std::string &v);
  static arma::cube array2cube(const std::string &v);

  static std::tuple<std::string, std::string> vec2array( const arma::vec  &v);
  static std::tuple<std::string, std::string> mat2array( const arma::mat  &v);
  static std::tuple<std::string, std::string> cube2array(const arma::cube &v);

  static IVEC  array2ivec( const std::string &v);
  static IMAT  array2imat( const std::string &v);
  static ICUBE array2icube(const std::string &v);

  static UVEC  array2uvec( const std::string &v);
  static UMAT  array2umat( const std::string &v);
  static UCUBE array2ucube(const std::string &v);

  static std::tuple<std::string, std::string> ivec2array( const IVEC  &v);
  static std::tuple<std::string, std::string> imat2array( const IMAT  &v);
  static std::tuple<std::string, std::string> icube2array(const ICUBE &v);

  static std::tuple<std::string, std::string> uvec2array( const UVEC  &v);
  static std::tuple<std::string, std::string> umat2array( const UMAT  &v);
  static std::tuple<std::string, std::string> ucube2array(const UCUBE &v);

  static Multi<IVEC > array2multiIVec( const std::map<std::vector<INT>, std::string> &v);
  static Multi<IMAT > array2multiIMat( const std::map<std::vector<INT>, std::string> &v);
  static Multi<ICUBE> array2multiICube(const std::map<std::vector<INT>, std::string> &v);

  static Multi<UVEC > array2multiUVec( const std::map<std::vector<INT>, std::string> &v);
  static Multi<UMAT > array2multiUMat( const std::map<std::vector<INT>, std::string> &v);
  static Multi<UCUBE> array2multiUCube(const std::map<std::vector<INT>, std::string> &v);

  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiIVec2array( const Multi<IVEC > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiIMat2array( const Multi<IMAT > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiICube2array(const Multi<ICUBE> &v);

  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiUVec2array( const Multi<UVEC > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiUMat2array( const Multi<UMAT > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiUCube2array(const Multi<UCUBE> &v);

  static Multi<arma::vec > array2multiVec( const std::map<std::vector<INT>, std::string> &v);
  static Multi<arma::mat > array2multiMat( const std::map<std::vector<INT>, std::string> &v);
  static Multi<arma::cube> array2multiCube(const std::map<std::vector<INT>, std::string> &v);

  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiVec2array( const Multi<arma::vec > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiMat2array( const Multi<arma::mat > &v);
  static std::tuple<std::string, std::map<std::vector<INT>, std::string> > multiCube2array(const Multi<arma::cube> &v);
};

#endif // IO_MSGP_H
