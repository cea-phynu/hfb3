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

#ifndef IO_AMEDEE_H
#define IO_AMEDEE_H

#include "global.h"
#include "basis.h"

/** \file
 *  \brief Headers for the IOamedee class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Interface to AMEDEE_DPFP9 results files.
 *
 * This class provides methods to load fort.88 files (results from the AMEDEE_DPFP9 code).
 */

class IOamedee
{
public :
  IOamedee(void);                                                      // #TEST#
  DataTree fromContent(const std::string &content);                    // #TEST#
  static bool checkFileType(const std::string &content);               // #TEST#
  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================

  /// The "customField" matrix [neutrons].
  arma::mat customFieldn;

  /// The "customField" matrix [protons].
  arma::mat customFieldp;

  //============================================================================
  //============================================================================
  //============================================================================

private :
  void readBasisRho(const std::string &content);

  //============================================================================
  //============================================================================
  //============================================================================

  /// Basis object built from the .theo file.
  Basis basis;

  /// The \f$\rho_p\f$ matrices.
  Multi<arma::mat> rho;

  /// The \f$\kappa_p\f$ matrices.
  Multi<arma::mat> kappa;


  /// The number of neutrons.
  INT N;

  /// The number of protons.
  INT Z;

  /// The energy contributions [neutrons].
  std::map<std::string, double> energiesNeut;

  /// The energy contributions [protons].
  std::map<std::string, double> energiesProt;
};

#endif // IOAMEDEE_H
