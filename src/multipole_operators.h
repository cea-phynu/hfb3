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

#ifndef MULTIPOLE_OPERATORS_H
#define MULTIPOLE_OPERATORS_H

/** \file
 *  \brief Headers for the MultipoleOperators class.
 */

#include "global.h"
#include "generic.h"
#include "multi.h"
#include "state.h"

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define MULTIPOLE_OPERATORS_VALID_KEYS \
{ "multipoleOperators/qlmObs", "Qlm mean values"             , "", "MD"}, \
{ "multipoleOperators/beta"  , "beta deformation mean values", "", "V" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief This class analytically calculates the multipole moment mean values.
 */

class MultipoleOperators : public Generic
{
public :

  // Constructors
  MultipoleOperators();                                                // #TEST#
  MultipoleOperators(State);                                           // #TEST#

  //============================================================================

  void calcQlmHO(void);                                                // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  const std::string getNiceInfo(void);                                 // #TEST#

  void calcQlmObs(const Multi<arma::mat> &rho);                        // #TEST#
  void calcBeta(void);                                                 // #TEST#
  void calcNPart(const Multi<arma::mat> &rho);                         // #TEST#
  static double getBetaFromQ20(double parts, double A, double vq20);   // #TEST#
  static double getQ20FromBeta(double parts, double A, double vBeta);  // #TEST#
  void calcRk(void);                                                   // #TEST#
  void calcRkn(void);                                                  // #TEST#
  void calcZk(void);                                                   // #TEST#
  DataTree getDataTree(void);                                          // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// The initial State instance.
  State state;

  /// The \f$\langle m,n,n_z,d,s|\hat{Q}_{l,0}|m',n',n_z',d',s'\rangle\f$ matrices.
  Multi<arma::mat> qlmHO;

  /// The integrated number of particles.
  arma::vec nPart = {0.0, 0.0, 0.0};

  /// The \f$\psi|r_\perp^{2l}|\psi\rangle\f$ values.
  Multi<arma::vec> matR;

  /// The \f$\psi|z^k|\psi\rangle\f$ values.
  Multi<arma::mat> matZ;

  /// The \f$\psi|\hat{Q}_{\lambda,0}|\psi\rangle\f$ values.
  Multi<double> qlmObs;

  /// The \f$\beta\f$ deformation parameter values.
  arma::vec beta;

  /// The \f$\langle m,n|z^k|m',n'\rangle\f$ matrices.
  Multi<arma::mat> matZk;

  /// The \f$\langle m,n|r^{2k}|m',n'\rangle\f$ matrices.
  Multi<arma::mat> matR2k;

  /// The \f$\langle s|s'\rangle\f$ arma::mat.
  arma::mat matS0;

  /// TODO
  Multi<MAT> matRkn;

};

#endif // MULTIPOLE_OPERATORS_H
