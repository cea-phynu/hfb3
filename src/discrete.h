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

#ifndef DISCRETE_H
#define DISCRETE_H

/** \file
 *  \brief Headers for the Discrete class.
 */

#include "state.h"
#include "global.h"
#include "generic.h"
#include "mesh.h"
#include "multi.h"
#include "basis.h"

/** \brief Numerical calculations.
 *
 * This class allows the numerical calculation of some quantities.
 */

class Discrete : public Generic
{
public :

  Discrete(void);                                                      // #TEST#
  explicit Discrete(Basis *_basis, const Mesh &_mesh = Mesh());        // #TEST#

  void calcWaveFunctions(void);                                        // #TEST#

  arma::mat  getWaveFunctionXZ(UINT id);                               // #TEST#
  arma::cube getWaveFunctionXYZ(UINT id);                              // #TEST#

  arma::mat  getLocalXZ(const arma::mat &rho,                          // #TEST#
                        bool ensurePositive);

  arma::cube getLocalXYZ(const arma::mat &rho);                        // #TEST#

  arma::cube getHFmXZ(const arma::mat &matD, UINT id);                 // #TEST#

  arma::mat  getHOXZ(UINT id);                                         // #TEST#

  void calcDensit(const arma::mat &rhonOrig,                           // #TEST#
                  const arma::mat &rhopOrig,
                  UINT ngla,
                  UINT nghe);      // redundant with getLocalXZ ?

  void clear(void);                                                    // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  static arma::vec getLocalDensity(State &state, double x, double y, double z);

  //============================================================================

  /// Auxiliary objects for the spin-orbit field. TODO: rename / move / remove ?
  /// @{
  Multi<arma::mat> Densit;
  Multi<arma::mat> SODensit;
  Multi<arma::mat> BDensit;
  Multi<arma::mat> BSODensit;
  Multi<arma::mat> BSODensit1;
  /// @}

  /// Flag set if density contains negative values.
  bool isNegative = false;

  /// The spatial mesh used for the wavefunctions evaluation.
  Mesh mesh;

  //============================================================================
  //============================================================================
  //============================================================================

  /// A pointer to a Basis instance.
  Basis *basis;

  /// Values of the \f$z\f$-part of the wave functions on the mesh.
  Multi<arma::vec> zVals;

  /// Values of the \f$r_\perp\f$-part of the wave functions on the \f$X\f$-mesh (Cylindrical Solution only).
  Multi<arma::vec> rVals;

  /// Values of the \f$r_\perp\f$-part of the wave functions on the \f$XY\f$-mesh (Cylindrical Solution only).
  Multi<arma::mat> rpVals;

  /// Pre-calculations for getLocalXZ.
  Multi<arma::cube> waveZ;
};

#endif // DISCRETE_H
