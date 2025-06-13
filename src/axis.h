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

#ifndef AXIS_H
#define AXIS_H

/** \file
 *  \brief Headers for the Axis class.
 */

#include "global.h"
#include "generic.h"

/** \brief Define an integration scheme for an axis, can be using Regular, Gauss-Legendre, Gauss-Laguerre, or Gauss-Hermite quadratures.
 *
 * This class contains methods to perform numerical integrations using Regular, Gauss-Legendre, Gauss-Laguerre, or Gauss-Hermite quadratures.
 */

class Axis : public Generic
{
public:

  // Constructors
  Axis(INT _type = REGULAR,                                            // #TEST#
       UINT _nb = 10,
       double vmin = -1.0,
       double vmax = 1.0);

  // Misc
  bool operator==(Axis const &other) const;                            // #TEST#

  inline bool operator!=(Axis const &other) const                      // #TEST#
  {
    return !(*this == other);
  }

  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //static Axis AddDiff(Axis axis, double diff);

  //============================================================================
  //============================================================================
  //============================================================================

  /// Possible quadrature types.
  enum
  {
    GAUSS_LEGENDRE,
    GAUSS_LAGUERRE,
    GAUSS_HERMITE,
    REGULAR,
  };

  /// The quadrature type used.
  INT type;

  /// Type names.
  static const std::vector<std::string> typeName;

  /// The number of nodes.
  UINT nb;

  /// The weights.
  arma::vec w;

  /// The nodes.
  arma::vec p;

  /// The weights multiplied by the exponentials.
  arma::vec we;

};

#endif // AXIS_H
