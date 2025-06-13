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

#ifndef MESH_H
#define MESH_H

/** \file
 *  \brief Headers for the Mesh class.
 */

#include "global.h"
#include "generic.h"
#include "axis.h"

/** \brief Define a spatial mesh.
 *
 * This class defines a spatial mesh, which can be adapted for Gauss-Legendre / -Laguerre / -Hermite integrations.
 */

class Mesh : public Generic
{
public :

  Mesh(void);                                                          // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  static Mesh regular(double,                                          // #TEST#
                      double,
                      double,
                      double,
                      double,
                      double,
                      UINT,
                      UINT,
                      UINT);
  static Mesh gaussLaguerreHermite(UINT, UINT);                        // #TEST#
  static Mesh df3(void);                                               // #TEST#

  bool operator==(Mesh const &other) const;                            // #TEST#

  /** Inequality test.
   *
   *  \param other Other Mesh.
   */

  inline bool operator!=(Mesh const &other) const                      // #TEST#
  {
    return !(*this == other);
  }

  //============================================================================
  //============================================================================
  //============================================================================

  /// The quadrature for \f$x\f$-axis.
  Axis ax;

  /// The quadrature for \f$y\f$-axis.
  Axis ay;

  /// The quadrature for \f$z\f$-axis.
  Axis az;

private:
  void setRegular(double, double, double, double, double, double, UINT, UINT, UINT);
};

#endif // MESH_H
