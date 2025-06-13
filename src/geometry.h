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

#ifndef GEOMETRY_H
#define GEOMETRY_H

/** \file
 *  \brief Headers for the Geometry class.
 */

#include "global.h"
#include "generic.h"
#include "state.h"
#include "mesh.h"

/** \brief Compute some geometrical properties of a density, like multipole moments, fragment properties, etc...
 *
 * This class allows to compute some geometrical properties of a density, like multipole moments, fragment properties,
 * by using numerical interpolation. For analytical calculation of the multipole moments, see the Common class.
 */

class Geometry : public Generic
{
public :
  Geometry(void);                                                      // #TEST#
  Geometry(State &_state);                                             // #TEST#
  Geometry(Mesh _mesh, arma::mat _dens, System _system, INT _izNeck = -1);

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  double getOldChargeRMS(double Z, double N);                          // #TEST#
  static arma::cube getDensityMultipoleExpansion(State&, const arma::vec & rVals = arma::linspace(0., 15.0, 151));
  void calcBeta(void);

  //============================================================================
  //============================================================================
  //============================================================================

  /// An instance of a System object.
  System system;

  /// Possible integration types.
  enum
  {
    INT_DENS_L,
    INT_DENS_R,
    INT_DENS_T,
    INT_Q00,
    INT_Q10,
    INT_Q20,
    INT_Q30,
    INT_Q40,
    INT_Q50,
    INT_Q60,
    INT_Z,
    INT_Z2,
    INT_RP2,
    INT_R2,
    INT_R,
    INT_Q_NECK,
  };

  /// The id of the neck.
  INT izNeck = -1;

  /// The position of the neck [fm].
  double neckPos;

  /// The density in the neck [fm^-3].
  double neckDens;

  /// The integral over the left (z < neckPos) subspace.
  double intLeft;

  /// The integral over the right (z > neckPos) subspace.
  double intRight;

  /// The integral over the full space.
  double intTotal;

  /// The \f$\beta\f$ deformation parameter value.
  double beta;

  /// The multipole moments.
  arma::vec qlm;

  /// The mean value of \f$r\f$.
  double int_r;

  /// The mean value of \f$r^2\f$.
  double int_r2;

  /// The mean value of \f$r_\perp^2\f$.
  double int_rp2;

  /// The root mean square (\f$\sqrt{\frac{\INT r^2 \rho d\tau}{\INT \rho d\tau}}\f$.
  double rms;

  /// The mean radius (\f$\sqrt{\frac{\INT r \rho d\tau}{\INT \rho d\tau}}\f$.
  double radius;

  /// The nuclear charge radius (\f$\sqrt{<r_p^2>+<r^2>_p+\frac{3}{4m^2}+\frac{N}{Z}}<r^2>_p)\f$.
  double rc;

  /// The mean value of \f$\hat{Q}_{N} = e^{-\left(\frac{z-z_N}{a_N}\right)^2}\f$.
  double q_neck;

  /** The charge RMS value.
  *
  * J. L. Friar and J. W. Negele
  * Adv. Nucl. Phys. 8 219 (1975)
  */
  double chargeRMS;

  void calc_neck_abcissa(void);
  void calc_neck_abcissa_schunck(void);
  void calc_fragment_properties(void);
  void calcQlm(void);
  void calcAll(void);

  double integrate(INT, INT what = INT_DENS_T);

  /// The mesh used for the calculation.
  Mesh mesh;

  /// The density used for the calculation.
  arma::mat dens;
};

#endif // GEOMETRY_H
