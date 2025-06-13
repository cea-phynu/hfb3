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

#ifndef WSPOT_H
#define WSPOT_H

/** \file
 *  \brief Headers for the WSPot class.
 */

#include "global.h"
#include "mesh.h"
#include "multi.h"

/** \brief Calculate Woods-Saxon potentials (direct and Coulomb).
 *
 * This class calculates the direct WS potential and the Coulomb potential.
 * These potentials appear in some terms of the WS hamiltonian.
 */

class WSPot
{
public:
  // Constructors
  WSPot(void);                                                            // #TEST#
  WSPot(INT _nNeut, INT _nProt, Multi<double> _def, Mesh _mesh);    // #TEST#

  // Potential calculations
  void calcDirectPot(void);                                               // #TEST#

  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;

  //============================================================================

  /// Direct potential for neutrons
  arma::mat directPotentialNeut;

  /// Direct potential for protons
  arma::mat directPotentialProt;

  /// Mesh used for the direct Potential
  Mesh mesh;

  /// z-Axis used for the Coulomb Potential
  Axis az;


private:
  /*
   * WOODS-SAXON PARAMETERS FROM :
   * W. Koepf and P. Ring, Hadrons and Nuclei 339, 81--90 (1991)
   *
   */

  /// Strength of the Woods-Saxon potential
  arma::vec V = {-71.28, -71.28}; // MeV

  /// Real strength of the Woods-Saxon potential
  arma::vec V_0;

  /// Diffusivity of the Woods-Saxon potential
  arma::vec a = {0.615, 0.612}; // fm

  /// Nucleus 'canonical' radius from the liquid drop model
  arma::vec r_0 = {1.233, 1.250}; // fm

  /// Assymmetry of the potential
  arma::vec kappa = {0.462, 0.462}; // Dim.less

  /// Nucleus radius \f$ R = r_0*A**(1./3.)\f$
  arma::vec R;

  /// \f$ c(\alpha)\f$ volume correction coefficient
  double c_alpha;

  /// Center of mass correction
  double z_cm = 0;

  /// Deformations
  Multi<double> def;

  /// Number of neutrons
  INT nNeut;

  /// Number of protons
  INT nProt;

  /// Mass number
  INT nucleusMass;
};

#endif
