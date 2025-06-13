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

#include "wspot.h"
#include "tools.h"
#include "basis.h"
#include "axis.h"

/** \file
 *  \brief Methods of the WSPot class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constuctor.
 */

WSPot::WSPot(void)
{
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

WSPot::WSPot(INT _nNeut, INT _nProt, Multi<double> _def, Mesh _mesh) : nNeut(_nNeut), nProt(_nProt), nucleusMass(nNeut + nProt)
{
  DBG_ENTER;

  R = arma::vec(2, arma::fill::zeros);
  R(NEUTRON) = r_0(NEUTRON) * pow(nucleusMass, 0.333333333);
  R(PROTON ) = r_0(PROTON ) * pow(nucleusMass, 0.333333333);
  V_0 = arma::vec(2, arma::fill::zeros);
  V_0(NEUTRON) = V(NEUTRON) * (1. - kappa(NEUTRON) * (double)(nNeut - nProt) / nucleusMass);
  V_0(PROTON)  = V(PROTON)  * (1. + kappa(PROTON)  * (double)(nNeut - nProt) / nucleusMass);
  def = _def;
  mesh = _mesh;
  az = Axis(Axis::REGULAR, 257, -32.0, 32.0);

  calcDirectPot();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculation of the direct potential.
 */

void WSPot::calcDirectPot(void)
{
  DBG_ENTER;

  // Temporary objects
  arma::mat x = Tools::matFromCol(mesh.ax.p, mesh.az.nb);
  arma::mat z = Tools::matFromCol(mesh.az.p, mesh.ax.nb).t();
  arma::mat theta = arma::atan2(x, z);
  arma::mat cost1 = arma::cos(theta);
  arma::mat cost2 = cost1 % cost1;
  arma::mat cost3 = cost1 % cost2;
  arma::mat cost4 = cost1 % cost3;
  arma::mat cost5 = cost1 % cost4;
  arma::mat cost6 = cost1 % cost5;
  arma::mat r = arma::sqrt(x % x + z % z);
  arma::mat deformedRadius = arma::ones(mesh.ax.nb, mesh.az.nb);

  // Radial part of spherical harmonics Y_{lm} for l >= 1
  for (auto &key : def.getKeys())
  {
    if (key.front() == 0) deformedRadius += def(key);
    if (key.front() == 1) deformedRadius += def(key) * 1.0 /  2.0 * std::sqrt( 3. / PI) * cost1;
    if (key.front() == 2) deformedRadius += def(key) * 1.0 /  4.0 * std::sqrt( 5. / PI) * (  3.0 * cost2 -   1.0);
    if (key.front() == 3) deformedRadius += def(key) * 1.0 /  4.0 * std::sqrt( 7. / PI) * (  5.0 * cost3 -   3.0 * cost1);
    if (key.front() == 4) deformedRadius += def(key) * 0.3 / 16.0 * std::sqrt( 1. / PI) * ( 35.0 * cost4 -  30.0 * cost2 + 3.0);
    if (key.front() == 5) deformedRadius += def(key) * 1.0 / 16.0 * std::sqrt(11. / PI) * ( 63.0 * cost5 -  70.0 * cost3 +  15.0 * cost1);
    if (key.front() == 6) deformedRadius += def(key) * 1.0 / 32.0 * std::sqrt(13. / PI) * (231.0 * cost6 - 315.0 * cost4 + 105.0 * cost2 - 5.0);
  }

  arma::mat surfNeut = R(NEUTRON) * deformedRadius;
  arma::mat surfProt = R(PROTON ) * deformedRadius;
  directPotentialNeut = V_0(NEUTRON) / ( 1.0 + arma::exp((r - surfNeut) / a(NEUTRON)));
  directPotentialProt = V_0(PROTON ) / ( 1.0 + arma::exp((r - surfProt) / a(PROTON )));

  DBG_LEAVE;
}
