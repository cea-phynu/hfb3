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

#include "mesh.h"
#include "tools.h"

/** \file
 *  \brief Methods of the Mesh class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

Mesh::Mesh(void)
{
  DBG_ENTER;
  // default mesh
  setRegular(0, 0, -20, 10, 0, 20, 41, 1, 161);
  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Build a regular mesh.
 *
 * \param min_x The minimum value for \f$x\f$.
 * \param min_y The minimum value for \f$y\f$.
 * \param min_z The minimum value for \f$z\f$.
 * \param max_x The maximum value for \f$x\f$.
 * \param max_y The maximum value for \f$y\f$.
 * \param max_z The maximum value for \f$z\f$.
 * \param nb_x The number of points on the \f$x\f$ axis.
 * \param nb_y The number of points on the \f$y\f$ axis.
 * \param nb_z The number of points on the \f$z\f$ axis.
 */

void Mesh::setRegular(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, UINT nb_x, UINT nb_y, UINT nb_z)
{
  ax = Axis(Axis::REGULAR, nb_x, min_x, max_x);
  ay = Axis(Axis::REGULAR, nb_y, min_y, max_y);
  az = Axis(Axis::REGULAR, nb_z, min_z, max_z);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a regular mesh.
 *
 * \param min_x The minimum value for \f$x\f$.
 * \param min_y The minimum value for \f$y\f$.
 * \param min_z The minimum value for \f$z\f$.
 * \param max_x The maximum value for \f$x\f$.
 * \param max_y The maximum value for \f$y\f$.
 * \param max_z The maximum value for \f$z\f$.
 * \param nb_x The number of points on the \f$x\f$ axis.
 * \param nb_y The number of points on the \f$y\f$ axis.
 * \param nb_z The number of points on the \f$z\f$ axis.
 */

Mesh Mesh::regular(double min_x, double min_y, double min_z, double max_x, double max_y, double max_z, UINT nb_x, UINT nb_y, UINT nb_z)
{
  DBG_ENTER;

  Mesh result;
  result.setRegular(min_x, min_y, min_z, max_x, max_y, max_z, nb_x, nb_y, nb_z);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create a mesh for a POV-Ray file.
 */

Mesh Mesh::df3(void)
{
  DBG_ENTER;

  Mesh result;
  result.setRegular(-10, -10, -20, 10, 10, 20, 32, 32, 64);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Create a (Gauss-Laguerre, 1, Gauss-Hermite) mesh.
 */

Mesh Mesh::gaussLaguerreHermite(UINT nbx, UINT nbz)
{
  DBG_ENTER;

  Axis ax(Axis::GAUSS_LAGUERRE, nbx);
  Axis ay(Axis::REGULAR, 1);
  Axis az(Axis::GAUSS_HERMITE, nbz);
  Mesh result;
  result.ax = ax;
  result.ay = ay;
  result.az = az;

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Test equality.
 */

bool Mesh::operator==(Mesh const &other) const
{
  return ((other.ax == ax) && (other.ay == ay) && (other.az == az));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Mesh::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"ax    ", ax.typeName.at(ax.type)},
      {"ay    ", ay.typeName.at(ax.type)},
      {"az    ", az.typeName.at(ax.type)},
    }, isShort);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"Mesh", ""},
      {"ax    ", ax.info(true)},
      {"ay    ", ay.info(true)},
      {"az    ", az.info(true)},
    }, isShort);
  }

  DBG_RETURN(result);
}


