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

#include "global.h"

/** \file
 *  \brief Global instances.
 */

//==============================================================================
//==============================================================================
//==============================================================================

// Default values
bool useBokeh = false;
bool saveToFile = false;
bool useColors = true;
bool forceInvalidDataTree = false;
bool exitOnError = true;
INT  tableStyle = 1;
bool useUtf8 = true;

General general;

Logger logger;

std::set<INT> msgToOut = {MSG_INFO, MSG_MAIN, MSG_ERROR};

//==============================================================================

bool operator==(const UVEC &lhs, const UVEC &rhs)
{
  return arma::all(arma::operator==(lhs, rhs));
}

//==============================================================================

/// Equality operator between IVEC.

bool operator==(const IVEC &lhs, const IVEC &rhs)
{
  return arma::all(arma::operator==(lhs, rhs));
}

//==============================================================================

/// Equality operator between arma::vec.

bool operator==(const arma::vec &lhs, const arma::vec &rhs)
{
  return arma::approx_equal(lhs, rhs, "absdiff", 1e-14);
}

//==============================================================================

/// Equality operator between UMAT.

bool operator==(const UMAT &lhs, const UMAT &rhs)
{
  return arma::all(arma::all(arma::operator==(lhs, rhs)));
}

//==============================================================================

/// Equality operator between IMAT.

bool operator==(const IMAT &lhs, const IMAT &rhs)
{
  return arma::all(arma::all(arma::operator==(lhs, rhs)));
}

//==============================================================================

/// Equality operator between arma::mat.

bool operator==(const arma::mat &lhs, const arma::mat &rhs)
{
  return arma::approx_equal(lhs, rhs, "absdiff", 1e-14);
}

//==============================================================================

/// Equality operator between UCUBE.

bool operator==(const UCUBE &lhs, const UCUBE &rhs)
{
  if (lhs.n_slices != rhs.n_slices) return false;

  for (UINT s = 0; s < lhs.n_slices; s++)
    if (!arma::all(arma::all(arma::operator==(lhs.slice(s), rhs.slice(s))))) return false;

  return true;
}

//==============================================================================

/// Equality operator between ICUBE.

bool operator==(const ICUBE &lhs, const ICUBE &rhs)
{
  if (lhs.n_slices != rhs.n_slices) return false;

  for (UINT s = 0; s < lhs.n_slices; s++)
    if (!arma::all(arma::all(arma::operator==(lhs.slice(s), rhs.slice(s))))) return false;

  return true;
}

//==============================================================================

/// Equality operator between arma::cube.

bool operator==(const arma::cube &lhs, const arma::cube &rhs)
{
  return arma::approx_equal(lhs, rhs, "absdiff", 1e-14);
}

//==============================================================================

/// Inequality operator between UVEC.

bool operator!=(const UVEC &lhs, const UVEC &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between IVEC.

bool operator!=(const IVEC &lhs, const IVEC &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between arma::vec.

bool operator!=(const arma::vec &lhs, const arma::vec &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between UMAT.

bool operator!=(const UMAT &lhs, const UMAT &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between IMAT.

bool operator!=(const IMAT &lhs, const IMAT &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between arma::mat.

bool operator!=(const arma::mat &lhs, const arma::mat &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between UCUBE.

bool operator!=(const UCUBE &lhs, const UCUBE &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between ICUBE.

bool operator!=(const ICUBE &lhs, const ICUBE &rhs)
{
  return !(lhs == rhs);
}

//==============================================================================

/// Inequality operator between arma::cube.

bool operator!=(const arma::cube &lhs, const arma::cube &rhs)
{
  return !(lhs == rhs);
}

// //==============================================================================
//
// double pow(const double v, const INT o)
// {
//   return pow(v, double(o));
// }
//
// //==============================================================================
//
// double pow(const double v, const INT o)
// {
//   return pow(v, double(o));
// }
//
// //==============================================================================
//
// double sqrt(const INT v)
// {
//   return sqrt(double(v));
// }
//
// //==============================================================================
//
// double sqrt(const INT v)
// {
//   return sqrt(double(v));
// }
