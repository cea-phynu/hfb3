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

#include "axis.h"
#include "tools.h"
#include "quadratures.h"

/** \file
 *  \brief Methods of the Axis class.
 */

const std::vector<std::string> Axis::typeName = {"GAUSS_LEGENDRE", "GAUSS_LAGUERRE", "GAUSS_HERMITE", "REGULAR"};

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 * \param _type Type of the axis.
 * \param _nb Number of points.
 * \param vmin Minimal value.
 * \param vmax Maximal value.
 */

Axis::Axis(INT _type, UINT _nb, double vmin, double vmax) : type(_type), nb(_nb)
{
  DBG_ENTER;

  double *n_tab = NULL;
  double *w_tab = NULL;
  double *we_tab = NULL;

  switch (type)  // Be ready for some C++ magic...
  {
    case GAUSS_LEGENDRE:
      switch (nb)
      {
          GLE_ASSIGN // defined in quadratures.h

        default:
          ERROR("gle_p[%d] not found in `quadratures.h`. Please re-generate quadratures.h using quadratures.max or use a different number of nodes", nb);
      }

      break;

    case GAUSS_LAGUERRE:
      switch (nb)
      {
          GLA_ASSIGN // defined in quadratures.h

        default:
          ERROR("gla_p[%d] not found in `quadratures.h`. Please re-generate quadratures.h using quadratures.max or use a different number of nodes", nb);
      }

      break;

    case GAUSS_HERMITE:
      switch (nb)
      {
          GHE_ASSIGN // defined in quadratures.h

        default:
          ERROR("ghe_p[%d] not found in `quadratures.h`. Please re-generate quadratures.h using quadratures.max or use a different number of nodes", nb);
      }

      break;

    case REGULAR:
      break;

    default:
      ERROR("unknown axis type: %d", type);
  }   // It's a kind of magic !!!

  if (type == REGULAR)
  {
    p = arma::zeros(nb);
    w = arma::zeros(nb);
    we = arma::zeros(nb);
    w(0) = 1.0;
    double step = 1.0;

    if (nb != 1)
      step = (vmax - vmin) / (double)(nb - 1);

    for (auto i = 0; i < nb; i++)
    {
      p(i) = vmin + step * (double)i;
      w(i) = step;
      we(i) = step;
    }
  }
  else
  {
    p = arma::zeros(nb);
    w = arma::zeros(nb);
    we = arma::zeros(nb);

    for (auto i = 0; i < nb; i++)
    {
      p[i] = n_tab[i];
      w[i] = w_tab[i];
      we[i] = we_tab[i];
    }
  }

  if (type == GAUSS_LEGENDRE)
  {
    p = (p + 1.0) * (vmax - vmin) / 2.0 + vmin;
    w *= (vmax - vmin) / 2.0;
    we *= (vmax - vmin) / 2.0;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Test equality.
 */

bool Axis::operator==(Axis const &other) const
{
  return ((other.type == type) && (other.nb == nb) && (other.p == p) && (other.w == w) && (other.we == we));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Axis::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"type  ", typeName.at(type)},
      {"min   ", PF("%e", p.min())},
      {"max   ", PF("%e", p.max())},
      {"nb    ", PF("%d", nb)},
    }, isShort);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"Axis", ""},
      {"type  ", typeName.at(type)},
      {"min   ", PF("%e", p.min())},
      {"max   ", PF("%e", p.max())},
      {"nb    ", PF("%d", nb)},
    }, isShort);
  }

  DBG_RETURN(result);
}
