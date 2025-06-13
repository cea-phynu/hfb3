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

#include "field_kinetic.h"
#include "global.h"
#include "tools.h"
#include "state.h"
#include "interaction.h"

/** \file
 *  \brief Methods of the FieldKinetic class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > FieldKinetic::validKeys =
  {
  };


//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldKinetic::FieldKinetic(Field::Parameters _parameters, State *_state) : Field(_parameters, _state)
{
  DBG_ENTER;

  name = "kinetic";
  shortName = "Kinet.";
  energyFactor = 2.0;
  mustBeCleared = false;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void FieldKinetic::calcField(void)
{
  DBG_ENTER;

  calculatingLength = -1.0;

  if (!field(NEUTRON, DIRECT).empty() && !field(PROTON, DIRECT).empty()) DBG_LEAVE;

  double startTime = Tools::clock();

  // dependencies
  basis.calcTab();

  // initialization
  field(NEUTRON, DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON, DIRECT  ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  arma::mat &tab = field(NEUTRON, DIRECT);

  INT np;
  double b_r2 = basis.b_r * basis.b_r;
  double b_z2 = basis.b_z * basis.b_z;

  for (INT m = 0; m < basis.mMax; m++)
  {
    for (INT n = 0; n < basis.nMax(m); n++)
    {
      // ***** n' = n.
      // Eq. (F-9), p. F-2, PhD J.-F. Berger.
      np = n;

      for (INT n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
      {
        for (INT d = 0; d < basis.dMax; d++)
        {
          INT a = basis.HOqn.find({m, n, n_z, d, 0}); // spin 0

          for (INT n_zp = 0; n_zp < basis.n_zMax(m, n); n_zp++)
          {
            for (INT dp = 0; dp < basis.dMax; dp++)
            {
              INT b = basis.HOqn.find({m, np, n_zp, dp, 0}); // spin 0
              double Sab = basis.tabzd(n_z, n_zp, d, dp);
              double Sabp2 = basis.tabzd(n_z, n_zp + 2, d, dp);
              double Sabm2 = basis.tabzd(n_z, n_zp - 2, d, dp);
              tab(a, b) = (double(2 * n + m + 1) / b_r2 + (double(n_zp) + 0.5) / b_z2) * Sab;
              tab(a, b) -= 1.0 / (2.0 * b_z2) * (sqrt(n_zp * (n_zp - 1)) * Sabm2 + sqrt((n_zp + 1) * (n_zp + 2)) * Sabp2);

              if (m != 0)
              {
                tab(a + 1, b + 1) = tab(a, b); // spin (1, 1)
              }
            }
          }
        }
      }

      // ***** n' = n - 1
      // Eq. (F-10), p. F-2, PhD J.-F. Berger.
      np = n - 1;

      if (np >= 0)
      {
        for (INT n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
        {
          for (INT d = 0; d < basis.dMax; d++)
          {
            INT a = basis.HOqn.find({m, n, n_z, d, 0}); // spin 0

            for (INT n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++)
            {
              for (INT dp = 0; dp < basis.dMax; dp++)
              {
                INT b = basis.HOqn.find({m, np, n_zp, dp, 0}); // spin 0
                double Sab = basis.tabzd(n_z, n_zp, d, dp);
                tab(a, b) = 1.0 / b_r2 * (sqrt(n * (n + m)) * Sab);

                if (m != 0)
                {
                  tab(a + 1, b + 1) = tab(a, b); // spin (1, 1)
                }
              } // dp
            } // n_zp
          } // d
        } // n_z
      }

      // ***** n' = n + 1
      // Eq. (F-10), p. F-2, PhD J.-F. Berger.
      np = n + 1;

      if (np < basis.nMax(m))
      {
        for (INT n_z = 0; n_z < basis.n_zMax(m, n); n_z++)
        {
          for (INT d = 0; d < basis.dMax; d++)
          {
            INT a = basis.HOqn.find({m, n, n_z, d, 0}); // spin 0
            for (INT n_zp = 0; n_zp < basis.n_zMax(m, np); n_zp++)
            {
              for (INT dp = 0; dp < basis.dMax; dp++)
              {
                INT b = basis.HOqn.find({m, np, n_zp, dp, 0}); // spin 0
                double Sab = basis.tabzd(n_z, n_zp, d, dp);
                tab(a, b) = 1.0 / b_r2 * (sqrt((n + 1) * (n + m + 1)) * Sab);

                if (m != 0)
                {
                  tab(a + 1, b + 1) = tab(a, b); // spin (1, 1)
                }
              } // dp
            } // n_zp
          } // d
        } // n_z
      }
    } // n
  } // m

  // Center of mass correction
  tab *= 0.5 - 0.5 / double(state.sys.nNeut + state.sys.nProt);
  field(PROTON, DIRECT) = field(NEUTRON, DIRECT);

  // Mass factor
  switch(general.compatibility)
  {
    case General::COMPAT_NONE:
      field(NEUTRON, DIRECT) *= HBARC2 / MASS_N;
      field(PROTON , DIRECT) *= HBARC2 / MASS_P;
    break;
    case General::COMPAT_BERGER:
    case General::COMPAT_ROBLEDO:
      field(NEUTRON, DIRECT) *= 41.47;
      field(PROTON , DIRECT) *= 41.47;
    break;
    case General::COMPAT_HFBTHO:
      field(NEUTRON, DIRECT) *= 20.736676229315790 * 2.0;
      field(PROTON , DIRECT) *= 20.736676229315790 * 2.0;
    break;
    default:
      ERROR("Unknown compatibility value.");
  }

  // Store the calculating length
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

