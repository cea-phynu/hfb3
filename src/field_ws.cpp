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

#include "field_ws.h"
#include "tools.h"
#include "basis.h"
#include "interaction.h"
#include "discrete.h"

/** \file
 *  \brief Methods of the FieldWS class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 */

FieldWS::FieldWS(Field::Parameters fp, State *_state) : Field(fp, _state)
{
  DBG_ENTER;

  name = "woodsSaxon";
  shortName = "WoodSa";
  nGLA = 96;
  nGHE = 96;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Woods-Saxon field.
 */

void FieldWS::calcField(void)
{
  DBG_ENTER;

  calculatingLength = -1.0;
  double startTime = Tools::clock();
  Discrete discrete(&basis, Mesh::gaussLaguerreHermite(nGLA, nGHE));
  arma::mat localRhon;
  arma::mat localRhop;
  INT da = 0, dc = 0;
  double z_a0 = 0.0;
  double factor = 1.0;
  Mesh mesh;
  arma::vec RRw;
  arma::rowvec Zw;
  arma::mat RRwZZw;
  field(NEUTRON, DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  field(PROTON , DIRECT) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  wspot = WSPot(state.sys.nNeut, state.sys.nProt, def, mesh);

  // precalculation
  arma::mat zValsO = arma::zeros(discrete.mesh.az.nb, basis.n_zGlobalMax);
  arma::mat zValsM = arma::zeros(discrete.mesh.az.nb, basis.n_zGlobalMax);
  arma::mat zValsP = arma::zeros(discrete.mesh.az.nb, basis.n_zGlobalMax);
  arma::vec zPointsShiftedM = discrete.mesh.az.p - 0.5 * basis.d_0 / basis.b_z;
  arma::vec zPointsShiftedP = discrete.mesh.az.p + 0.5 * basis.d_0 / basis.b_z;

  for (INT n_z = 0; n_z < basis.n_zGlobalMax; n_z++)
  {
    zValsO.col(n_z) = basis.zPartNormReduced(discrete.mesh.az.p, n_z);
    zValsM.col(n_z) = basis.zPartNormReduced(zPointsShiftedM, n_z);
    zValsP.col(n_z) = basis.zPartNormReduced(zPointsShiftedP, n_z);
  }

  arma::mat zValsShiftedM;
  arma::mat zValsShiftedP;

  for (INT icase = 0; icase < basis.dMax * basis.dMax; icase++)
  {
    mesh = discrete.mesh;

    switch (icase)
    {
      case 0:
        da = 0;
        dc = 0;
        z_a0 = (0.5 - da) * basis.d_0;
        mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
        mesh.az.p = basis.b_z * discrete.mesh.az.p + z_a0;
        zValsShiftedP = zValsO;
        zValsShiftedM = zValsO;
        factor = 1.0;
        break;

      case 1:
        da = 1;
        dc = 1;
        z_a0 = (0.5 - da) * basis.d_0;
        mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
        mesh.az.p = basis.b_z * discrete.mesh.az.p + z_a0;
        zValsShiftedP = zValsO;
        zValsShiftedM = zValsO;
        factor = 1.0;
        break;

      case 2:
        da = 0;
        dc = 1;
        z_a0 = (0.5 - da) * basis.d_0;
        mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
        mesh.az.p = basis.b_z * discrete.mesh.az.p;
        zValsShiftedM = zValsM;
        zValsShiftedP = zValsP;
        factor = exp(-1.0 * z_a0 * z_a0 / (basis.b_z * basis.b_z));
        break;

      case 3:
        da = 1;
        dc = 0;
        z_a0 = (0.5 - da) * basis.d_0;
        mesh.ax.p = basis.b_r * arma::sqrt(discrete.mesh.ax.p);
        mesh.az.p = basis.b_z * discrete.mesh.az.p;
        zValsShiftedM = zValsP;
        zValsShiftedP = zValsM;
        factor = exp(-1.0 * z_a0 * z_a0 / (basis.b_z * basis.b_z));
        break;
    }

    arma::vec Rw;
    WSPot wspot(state.sys.nNeut, state.sys.nProt, def, mesh);

    // TODO
    // the next loop takes ~0.5s and is the bottleneck for the WS optimizations.
    // Maybe it is possible to speed it up ?

    for (INT m = 0; m < basis.mMax; m++)
    {
      for (INT na = 0; na < basis.nMax(m); na++)
      {
        Rw = factor * basis.rPartNormReduced(discrete.mesh.ax.p, m, na) % discrete.mesh.ax.w;

        for (INT nc = 0; nc < basis.nMax(m); nc++)
        {
          RRw =  Rw % basis.rPartNormReduced(discrete.mesh.ax.p, m, nc);

          for (INT n_za = 0; n_za < basis.n_zMax(m, na); n_za++)
          {
            INT ia = basis.HOqn.find({m, na, n_za, da, 0});
            Zw = (zValsShiftedM.col(n_za) % discrete.mesh.az.w).t();
            double totn;
            double totp;

            for (INT n_zc = 0; n_zc < basis.n_zMax(m, nc); n_zc++)
            {
              INT ic = basis.HOqn.find({m, nc, n_zc, dc, 0});

              if (ic < ia) continue;

              RRwZZw = RRw * (Zw % zValsShiftedP.col(n_zc).t());
              totn = arma::accu(RRwZZw % wspot.directPotentialNeut);
              totp = arma::accu(RRwZZw % wspot.directPotentialProt);
              field(NEUTRON, DIRECT)(ia, ic) = totn;
              field(PROTON , DIRECT)(ia, ic) = totp;

              if (m > 0)
              {
                field(NEUTRON, DIRECT)(ia + 1, ic + 1) = totn;
                field(PROTON , DIRECT)(ia + 1, ic + 1) = totp;
              }
            } // n_zc
          } // n_za
        } // nc
      } // na
    } // m
  } // icase

  field(NEUTRON, DIRECT) = arma::symmatu(field(NEUTRON, DIRECT));
  field(PROTON , DIRECT) = arma::symmatu(field(PROTON , DIRECT));
  calculatingLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

WSPot FieldWS::getWSPot(void)
{
  DBG_ENTER;

  DBG_RETURN(wspot);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Set the deformation parameters.
 */

void FieldWS::setDef(const Multi<double> &_def)
{
  DBG_ENTER;

  def = _def;

  DBG_LEAVE;
}

