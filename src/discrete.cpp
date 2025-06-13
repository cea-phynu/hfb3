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

#include "discrete.h"
#include "tools.h"
#include "basis.h"

/** \file
 *  \brief Methods of the Discrete class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

Discrete::Discrete()
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Construct a Discrete object from a Basis object.
 *
 *  \param _basis A pointer to a Basis object.
 *  \param _mesh A Mesh object.
 */

Discrete::Discrete(Basis *_basis, const Mesh &_mesh) :  mesh(_mesh), basis(_basis)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the wavefunctions on the mesh.
 */

void Discrete::calcWaveFunctions(void)
{
  DBG_ENTER;

  if (!zVals.empty()) DBG_LEAVE;

  for (INT n_z = 0; n_z < basis->n_zGlobalMax; n_z++)
  {
    for (INT d = 0; d < basis->dMax; d++)
    {
      zVals(n_z, d) = basis->zPart(mesh.az.p, n_z, d);
    }
  }

  INT DMax = (basis->dMax - 1) * 3 + 1;

  for (INT z = 0 ; z < mesh.az.nb ; z++)
  {
    arma::cube &wavecube = waveZ(z);
    wavecube = arma::zeros(basis->n_zGlobalMax, basis->n_zGlobalMax, DMax);

    for (INT d = 0; d < basis->dMax ; d++)
    {
      for (INT dp = 0; dp < basis->dMax ; dp++)
      {
        for (INT n_z = 0; n_z < basis->n_zGlobalMax ; n_z++)
        {

          for (INT n_zp = 0; n_zp < basis->n_zGlobalMax ; n_zp++)
          {
            wavecube(n_z, n_zp, d + 2 * dp) = zVals(n_z, d)(z) * zVals(n_zp, dp)(z);
          }
        }
      }
    }
  }

  INT M = basis->mMax - 1;
  INT N = basis->nGlobalMax - 1;
  INT P = N + (M + 1) / 2;

  for (INT m = 0; m < basis->mMax; m++)
  {
    for (INT n = 0; n < 2 * P + 2 * M + 1; n++)
    {
      rVals(m, n) = basis->rPart(mesh.ax.p, m, n);
    }
  }

  for (INT m = 0; m < basis->mMax; m++)
  {
    for (INT n = 0; n < basis->nMax(m); n++)
    {
      rpVals(m, n) = arma::zeros(mesh.ax.nb, mesh.ay.nb);

      for (INT ix = 0; ix < mesh.ax.nb; ix++)
      {
        for (INT iy = 0; iy < mesh.ay.nb; iy++)
        {
          double x = mesh.ax.p(ix);
          double y = mesh.ay.p(iy);
          rpVals(m, n)(ix, iy) = basis->rPartScalar(sqrt(x * x + y * y), m, n);
        }
      }
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a XY wavefunction on the full mesh.
 *
 *  \param id wavefunction index.
 */

arma::mat Discrete::getWaveFunctionXZ(UINT id)
{
  DBG_ENTER;

  if (zVals.empty()) calcWaveFunctions();

  arma::mat result(mesh.ax.nb, mesh.az.nb);

  INT m = basis->HOqn(0, id);
  INT n = basis->HOqn(1, id);
  INT nz = basis->HOqn(2, id);
  INT d = basis->HOqn(3, id);
  //    INT s = QN(id, 4);
  result = rVals(m, n) * zVals(nz, d).t();

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a XYZ wavefunction on the full mesh.
 *
 *  \param id wavefunction index.
 */

arma::cube Discrete::getWaveFunctionXYZ(UINT id)
{
  DBG_ENTER;

  if (zVals.empty()) calcWaveFunctions();

  arma::cube result(mesh.ax.nb, mesh.ay.nb, mesh.az.nb);

  INT m = basis->HOqn(0, id);
  INT n = basis->HOqn(1, id);
  INT nz = basis->HOqn(2, id);
  INT d = basis->HOqn(3, id);
  //    INT s = QN(id, 4);
  result = Tools::matTimesVec(rpVals(m, n), zVals(nz, d));

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::mat Discrete::getLocalXZ(const arma::mat &rho, bool ensurePositive)
{
  DBG_ENTER;

  if (rho.empty())
  {
    Tools::warning("Empty rho matrix in Discrete::getLocalXZ.");
    DBG_RETURN(arma::mat());
  }

  //dependencies
  calcWaveFunctions();

  arma::mat localRho;

  localRho = arma::zeros(mesh.ax.nb, mesh.az.nb);

  //Builds rho Multi with explicit quantum numbers.
  Multi<arma::cube> Rho; //Rho (m,n_0,n_1)(n_z0,n_z1,d_0+2*d_1)
  INT DMax = (basis->dMax - 1) * 3 + 1;
  Qnumbers &HOqn = basis->HOqn;
  UINT nbBlocks =  basis->HOqn.calcBlocks({0, 1, 3, 4});

  for (INT m = 0; m < basis->mMax; m++)
  {
    for (INT n = 0; n < basis->nMax(m); n++)
    {
      for (INT np = 0; np < basis->nMax(m); np++)
      {
        Rho(m, n, np) = arma::zeros(basis->n_zMax(m, n), basis->n_zMax(m, np), DMax);
      }//np
    }//n
  }//m

  for (UINT i = 0; i < nbBlocks; i++)
  {
    Qnumbers &bras = HOqn.blocks[i];
    INT m = bras(0, 0);
    INT n = bras(1, 0);
    INT d = bras(3, 0);
    INT s = bras(4, 0);

    for (UINT j = 0; j < nbBlocks; j++)
    {
      Qnumbers &kets = HOqn.blocks[j];
      INT mp = kets(0, 0);
      INT np = kets(1, 0);
      INT dp = kets(3, 0);
      INT sp = kets(4, 0);

      if (mp == m && s == sp)
      {
        if (np == n)
        {
          Rho(m, n, np).slice(d + 2 * dp) += 2 * rho.submat(bras.filter, kets.filter);
        }

        if (np != n)
        {
          Rho(m, n, np).slice(d + 2 * dp) += 4 * rho.submat(bras.filter, kets.filter);
        }
      }//if
    }//j
  }//i

  arma::vec zVec;
  arma::vec rVec;

  for (INT m = 0 ; m < basis->mMax ; m++)
  {
    for (INT n = 0 ; n < basis->nMax(m) ; n++)
    {
      for (INT np = 0 ; np < n + 1 ; np++) //Symmetry
      {
        zVec = arma::zeros(mesh.az.nb);
        rVec = rVals(m, n) % rVals(m, np);

        for (INT z = 0 ; z < mesh.az.nb ; z++)
        {
          zVec(z) = arma::accu(waveZ(z).subcube(0, 0, 0, basis->n_zMax(m, n) - 1, basis->n_zMax(m, np) - 1, DMax - 1) % Rho(m, n, np));
        }//z

        localRho += rVec * zVec.t();
      }//np
    }//n
  }//m

  if (ensurePositive)
  {
    if (localRho.min() < 0.0)
    {
      Tools::debug("Negative local density calculated in Discrete::getLocalXZ() ! val: " + std::to_string(localRho.min()));
      isNegative = true;
      localRho = arma::abs(localRho);
    }
  }

  if (localRho.has_nan())
  {
    Tools::warning("local density contains NaN");
  }

  DBG_RETURN(localRho);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and return the local density in the full \f$XYZ\f$ space.
 *
 *  cf. \ref localDensityOptim for optimization details.
 */

arma::cube Discrete::getLocalXYZ(const arma::mat &rho)
{
  DBG_ENTER;

  if (zVals.empty()) calcWaveFunctions();

  arma::cube result;

  result = arma::zeros(mesh.ax.nb, mesh.ay.nb, mesh.az.nb);
  arma::cube temp0(mesh.ax.nb, mesh.ay.nb, mesh.az.nb);
  arma::cube temp1(mesh.ax.nb, mesh.ay.nb, mesh.az.nb);
  arma::cube temp2(mesh.ax.nb, mesh.ay.nb, mesh.az.nb);

  double vrho = 0.0;

  for (INT m = 0; m < basis->mMax; m++)
  {
    for (INT n = 0; n < basis->nMax(m); n++)
    {
      for (INT n_z = 0; n_z < basis->n_zMax(m, n); n_z++)
      {
        for (INT d = 0; d < basis->dMax; d++)
        {
          temp0 = Tools::matTimesVec(rpVals(m, n), zVals(n_z, d)); // cube  = mat * col

          for (INT np = n; np < basis->nMax(m); np++)
          {
            bool n_equal_np = (n == np);
            temp1 = Tools::cubeTimesMat(temp0, rpVals(m, np));           // cube = cube * mat

            for (INT n_zp = 0; n_zp < basis->n_zMax(m, np); n_zp++)
            {
              bool n_z_equal_n_zp = (n_z == n_zp);

              if ((n_equal_np) && (n_zp < n_z)) continue;                          // skip lower triangular part

              for (INT dp = 0; dp < basis->dMax; dp++)
              {
                if ((n_equal_np) && (n_z_equal_n_zp) && (dp < d)) continue;        // skip lower triangular part

                INT a0 = basis->HOqn.find({m, n, n_z, d, 0});
                INT b0 = basis->HOqn.find({m, np, n_zp, dp, 0});
                vrho = rho(a0, b0) * 4;                     // factor 4 for the lower triangular part and the negative time-reversal states

                INT a1 = basis->HOqn.find({m, n, n_z, d, 1});
                INT b1 = basis->HOqn.find({m, np, n_zp, dp, 1});

                if ((a1 >= 0) && (b1 >= 0))
                {
                  vrho += rho(a1, b1) * 4;                     // factor 4 for the lower triangular part and the negative time-reversal states
                }

                if ((n_equal_np) && (n_z_equal_n_zp) && (dp == d))                 // diagonal term has been counted twice
                {
                  vrho /= 2;
                }

                temp2 = Tools::cubeTimesVec(temp1, zVals(n_zp, dp));         // cube = cube * vec
                result += temp2 * vrho;                                            // cube += cube * dbl
              }
            }
          }
        }
      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and return the HO wave function in the \f$XZ\f$-plane.
 */

arma::mat Discrete::getHOXZ(UINT id)
{
  DBG_ENTER;

  arma::mat dens = getWaveFunctionXZ(id);

  DBG_RETURN(dens);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and return the \f$m\f$-contributions to a HF wave function in the \f$XZ\f$-plane.
 */

arma::cube Discrete::getHFmXZ(const arma::mat &matD, UINT id)
{
  DBG_ENTER;

  arma::cube dens = arma::cube(0, 0, 0);

  dens = arma::zeros(mesh.ax.nb, mesh.az.nb, basis->mMax * 2);

  for (UINT a = 0; a < basis->HOqn.nb; a++)
  {
    INT m = basis->HOqn(0, a);
    INT s = basis->HOqn(4, a);
    dens.slice(m * 2 + s) += getWaveFunctionXZ(a) * matD(id, a);
  }

  DBG_RETURN(dens);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the local and spin-orbit densities.
 */

void Discrete::calcDensit(const arma::mat &rhonOrig, const arma::mat &rhopOrig, UINT ngla, UINT nghe)
{
  DBG_ENTER;

  const arma::mat &rhon = rhonOrig;
  const arma::mat &rhop = rhopOrig;

  Mesh mesh = Mesh::gaussLaguerreHermite(ngla, nghe);
  INT Mmax = basis->mMax;
  arma::vec eta = mesh.ax.p; // radial coordinate
  arma::vec zeta = mesh.az.p; // axial coordinate

  arma::mat zetashift = arma::zeros(mesh.az.nb, basis->dMax);

  for (INT dd = 0; dd < basis->dMax; dd++)
  {
    zetashift.col(dd) = zeta - (0.5 - double(dd)) * basis->d_0 / basis->b_z;
  }

  for (INT iiso = 0; iiso < 2; iiso ++)
  {
#ifndef NEW_DENSIT
    BSODensit(iiso)  = arma::zeros(ngla, nghe);
    BDensit(iiso)    = arma::zeros(ngla, nghe);
#endif
    BSODensit1(iiso) = arma::zeros(ngla, nghe);
  }

  double fac2 = 2.*2. / (pow(basis->b_r, 3) * pow(basis->b_z, 2));

#ifndef NEW_DENSIT
  //--------------- Calculation of Densit + SODensit compo 00 -------------------
  double facnorm = pow(basis->b_r, 2) * basis->b_z / 2.;

  arma::mat DensSO_00n = arma::zeros(ngla, nghe);
  arma::mat DensSO_00p = arma::zeros(ngla, nghe);

  for (INT ma = 0; ma < Mmax; ma++ )
  {
    INT mc = ma;
    INT Nmax = basis->nMax(ma);

    for (INT na = 0; na < Nmax; na++ )
    {
      arma::vec Ra = basis->rPartNorm(eta, ma, na);

      for (INT nc = 0; nc < Nmax; nc++ )
      {
        arma::vec Rc = basis->rPartNorm(eta, mc, nc);

        arma::vec Rac = Ra % Rc;
        arma::vec Racm = ma * Rac;
        INT NZmaxa = basis->n_zMax(ma, na);
        INT NZmaxc = basis->n_zMax(mc, nc);

        for (INT da = 0; da < basis->dMax; da++ )
        {
          arma::vec zda = zetashift.col(da);

          for (INT dc = 0; dc < basis->dMax; dc++ )
          {
            arma::vec zdc = zetashift.col(dc);

            for (INT n_za = 0; n_za < NZmaxa; n_za++ )
            {
              UINT ia = basis->HOqn.find({ma, na, n_za, da, 0}); // spin up
              arma::vec Za = basis->zPartNorm(zda, n_za);

              arma::vec Zc;
              arma::vec Zac;
              UINT ic;

              for (INT n_zc = 0; n_zc < NZmaxc; n_zc ++ )
              {
                ic = basis->HOqn.find({mc, nc, n_zc, dc, 0}); // spin up

                if (ic < ia) continue;

                Zac = Za % basis->zPartNorm(zdc, n_zc);

                if (ic != ia) Zac *= 2.0;

                BDensit(0) += ( Rac  * Zac.t() ) * (rhon(ia, ic) + rhon(ia + 1, ic + 1));
                BDensit(1) += ( Rac  * Zac.t() ) * (rhop(ia, ic) + rhop(ia + 1, ic + 1));
                DensSO_00n += ( Racm * Zac.t() ) * (rhon(ia, ic) - rhon(ia + 1, ic + 1));
                DensSO_00p += ( Racm * Zac.t() ) * (rhop(ia, ic) - rhop(ia + 1, ic + 1));

                arma::mat toto = Racm * Zac.t();
              } // n_zc
            } // n_za
          } // dc
        } // da
      } // nc
    } // na
  } // ma


  for (INT iiso = 0; iiso < 2; iiso ++)
  {
    BDensit(iiso) = 2.0 * BDensit(iiso) / facnorm;

    if (BDensit(iiso).min() < 0.0)
    {
      Tools::debug("Negative local density calculated in Discrete::calcDensit()   ! val: " + std::to_string(BDensit(iiso).min()));
      isNegative = true;
      BDensit(iiso) = arma::abs(BDensit(iiso));
    }
  }

  BSODensit(0) = 2.*DensSO_00n / 2. / facnorm;
  BSODensit(1) = 2.*DensSO_00p / 2. / facnorm;
#endif

  //------------------ Calculation of SODensit compo -+ ------------------------
  arma::mat DensSO_mpn = arma::zeros(ngla, nghe);
  arma::mat DensSO_mpp = arma::zeros(ngla, nghe);

  for (INT ma = 1; ma < Mmax; ma++ )
  {
    INT mc = ma - 1;
    UINT Nmaxa = basis->nMax(ma);
    UINT Nmaxc = basis->nMax(mc);

    for (INT na = 0; na < Nmaxa; na++)
    {
      arma::vec Ra = basis->rPartNorm(eta, ma, na);
      arma::vec L0a = basis->rPartL0(eta, ma, na);
      arma::vec Lavecma = basis->rPartLavecm(eta, ma, na);

      for (INT nc = 0; nc < Nmaxc; nc++ )
      {
        arma::vec Rc = basis->rPartNorm(eta, mc, nc);
        arma::vec L0c = basis->rPartL0(eta, mc, nc);
        arma::vec Lavecmc = basis->rPartLavecm(eta, mc, nc);
        arma::vec A_r1 = L0a % Rc;
        arma::vec A_r2 = Ra % L0c;
        arma::vec S_r1 = Lavecma % Rc;
        arma::vec S_r2 = Ra % Lavecmc;
        UINT NZmaxa = basis->n_zMax(ma, na);
        UINT NZmaxc = basis->n_zMax(mc, nc);

        for (INT da = 0; da < basis->dMax; da++ )
        {
          arma::vec zda = zetashift.col(da);

          for (INT dc = 0; dc < basis->dMax; dc++ )
          {
            arma::vec zdc = zetashift.col(dc);

            for (INT n_za = 0; n_za < NZmaxa; n_za++ )
            {
              UINT ia = basis->HOqn.find({ma, na, n_za, da, 1});
              arma::vec Za = basis->zPartNorm(zda, n_za);
              arma::vec dZa = basis->zPartNormd(zda, n_za);

              for (INT n_zc = 0; n_zc < NZmaxc; n_zc ++ )
              {
                UINT ic = basis->HOqn.find({mc, nc, n_zc, dc, 0}); // spin up
                arma::vec Zc = basis->zPartNorm(zdc, n_zc);
                arma::vec dZc = basis->zPartNormd(zdc, n_zc);
                arma::mat Aac = A_r1 * (Za % dZc).t() - A_r2 * (dZa % Zc).t();
                arma::mat Sac = S_r1 * (Za % dZc).t() + S_r2 * (dZa % Zc).t();
                DensSO_mpn += (Sac + Aac) * rhon(ia, ic);
                DensSO_mpp += (Sac + Aac) * rhop(ia, ic);
              } // n_zc
            } // n_za
          } // dc
        } // da
      } // nc
    } // na
  } // ma

  BSODensit1(0) = -fac2 * DensSO_mpn;
  BSODensit1(1) = -fac2 * DensSO_mpp;


  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Clear the object.
 */

void Discrete::clear(void)
{
  DBG_ENTER;

  zVals.clear();
  rVals.clear();
  rpVals.clear();
  Densit.clear();
  SODensit.clear();
  BDensit.clear();
  BSODensit.clear();
  BSODensit1.clear();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Discrete::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  result += Tools::treeStr(
  {
    {"Discrete", ""},
    {"mesh  ", mesh.info(true)},
    {"basis ", basis->info(true)},
  }, isShort);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the local density at a given point in space.
 */

arma::vec Discrete::getLocalDensity(State &state, double x, double y, double z)
{
  DBG_ENTER;

  arma::vec result = arma::zeros(2);
  arma::vec vrho;
  double temp0, temp1, temp2;
  INT a0;
  INT b0;

  double rp = sqrt(x * x + y * y);

  Basis *basis = &(state.basis);

  for (INT m = 0; m < basis->mMax; m++)
  {
    for (INT n = 0; n < basis->nMax(m); n++)
    {
      for (INT n_z = 0; n_z < basis->n_zMax(m, n); n_z++)
      {
        for (INT d = 0; d < basis->dMax; d++)
        {
          a0 = basis->HOqn.find({m, n, n_z, d, 0});

          temp0 = basis->rPartScalar(rp, m, n) * basis->zPartScalar(z, n_z, d);

          for (INT np = n; np < basis->nMax(m); np++)
          {
            bool n_equal_np = (n == np);
            temp1 = temp0 * basis->rPartScalar(rp, m, np);

            for (INT n_zp = 0; n_zp < basis->n_zMax(m, np); n_zp++)
            {
              bool n_z_equal_n_zp = (n_z == n_zp);

              if ((n_equal_np) && (n_zp < n_z)) continue;                          // skip lower triangular part

              for (INT dp = 0; dp < basis->dMax; dp++)
              {
                if ((n_equal_np) && (n_z_equal_n_zp) && (dp < d)) continue;        // skip lower triangular part

                b0 = basis->HOqn.find({m, np, n_zp, dp, 0});

                vrho = arma::vec({state.rho(NEUTRON)(a0, b0), state.rho(PROTON)(a0, b0)}) * 4.0;  // factor 4 for the lower triangular part and the negative time-reversal states

                INT a1 = basis->HOqn.find({m, n, n_z, d, 1});
                INT b1 = basis->HOqn.find({m, np, n_zp, dp, 1});

                if ((a1 >= 0) && (b1 >= 0))
                {
                  vrho += arma::vec({state.rho(NEUTRON)(a1, b1), state.rho(PROTON )(a1, b1)}) * 4.0;
                }

                if ((n_equal_np) && (n_z_equal_n_zp) && (dp == d))                 // diagonal term has been counted twice
                {
                  vrho /= 2;
                }

                temp2 = temp1 * basis->zPartScalar(z, n_zp, dp);
                result += temp2 * vrho;
              }
            }
          }
        }
      }
    }
  }

  DBG_RETURN(result);
}
