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

#include "basis.h"
#include "quadratures.h"
#include "tools.h"

/** \file
 *  \brief Methods of the Basis class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > Basis::validKeys =
  {
    { "basis/d_0"         , "d_0 basis deformation parameter"                , "10.0", "D" },
    { "basis/b_r"         , "b_r basis deformation parameter"                , "1.9" , "D" },
    { "basis/b_z"         , "b_z basis deformation parameter"                , "1.9" , "D" },
    { "basis/nOscil"      , "Number of major HO shells in the basis"         , "11"  , "I" },
    { "basis/n_zMax"      , "Maximum value of n_z in the basis"              , "24"  , "I" },
    { "basis/g_q"         , "g_q basis truncation parameter"                 , "1.0" , "D" },
    { "state/basis/d_0"   , "Previous d_0 basis deformation parameter"       , ""    , "D" },
    { "state/basis/b_r"   , "Previous b_r basis deformation parameter"       , ""    , "D" },
    { "state/basis/b_z"   , "Previous b_z basis deformation parameter"       , ""    , "D" },
    { "state/basis/nOscil", "Previous Number of major HO shells in the basis", ""    , "I" },
    { "state/basis/n_zMax", "Previous Maximum value of n_z in the basis"     , ""    , "I" },
    { "state/basis/g_q"   , "Previous g_q basis truncation parameter"        , ""    , "D" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

Basis::Basis(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constructor from a filename.
 *
 *  \param filename name of the file to be used to instantiate a DataTree and then a Basis instance.
 */

Basis::Basis(const std::string &filename) : Basis(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constructor from a DataTree.
 *
 *  \param dataTree instance of a DataTree to be used to instantiate a Basis instance.
 *  \param prefix prefix to be used when looking for basis parameters in the DataTree instance.
 */

Basis::Basis(const DataTree &dataTree, const std::string &prefix)
{
  DBG_ENTER;

  ASSERT((prefix == "")||(prefix == "state/"), PF("Wrong prefix '%s'", prefix.c_str()));

  dataTree.get(d_0,           prefix + "basis/d_0",    true);
  dataTree.get(b_r,           prefix + "basis/b_r",    true);
  dataTree.get(b_z,           prefix + "basis/b_z",    true);
  dataTree.get(nOscil,        prefix + "basis/nOscil", true);
  dataTree.get(n_zMaxImposed, prefix + "basis/n_zMax", true);
  dataTree.get(g_q,           prefix + "basis/g_q",    true);

  alfBerger = 1 / pow(b_z, 2);
  betBerger = 1 / pow(b_r, 2);
  qBerger   = alfBerger / betBerger;
  hwrBerger = betBerger * 41.47;
  hwzBerger = alfBerger * 41.47;
  hwBerger  = pow((pow(hwrBerger, 2) * hwzBerger), 1.0 / 3);
  cylTruncate();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constructor from basis parameters.
 *
 *  \param _d_0 \f$d_0\f$ parameter.
 *  \param _b_r \f$b_\perp\f$ parameter.
 *  \param _b_z \f$b_z\f$ parameter.
 *  \param _nOscil number of major HO shells.
 *  \param _n_zMaxImposed imposed maximum value for \f$n_z\f$.
 *  \param _g_q \f$Q\f$ truncation parameter of the basis.
 *
 *  For a definition of the cylindrical basis states, cf. \ref basisFunctions. If the basis deformation parameter \f$d_0\f$ is strictly positive, the generated basis is a 1-center basis, otherwise it is a 2-center basis.
 */

Basis::Basis(double _d_0, double _b_r, double _b_z, INT _nOscil, INT _n_zMaxImposed, double _g_q)
{
  DBG_ENTER;

  if (_d_0 < 0) _d_0 = -1.0;

  d_0           = _d_0;
  b_r           = _b_r;
  b_z           = _b_z;
  nOscil        = _nOscil;
  n_zMaxImposed = _n_zMaxImposed;
  g_q           = _g_q;

  alfBerger = 1 / pow(b_z, 2);
  betBerger = 1 / pow(b_r, 2);
  qBerger   = alfBerger / betBerger;
  hwrBerger = betBerger * 41.47;
  hwzBerger = alfBerger * 41.47;
  hwBerger  = pow((pow(hwrBerger, 2) * hwzBerger), 1.0 / 3);
  cylTruncate();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a cylindrical basis from `BERGER2CT` parameters.
 *
 *  \param _nOscil number of oscillator shells.
 *  \param _g_q \f$Q\f$ truncation parameter of the basis.
 *  \param _hw \f$\hbar \omega\f$ parameter.
 *  \param _qzr \f$q_{zr}\f$ parameter.
 *  \param _d_0 \f$d_0\f$ parameter.
 *  \param _n_zMaxImposed imposed maximum value for \f$n_z\f$.
 *
 *  For a definition of the cylindrical basis states, cf. \ref basisStates.
 *  If the basis deformation parameter \f$d_0\f$ is strictly positive, the
 *  generated basis is a 1-center basis, otherwise it is a 2-center basis.
 *  This static method is proposed as a convenient way to instanciate a Basis
 *  object from parameters coming from the `BERGER2CT` solver.
 */

Basis Basis::fromBerger2ct(INT _nOscil, double _g_q, double _hw, double _qzr, double _d_0,
                           INT _n_zMaxImposed)
{
  DBG_ENTER;

  Basis result;

  double betBerger = _hw / (41.47 * pow(_qzr, 1.0 / 3));
  double alfBerger = betBerger * _qzr;
  double _b_r      = 1.0 / sqrt(betBerger);
  double _b_z      = 1.0 / sqrt(alfBerger);

  DBG_RETURN(Basis(_d_0, _b_r, _b_z, _nOscil, _n_zMaxImposed, _g_q));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the cylindrical basis truncation (empirical formula by J.F. Berger):
 *
 * \f[\forall i\ge 1, nzMax[i]=(nn+2).(g_q)^{2/3}+0.5-(g_q).i\f]
 * \f[mx2=\textrm{sup}\left\{i:nzMax[i]\ge 1\right\}\f]
 *
 *  The basis then has the following quantum numbers:
 *
 * \f{eqnarray*}{
 * 0 \le &m& \le mx2-1\\
 * 0 \le &n& \le (mx2-m-1)/2\\
 * 0 \le &n_z& \le nzMax[m+2n]-1
 * \f}
 */

void Basis::cylTruncate(void)
{
  DBG_ENTER;

  mMax         = 0;
  sMax         = 0;
  INT        i = 0;
  IVEC tempN_zMaxBerger;

  while (true)
  {
    i++;
    INT nzm = (INT)(double(nOscil + 2) * pow(g_q, 2.0 / 3) + 0.5 - g_q * double(i));

    if (nzm < 1) break;

    auto size = tempN_zMaxBerger.n_rows;
    tempN_zMaxBerger.resize(size + 1);
    tempN_zMaxBerger(size) = MIN(n_zMaxImposed, nzm);
  }

  mMax            = i - 1;
  nMax            = IVEC(mMax, arma::fill::zeros);
  n_zMaxBerger    = IVEC(tempN_zMaxBerger);
  n_zGlobalMax    = 0;
  nGlobalMax      = 0;
  INT mnGlobalMax = 0;

  for (INT m = 0; m < mMax; m++)
  {
    nMax(m) = (mMax + 1 - m) / 2;

    for (INT n = 0; n < nMax(m); n++)
    {
      if (n >= nGlobalMax) nGlobalMax = n + 1;

      if (m + n >= mnGlobalMax) mnGlobalMax = m + n + 1;
    }
  }

  n_zMax = IMAT(mMax, nGlobalMax, arma::fill::zeros);
  n_zMax.zeros();

  for (INT m = 0; m < mMax; m++)
  {
    for (INT n = 0; n < nMax(m); n++)
    {
      n_zMax(m, n) = n_zMaxBerger(m + 2 * n);

      if (n_zMax(m, n) > n_zGlobalMax) n_zGlobalMax = (int)(n_zMax(m, n));
    }
  }

  dMax = 1;

  if (d_0 > 0.0)
    dMax = 2;
  else
    d_0 = 0.0;

  sMax = 2;
  HOqn = Qnumbers(5);
  std::vector<std::string> names(5);
  names[0] = "m";
  names[1] = "n";
  names[2] = "n_z";
  names[3] = "d";
  names[4] = "s";
  HOqn.setNames(names);
  QNmnn_zd = Qnumbers(4);
  names    = std::vector<std::string>(4);
  names[0] = "m";
  names[1] = "n";
  names[2] = "n_z";
  names[3] = "d";
  QNmnn_zd.setNames(names);
  names    = std::vector<std::string>(2);
  names[0] = "n_z";
  names[1] = "d";
  QNmn     = Qnumbers(2);
  names    = std::vector<std::string>(2);
  names[0] = "m";
  names[1] = "n";
  QNmn.setNames(names);

  for (INT m = 0; m < mMax; m++)
  {
    for (INT n = 0; n < nMax(m); n++)
    {
      QNmn.append({m, n});

      for (INT n_z = 0; n_z < n_zMax(m, n); n_z++)
      {
        for (INT d = 0; d < dMax; d++)
        {
          QNmnn_zd.append({m, n, n_z, d});

          for (INT s = 0; s < sMax; s++)
          {
            if (m - s < 0) continue;

            HOqn.append({m, n, n_z, d, s});
          }
        }
      }
    }
  }

  QNn_zMaxd = Qnumbers(2);
  names     = std::vector<std::string>(2);
  names[0]  = "n_z";
  names[1]  = "d";
  QNn_zMaxd.setNames(names);
  QNn_zMax = Qnumbers(1);
  names    = std::vector<std::string>(1);
  names[0] = "n_z";

  for (INT n_z = 0; n_z < n_zGlobalMax; n_z++)
  {
    for (INT d = 0; d < dMax; d++)
    {
      QNn_zMaxd.append({n_z, d});
    }

    QNn_zMax.append({n_z});
  }

  // Calculate the m block sizes
  mSize = arma::zeros<UVEC>(mMax);

  for (INT m = 0; m < mMax; m++)
  {
    UINT nb = 0;

    for (INT n = 0; n < nMax(m); n++)
    {
      for (INT n_z = 0; n_z < n_zMax(m, n); n_z++)
      {
        for (INT d = 0; d < dMax; d++)
        {
          for (INT s = 0; s < sMax; s++)
          {
            if (m - s < 0) continue;

            nb++;
          }
        }
      }
    }

    mSize(m) = nb;
  }

  // Calculate the omega block indexes (same order as in BERGER2CT)
  for (INT omega = 0; omega < mMax; omega++) omegaIndexHO(omega) = UVEC();

  for (INT m = 0; m < mMax; m++)
  {
    for (INT n = 0; n < nMax(m); n++)
    {
      for (INT n_z = 0; n_z < n_zMax(m, n); n_z++)//J'ai inversé n_z et d...
      {
        for (INT d = 0; d < dMax; d++)
        {
          for (INT s = 0; s < sMax; s++)
          {
            if (m - s < 0) continue;

            INT a     = HOqn.find({m, n, n_z, d, s});
            INT omega = m - s;

            if ((omega >= 0) and (omega < mMax))
            {
              auto size = omegaIndexHO(omega).n_rows;
              omegaIndexHO(omega).resize(size + 1);
              omegaIndexHO(omega)(size) = a;
            }
          }
        }
      }
    }
  }

  DBG_LEAVE;
  return;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree object with the basis parameters.
 */

DataTree Basis::getDataTree(const std::string &prefix) const
{
  DBG_ENTER;

  ASSERT((prefix == "")||(prefix == "state/"), PF("Wrong prefix '%s'", prefix.c_str()));

  DataTree result;

  result.set(prefix + "basis/nOscil", nOscil);
  result.set(prefix + "basis/d_0", d_0);
  result.set(prefix + "basis/b_r", b_r);
  result.set(prefix + "basis/b_z", b_z);
  result.set(prefix + "basis/g_q", g_q);
  result.set(prefix + "basis/n_zMax", n_zMaxImposed);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return minimal information on the basis deformation parameters.
 */

const std::string Basis::niceStr(void) const
{
  DBG_ENTER;

  std::string result = "(";

  result += PF("%sb_r:%s%9.2e%s", Tools::color("yellow").c_str(), Tools::color("blue").c_str(), b_r,
               Tools::color().c_str());
  result += ", ";
  result += PF("%sb_z:%s%9.2e%s", Tools::color("yellow").c_str(), Tools::color("blue").c_str(), b_z,
               Tools::color().c_str());

  if (dMax > 1)
  {
    result += ", ";
    result += PF("%sd_0:%s%9.2e%s", Tools::color("yellow").c_str(), Tools::color("blue").c_str(),
                 d_0, Tools::color().c_str());
  }

  result += ")";

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Basis::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr({{"nOscil", PF("%d", nOscil)},
      {"b_r", PF("%.3f", b_r)},
      {"b_z", PF("%.3f", b_z)},
      {"d_0", PF("%.3f", d_0)},
      {"g_q", PF("%.1f", g_q)},
      {"states", PF("%d", HOqn.nb)},
      {"max_nz", PF("%d", n_zGlobalMax)},
    }, true);
  }
  else
  {
    result += Tools::treeStr({{"Basis", ""},
      {"nOscil", PF("%9d", nOscil)},
      {"b_r   ", PF("%9.6f", b_r)},
      {"b_z   ", PF("%9.6f", b_z)},
      {"d_0   ", PF("%9.6f", d_0)},
      {"g_q   ", PF("%9.6f", g_q)},
      {"States", PF("%9d", HOqn.nb)},
      {"max_nz", PF("%9d", n_zGlobalMax)},

/* #ifdef PRINT_BERGER_BASIS_PARAMETERS */
      {"hw0_Bg", PF("%9.6f", hwBerger)},
      {"qzr_Bg", PF("%9.6f", qBerger)},
      {"alfaBg", PF("%9.6f", alfBerger)},
      {"betaBg", PF("%9.6f", betBerger)},
/* #endif */

    }, false);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the z-part of a normalized wave function without the exponential term:
 */

arma::vec Basis::zPartNormReduced(const arma::vec &zeta, INT n_z)
{
  double fac = 1.0 / sqrt(pow(2, n_z) * fact[n_z] * sqrt(PI));
  return fac * hermite(n_z, zeta);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the z-part of a normalized wave function.
 *
 * For the definition, cf. \ref ssec01_02.
 */

arma::vec Basis::zPartNorm(const arma::vec &zeta, INT n_z)
{
  return zPartNormReduced(zeta, n_z) % arma::exp(-arma::square(zeta) / 2);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the first derivative of the normalized z-part wave function \f$Z_n(\zeta)\f$:
 *
 *  \f{eqnarray*}{
 * \frac{d}{d\zeta} Z_n(\zeta) &=& \sqrt{n/2}Z_{n-1}(\zeta) - \sqrt{(n+1)/2}Z_{n+1}(\zeta) \\
 *                             &=& \sqrt{2n}Z_{n-1}(\zeta) - \zeta Z_n(\zeta)              \\
 *                             &=& \zeta Z_n(\zeta) - \sqrt{2(n+1)}Z_{n+1}(\zeta)
 *  \nonumber
 *  \f}
 *
 * Only the last expression holds for \f$n=0\f$.
 */

arma::vec Basis::zPartNormd(const arma::vec &zeta, INT n_z)
{
  return zeta % zPartNorm(zeta, n_z) - sqrt(2. * (double(n_z) + 1.)) * zPartNorm(zeta, n_z + 1);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the z-part of a wave function \f$\eqref{eq:zPart}\f$ (vectorized version):
 *
 * \f{equation}{
 * Z(z, n_z, d)
 * \equiv
 * \phi_{n_z}(z - z_d)
 * = \frac{1}{\sqrt{b_z}}
 *          \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}}
 *          e^{-\frac{1}{2}(z - z_d)^2/b_z^2}H_{n_z}((z-z_d)/b_z)
 * \f}
 */

arma::vec Basis::zPart(const arma::vec &z, INT n_z, INT d)
{
  double    fac     = 1.0 / sqrt(b_z);
  arma::vec tempVec = (z - (0.5 - double(d)) * d_0) / b_z;
  arma::vec result  = fac * zPartNorm(tempVec, n_z);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the z-part of a wave function \f$\eqref{eq:zPart}\f$ (scalar version):
 *
 * \f{equation}{
 * Z(z, n_z, d)
 * \equiv
 * \phi_{n_z}(z - z_d)
 * = \frac{1}{\sqrt{b_z}}
 *          \frac{1}{\sqrt{2^{n_z} \sqrt{\pi}n_z!}}
 *          e^{-\frac{1}{2}(z - z_d)^2/b_z^2}H_{n_z}((z-z_d)/b_z)
 * \f}
 */

double Basis::zPartScalar(double z, INT n_z, INT d)
{
  double fac1 = 1.0 / sqrt(b_z);
  double fac2 = 1.0 / sqrt(pow(2, n_z) * fact[n_z] * sqrt(PI));
  double zeta = 0.0;

  if (d == -1)
    zeta = z / b_z;
  else
    zeta = (z - (0.5 - double(d)) * d_0) / b_z;

  double herm   = hermite(n_z, zeta);
  double expo   = exp(-zeta * zeta / 2.0);
  double result = fac1 * fac2 * herm * expo;
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the r-part of a normalized wave function WITHOUT THE EXPONENTIAL TERM.
 */

arma::vec Basis::rPartNormReduced(const arma::vec &eta, INT m, INT n)
{
  double fac = sqrt(fact[n] / fact[n + abs(m)]);
  return fac * laguerre(abs(m), n, eta) % arma::pow(arma::sqrt(eta), fabs(double(m)));
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the r-part of a normalized wave function.
 *
 * For the definition, cf. \ref ssec01_02.
 */

arma::vec Basis::rPartNorm(const arma::vec &eta, INT m, INT n)
{
  return rPartNormReduced(eta, abs(m), n) % arma::exp(-eta / 2);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the r-part of a wave function \f$\eqref{eq:rPart}\f$ (vectorized version).
 */

arma::vec Basis::rPart(const arma::vec &r, INT m, INT n)
{
  double    fac     = 1.0 / (b_r * sqrt(PI));
  arma::vec tempVec = arma::square(r / b_r);
  arma::vec result  = fac * rPartNorm(tempVec, abs(m), n);
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the r-part of a wave function \f$\eqref{eq:rPart}\f$ (scalar version).
 */

double Basis::rPartScalar(double r, INT m, INT n)
{
  double fac1   = 1.0 / (b_r * sqrt(PI));
  double fac2   = sqrt(fact[n] / fact[n + abs(m)]);
  double eta    = pow(r / b_r, 2.0);
  double lagu   = laguerre(abs(m), n, eta);
  double powe   = pow(r / b_r, abs(m));
  double expo   = exp(-eta / 2.0);
  double result = fac1 * fac2 * lagu * powe * expo;
  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the 1st derivative of the 2D (r,phi) part of a wave function.
 *
 * TODO: update LaTeX
 * \f[
 * If m>0 and n>0
 *  L^+mn = (2 sqrt{eta}d/d\eta + m/sqrt{eta}) Rmn
 *        = sqrt(n+m) Rm-1,n + sqrt(n+1)Rm-1,n+1
 *
 * L^-mn = (2 sqrt{eta}d/d\eta - m/sqrt{eta}) Rmn
 *       = sqrt(n+m)Rm-1,n - sqrt(n+1)Rm-1,n+1 - 2sqrt(n+m+1)Rm+1,n
 *
 * If n>0
 * L^+-0n = -sqrt(n+1)R1n+sqrt(n)R1,n-1
 *
 * Else
 * L^+-00 = R1,0
 *
 * ipm = 0 => +
 * ipm = 1 => -
 * \f]
 */

arma::vec Basis::rPartNormd(const arma::vec &eta, INT m, INT n, INT ipm)
{
  arma::vec result = arma::zeros(arma::size(eta));

  double vn = double(n);
  double vm = double(m);

  if (m == 0 && n == 0)
  {
    result = rPartNorm(eta, 1, 0) * (-1.);
  }
  else if (m == 0 && n > 0)
  {
    result = rPartNorm(eta, 1, n) * (-sqrt(vn + 1.)) + rPartNorm(eta, 1, n - 1) * sqrt(vn);
  }
  else if (m > 0 && ipm == 0)
  {
    result = rPartNorm(eta, m - 1, n) * sqrt(vn + vm) + rPartNorm(eta, m - 1, n + 1) * sqrt(vn + 1.);
  }
  else if (m > 0 && ipm == 1)
  {
    result = rPartNorm(eta, m - 1, n) * sqrt(vn + vm) - rPartNorm(eta, m - 1, n + 1) * sqrt(vn + 1.) -
             2. * sqrt(vn + vm + 1.) * rPartNorm(eta, m + 1, n);
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** TODO: write documentation.
 */

arma::vec Basis::rPartL0(const arma::vec &eta, INT m, INT n)
{
  arma::vec result = arma::zeros(arma::size(eta));

  double vn = double(n);
  double vm = double(m);

  if (m == 0 && n == 0)
  {
    result = -rPartNorm(eta, 1, 0);
  }
  else if (m == 0 && n > 0)
  {
    result = rPartNorm(eta, 1, n) * (-sqrt(vn + 1.)) - rPartNorm(eta, 1, n - 1) * sqrt(vn);
  }
  else if (m > 0)
  {
    result = rPartNorm(eta, m - 1, n) * sqrt(vn + vm) - rPartNorm(eta, m + 1, n) * sqrt(vn + vm + 1.);
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the r-part of a wave function.
 */

arma::vec Basis::rPartLavecm(const arma::vec &eta, INT m, INT n)
{
  arma::vec result = arma::zeros(arma::size(eta));

  double vn = double(n);
  double vm = double(m);

  if (m != 0)
  {
    result =
      rPartNorm(eta, m - 1, n + 1) * sqrt(vn + 1.) + rPartNorm(eta, m + 1, n) * sqrt(vn + vm + 1.);
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Put a field in msBlock as msBlock(m,s) is a list of every HO states ident in b
 *  of spin s and M = m
 *
 * \param msBlock field where the result is saved
 * \param b Basis of HO1ct or HO2ct states
 * \return 0 if all is ok.
 */

Multi<UVEC> Basis::calcMSBlocks(Basis &b)
{
  DBG_ENTER;

  Multi<UVEC> result;

  for (UINT i = 0; i < b.HOqn.nb; i++)
  {
    INT  s    = b.HOqn(4, i);
    INT  m    = b.HOqn(0, i);
    auto size = result(m, s).n_rows;
    result(m, s).resize(size + 1);
    result(m, s)(size) = i;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the basis distance between two bases.
 */

double Basis::getBasisDistance(const Basis &otherBasis)
{
  DBG_ENTER;

  double sum      = 0;
  double g_q_norm = 1.3;
  double d_0_norm = 10;
  double b_r_norm = 1;
  double b_z_norm = 1;
  sum += pow(g_q - otherBasis.g_q, 2) / pow(g_q_norm, 2);
  sum += pow(d_0 - otherBasis.d_0, 2) / pow(d_0_norm, 2);
  sum += pow(b_r - otherBasis.b_r, 2) / pow(b_r_norm, 2);
  sum += pow(b_z - otherBasis.b_z, 2) / pow(b_z_norm, 2);

  DBG_RETURN(sum);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$W_{ab}\f$ and \f$D_{ab}\f$ matrices.
 *
 *  For a definition of the \f$W_{ab}\f$ and \f$D_{ab}\f$ matrices, cf. \ref sec02_01.
 */

void Basis::calcWDN(void)
{
  DBG_ENTER;

  // dependencies
  calcTab();
  calcHharmo2ct();

  if (!W.empty() && !S.empty() && !N.empty() && !ORtoHO.empty() && !HOtoOR.empty())
  {
    DBG_LEAVE;
  }

  arma::mat eigvec;
  arma::vec eigval;
  double    emaxx = 0.0005;

  if (dMax == 1)
  {
    for (INT m = 0; m < mMax; m++)
    {
      INT n = 0;
      W(m)  = arma::mat(n_zMax(m, n), n_zMax(m, n));
      N(m)  = arma::mat(n_zMax(m, n), n_zMax(m, n));
      H(m)  = arma::mat(n_zMax(m, n), n_zMax(m, n));
      W(m).eye();
      N(m).eye();
      H(m).eye();
    }

    ORqn   = HOqn;
    HOtoOR = arma::mat(ORqn.nb, HOqn.nb);
    HOtoOR.eye();
    ORtoHO = arma::mat(HOqn.nb, ORqn.nb);
    ORtoHO.eye();
    S = arma::zeros(HOqn.nb, HOqn.nb);
    S.eye();
  }
  else     // ----- 2 center basis -----
  {
    UVEC nbn(mMax);
    UVEC nbp(mMax);

    for (INT m = 0; m < mMax; m++)
    {
      INT       n    = 0;
      arma::mat tabj = Tab.submat(0, 0, n_zMax(m, n) - 1, n_zMax(m, n) - 1);

      for (INT i = 1; i < n_zMax(m, n); i += 2) tabj.row(i) *= -1;

      if (!Tools::checkSymmetry(tabj, "tabj in " + std::string(__PRETTY_FUNCTION__)))
        tabj = arma::symmatu(tabj);

      Tools::eig_sym(eigval, eigvec, tabj);

      arma::vec val =
        arma::flipud(eigval);     // eigenvalues: increasing order -> decreasing order
      arma::mat vec = arma::fliplr(eigvec);
      arma::mat u   = vec;

      for (INT i = 1; i < n_zMax(m, n); i += 2) u.row(i) *= -1;

      for (INT i = 0; i < n_zMax(m, n); i++)
        if (val(i) < 0) u.col(i) *= -1;

      arma::mat w(n_zMax(m, n) * 2, 0);
      // I have no idea of what is going on in the following 53 lines,
      // I have simply implemented the exact same algorithm as in HFB2CT.
      // Feel free to document this. N. Dubray
      arma::mat nstates(n_zMax(m, n) * 2, 0);
      arma::mat pstates(n_zMax(m, n) * 2, 0);

      for (INT i = 0; i < n_zMax(m, n); i++)
      {
        double a = (1.0 - fabs(val(i))) * 2.0;

        if (a > emaxx)
        {
          double    s1 = sqrt(a);
          arma::vec tutu(n_zMax(m, n) * 2);
          tutu.subvec(0, u.n_rows - 1)            = u.col(i) / s1;
          tutu.subvec(u.n_rows, u.n_rows * 2 - 1) = -1 * vec.col(i) / s1;

          if (tutu(0) * tutu(n_zMax(m, n)) < 0)
          {
            nstates.resize(nstates.n_rows, nstates.n_cols + 1);
            nstates.col(nstates.n_cols - 1) = tutu;
          }
          else
          {
            pstates.resize(pstates.n_rows, pstates.n_cols + 1);
            pstates.col(pstates.n_cols - 1) = tutu;
          }
        }

        double    s1 = sqrt((1.0 + fabs(val(i))) * 2.0);
        arma::vec tutu(n_zMax(m, n) * 2);
        tutu.subvec(0, u.n_rows - 1)            = u.col(i) / s1;
        tutu.subvec(u.n_rows, u.n_rows * 2 - 1) = vec.col(i) / s1;

        if (tutu(0) * tutu(n_zMax(m, n)) < 0)
        {
          nstates.resize(nstates.n_rows, nstates.n_cols + 1);
          nstates.col(nstates.n_cols - 1) = tutu;
        }
        else
        {
          pstates.resize(pstates.n_rows, pstates.n_cols + 1);
          pstates.col(pstates.n_cols - 1) = tutu;
        }
      }

      w      = arma::join_rows(pstates, nstates);
      nbp(m) = pstates.n_cols;
      nbn(m) = nstates.n_cols;
      W(m)   = w;

      if (!Tools::checkSymmetry(H(m), "H(m) in " + std::string(__PRETTY_FUNCTION__)))
        H(m) = arma::symmatu(H(m));

      arma::mat ham = W(m).t() * H(m) * W(m);

      // INFO("m: %d", m);
      // Tools::info("H(m)", H(m), true);
      // Tools::info("W(m)", W(m), true);
      // Tools::info("ham", ham, true);

      if (!Tools::checkSymmetry(ham, "ham in " + std::string(__PRETTY_FUNCTION__)))
        ham = arma::symmatu(ham);

      Tools::eig_sym(eigval, eigvec, ham);

      // Tools::info("ham", ham);

      arma::vec val2 =
        arma::flipud(eigval);     // eigenvalues: increasing order -> decreasing order
      arma::mat vec2 = arma::fliplr(eigvec);

      N(m)          = vec2;
      arma::mat WDN = W(m) * N(m);
    }

    ORqn = Qnumbers();
    ORqn.appendName("m");
    ORqn.appendName("n");
    ORqn.appendName("o");
    ORqn.appendName("p");
    ORqn.appendName("s");

    for (INT m = 0; m < mMax; m++)
    {
      for (INT n = 0; n < nMax(m); n++)
      {
        // search for an n_zMax equivalent m block.
        INT mEquiv = -1;

        for (INT m2 = 0; m2 < mMax; m2++)
        {
          if (n_zMax(m2, 0) == n_zMax(m, n))
          {
            mEquiv = m2;     // found !
            break;
          }
        }

        ASSERT(mEquiv != -1, "equivalent m not found");

        arma::mat WDN = W(mEquiv) * N(mEquiv);

        for (INT io = 0; io < nbp(mEquiv); io++)
        {
          ORqn.append({m, n, io, 1, 0});

          if (m != 0)
          {
            ORqn.append({m, n, io, 1, 1});
          }
        }

        for (INT io = 0; io < nbn(mEquiv); io++)
        {
          ORqn.append({m, n, io, 0, 0});

          if (m != 0)
          {
            ORqn.append({m, n, io, 0, 1});
          }
        }
      }
    }

    HOtoOR = arma::zeros(ORqn.nb, HOqn.nb);
    ORtoHO = arma::zeros(HOqn.nb, ORqn.nb);

    for (UINT a = 0; a < ORqn.nb; a++)
    {
      INT ma     = ORqn(0, a);
      INT na     = ORqn(1, a);
      INT oa     = ORqn(2, a);
      INT pa     = ORqn(3, a);
      INT sa     = ORqn(4, a);
      INT mEquiv = -1;

      for (INT m2 = 0; m2 < mMax; m2++)
      {
        if (n_zMax(m2, 0) == n_zMax(ma, na))
        {
          mEquiv = m2;     // found !
          break;
        }
      }

      ASSERT(mEquiv != -1, "equivalent m not found");

      arma::mat WDN = W(mEquiv) * N(mEquiv);
      // std::cout << mEquiv << std::endl;
      // Tools::info("W", W(mEquiv));
      // Tools::info("N", N(mEquiv), true);
      // std::cout << std::endl;
      INT d = (pa == 1) ? oa : oa + int(nbp(mEquiv));

      for (UINT b = 0; b < HOqn.nb; b++)
      {
        INT mb = HOqn(0, b);
        INT nb = HOqn(1, b);
        INT sb = HOqn(4, b);

        if (ma != mb) continue;

        if (na != nb) continue;

        if (sa != sb) continue;

        INT n_zb     = HOqn(2, b);
        INT db       = HOqn(3, b);
        auto c       = n_zb + db * n_zMax(mb, nb);
        HOtoOR(a, b) = WDN(c, d);
      }
    }

    S = arma::zeros(HOqn.nb, HOqn.nb);

    for (UINT a = 0; a < HOqn.nb; a++)
    {
      INT ma   = HOqn(0, a);
      INT na   = HOqn(1, a);
      INT n_za = HOqn(2, a);
      INT da   = HOqn(3, a);
      INT sa   = HOqn(4, a);

      for (UINT b = 0; b < HOqn.nb; b++)
      {
        INT mb   = HOqn(0, b);
        INT nb   = HOqn(1, b);
        INT n_zb = HOqn(2, b);
        INT db   = HOqn(3, b);
        INT sb   = HOqn(4, b);

        if (ma != mb) continue;

        if (na != nb) continue;

        if (sa != sb) continue;

        S(a, b) = tabzd(n_za, n_zb, da, db);
      }
    }

    /// Fixing the convention for the sign of each rows
    for (UINT irow = 0; irow < HOtoOR.n_rows; irow++)
    {
      if (HOtoOR.row(irow).index_min() < HOtoOR.row(irow).index_max()) HOtoOR.row(irow) *= -1.;
    }

    // Tools::info("HOtoOR", HOtoOR);
    // Inverse transformation
    ORtoHO = S * HOtoOR.t();
  }

  // Calculate the omega block indices

  omegaIndexOR.clear();

  for (UINT a = 0; a < ORqn.nb; a++)
  {
    INT m     = ORqn(0, a);
    INT s     = ORqn(4, a);
    INT omega = m - s;

    if ((omega >= 0) and (omega < mMax))
    {
      auto size = omegaIndexOR(omega).n_rows;
      omegaIndexOR(omega).resize(size + 1);
      omegaIndexOR(omega)(size) = a;
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$H_{ab}\f$.
 *
 *  For a definition of the \f$H_{ab}\f$ matrices, cf. \ref sec02_01.
 */

void Basis::calcHharmo2ct(void)
{
  DBG_ENTER;

  // dependencies
  calcTab();
  calcRecouv();

  if (!H.empty())
  {
    DBG_LEAVE;
  }

  if (dMax != 2)
  {
    DBG_LEAVE;
  }

  for (INT m = 0; m < mMax; m++)
  {
    INT n          = 0;
    H(m)           = arma::mat(n_zMax(m, n) * 2, n_zMax(m, n) * 2);
    arma::mat &mat = H(m);

    for (INT nz = 0; nz < n_zMax(m, n); nz++)
    {
      for (INT nzp = 0; nzp < n_zMax(m, n); nzp++)
      {
        for (INT d = 0; d < dMax; d++)
        {
          for (INT dp = 0; dp < dMax; dp++)
          {
            double di     = (0.5 - double(d)) * d_0;
            double dj     = (0.5 - double(dp)) * d_0;
            double hzElem = (((2 * double(n) + double(m) + 1) / pow(b_r, 2) + (double(nzp) + 0.5) / pow(b_z, 2) +
                              (d_0 * 0.5 + di) * (d_0 * 0.5 + dj) / (2 * pow(b_z, 4))) *
                             tabzd(nz, nzp, d, dp) +
                             (d_0 * 0.5 + dj) / (pow(b_z, 3) * sqrt(2)) *
                             (sqrt(nz) * tabzd(nz - 1, nzp, d, dp) +
                              sqrt(nzp) * tabzd(nz, nzp - 1, d, dp)) -
                             d_0 / pow(b_z, 4) * Recouvz(nz * 2 + d, nzp * 2 + dp)) *
                            2 * 20.735;
            mat(d * n_zMax(m, n) + nz, dp * n_zMax(m, n) + nzp) = hzElem;
          }
        }
      }
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the overlaps.
 */

void Basis::calcRecouv(void)
{
  DBG_ENTER;

  if (!Recouv.empty() && !Recouvz.empty())
  {
    DBG_LEAVE;
  }

  if (dMax == 1)
  {
    DBG_LEAVE;
  }

  Recouv = arma::mat(QNn_zMaxd.nb, QNn_zMaxd.nb);     // Eq. (E-17), p. E-3, PhD J.-F. Berger.

  for (INT d = 0; d < dMax; d++)
  {
    for (INT dp = 0; dp < dMax; dp++)
    {
      double di     = (0.5 - double(d)) * d_0;
      double dj     = (0.5 - double(dp)) * d_0;
      Recouv(d, dp) = 0.5 * exp(-1.0 * pow((di - dj) / (2.0 * b_z), 2)) *
                      erfc(-1.0 * (di + dj) /
                           (2.0 * b_z));     // Modified !? Eq. (E-18), p. E-3, PhD J.-F. Berger.

      for (INT nzj = 1; nzj < n_zGlobalMax; nzj++)
      {
        Recouv(d, nzj * 2 + dp) = 0;
        Recouv(d, nzj * 2 + dp) += (di - dj) / (b_z * sqrt(2)) * Recouv(d, (nzj - 1) * 2 + dp);
        Recouv(d, nzj * 2 + dp) += b_z / sqrt(2) * zPartScalar(0.0, 0, d) * zPartScalar(0.0, nzj - 1, dp);
        Recouv(d, nzj * 2 + dp) /= sqrt(nzj);
      }
    }
  }

  for (INT d = 0; d < dMax; d++)
  {
    for (INT dp = 0; dp < dMax; dp++)
    {
      for (INT nzj = 1; nzj < n_zGlobalMax; nzj++)
      {
        Recouv(nzj * 2 + dp, d) = Recouv(d, nzj * 2 + dp);
      }
    }
  }

  for (INT d = 0; d < dMax; d++)
  {
    for (INT dp = 0; dp < dMax; dp++)
    {
      double di = (0.5 - double(d)) * d_0;
      double dj = (0.5 - double(dp)) * d_0;

      for (INT nzi = 1; nzi < n_zGlobalMax; nzi++)
      {
        for (INT nzj = 1; nzj < n_zGlobalMax; nzj++)
        {
          Recouv(nzi * 2 + d, nzj * 2 + dp) = 0;
          Recouv(nzi * 2 + d, nzj * 2 + dp) +=
            (di - dj) / (b_z * sqrt(2)) * Recouv(nzi * 2 + d, (nzj - 1) * 2 + dp);
          Recouv(nzi * 2 + d, nzj * 2 + dp) +=
            b_z / sqrt(2) * zPartScalar(0.0, nzi, d) * zPartScalar(0.0, nzj - 1, dp);
          Recouv(nzi * 2 + d, nzj * 2 + dp) +=
            sqrt(nzi) * Recouv((nzi - 1) * 2 + d, (nzj - 1) * 2 + dp);
          Recouv(nzi * 2 + d, nzj * 2 + dp) /= sqrt(nzj);
        }
      }
    }
  }

  for (INT d = 0; d < dMax; d++)
  {
    for (INT dp = 0; dp < dMax; dp++)
    {
      for (INT nzi = 1; nzi < n_zGlobalMax; nzi++)
      {
        for (INT nzj = 1; nzj < n_zGlobalMax; nzj++)
        {
          Recouv(nzj * 2 + dp, nzi * 2 + d) = Recouv(nzi * 2 + d, nzj * 2 + dp);
        }
      }
    }
  }

  Recouvz = arma::mat(QNn_zMaxd.nb, QNn_zMaxd.nb);     // Eq. (E-16), p. E-2, PhD J.-F. Berger.

  for (INT d = 0; d < dMax; d++)
  {
    for (INT dp = 0; dp < dMax; dp++)
    {
      double di = (0.5 - double(d)) * d_0;
      double dj = (0.5 - double(dp)) * d_0;

      for (INT nzi = 0; nzi < n_zGlobalMax; nzi++)
      {
        for (INT nzj = 0; nzj < n_zGlobalMax; nzj++)
        {
          Recouvz(nzi * 2 + d, nzj * 2 + dp) = 0;

          if (nzi > 0)
            Recouvz(nzi * 2 + d, nzj * 2 + dp) +=
              sqrt(nzi) * Recouv((nzi - 1) * 2 + d, nzj * 2 + dp);

          if (nzj > 0)
            Recouvz(nzi * 2 + d, nzj * 2 + dp) +=
              sqrt(nzj) * Recouv(nzi * 2 + d, (nzj - 1) * 2 + dp);

          Recouvz(nzi * 2 + d, nzj * 2 + dp) +=
            (di + dj) / sqrt(2) / b_z * Recouv(nzi * 2 + d, nzj * 2 + dp);
          Recouvz(nzi * 2 + d, nzj * 2 + dp) +=
            zPartScalar(0.0, nzi, d) * zPartScalar(0.0, nzj, dp) * b_z / sqrt(2);
          Recouvz(nzi * 2 + d, nzj * 2 + dp) /= sqrt(2) / b_z;
        }
      }
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convenience function.
 */

double Basis::tabzd(INT nza, INT nzb, INT da, INT db) const
{
  if (nza < 0) return 0;

  if (nzb < 0) return 0;

  if (da == db)
  {
    return (nza == nzb) ? 1 : 0;
  }

  return (da < db) ? Tab(nza, nzb) : Tab(nzb, nza);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check equality between Basis instances.
 */

bool Basis::operator==(const Basis &other) const
{
  double maxErr = 1e-8;

  if (fabs(b_r - other.b_r) > maxErr) return false;

  if (fabs(b_z - other.b_z) > maxErr) return false;

  if (fabs(d_0 - other.d_0) > maxErr) return false;

  if (fabs(g_q - other.g_q) > maxErr) return false;

  if (nOscil != other.nOscil) return false;

  if (n_zMaxImposed != other.n_zMaxImposed) return false;

  return true;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Talman coefficients for the \f$z\f$-part.
 *
 * Appendix D, PhD J.-F. Berger.
 * J.-P. Ebran's derivations, p. 69, Eq. (IV-A.29).
 * Eq. (2.56), p. 16, P. Descourt's report.
 */

void Basis::calcTalmanz(void)
{
  DBG_ENTER;
  // dependencies
  calcTab();

  if (!talmanz.empty())
  {
    DBG_LEAVE;
  }

  // The mnz vector is generated by Maxima (better accuracy).

  for (INT n_za = 0; n_za < n_zGlobalMax + 1; n_za++) // TODO: check bounds (it was n_zGlobalMax + 2 at some point)
  {
    for (INT da = 0; da < dMax; da++)
    {
      for (INT n_zb = 0; n_zb < n_zGlobalMax + 1; n_zb++) // TODO: check bounds (it was n_zGlobalMax + 2 at some point)
      {
        for (INT db = 0; db < dMax; db++)
        {
          talmanz(n_za, da, n_zb, db) = arma::zeros(n_zGlobalMax * 2 + 2); // TODO: check bounds
          arma::vec &talvec           = talmanz(n_za, da, n_zb, db);

          INT n_zcMin = 0;
          INT n_zcMax = n_za + n_zb + 1;

          for (INT n_zc = n_zcMin; n_zc < n_zcMax; n_zc++)
          {
            INT    n_zdMin = MAX(n_zc - n_zb, 0);
            INT    n_zdMax = MIN(n_zc, n_za) + 1;
            double factor  = 1.0 / (mnz[n_za] * mnz[n_zb] * mnz[n_zc]);
            double sum     = 0;

            for (INT n_zd = n_zdMin; n_zd < n_zdMax; n_zd++)
            {
              sum += pow(mnz[n_zd] * mnz[n_zc - n_zd], 2) * mnz[n_za - n_zd] *
                     mnz[n_zb - n_zc + n_zd] * tabzd(n_za - n_zd, n_zb - n_zc + n_zd, da, db);
            }

            talvec(n_zc) = sum * factor;
          }     // n_zc
        }       // db
      }         // n_zb
    }           // da
  }             // n_za

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

// /** Reduced expression of the Moshinsky coefficients for the \f$r\f$-part.
// */
//
// double Basis::getMoshinskyr(INT m1, INT n1, INT m2, INT n2, INT N2) const
// {
//   if (m1 + m2 != 0) return 0.0;
//
//   if (2 * N2 != 2 * n1 + abs(m1) + 2 * n2 + abs(m2)) return 0.0;
//
//   INT    v1     = n1 + (abs(m1) - m1) / 2;
//   INT    v2     = n2 + (abs(m2) - m2) / 2;
//   INT    V2     = N2;
//   INT    e1     = 2 * n1 + abs(m1);
//   INT    e2     = 2 * n2 + abs(m2);
//   double sign0  = pow(-1.0, N2 - n1 - n2);
//   double sign1  = pow(-1.0, V2);
//   double fac0a  = fact[n1 + abs(m1)] * fact[n1];
//   double fac0b  = fact[n2 + abs(m2)] * fact[n2];
//   double fac0d  = fact[N2] * fact[N2];
//   double fac0   = sqrt(fac0a * fac0b / fac0d);
//   double sum    = 0.0;
//   double sdelta = sqrt(2.0) / 2.0;
//
//   for (INT g1 = 0; g1 < v1 + 1; g1++)
//   {
//     for (INT g2 = MAX(0, v2 - v1 + g1 - m1); g2 < v2 + 1; g2++)
//     {
//       for (INT b1 = 0; b1 < MIN(V2, g1) + 1; b1++)
//       {
//         for (INT b2 = 0; b2 < MIN(V2, g2) + 1; b2++)
//         {
//           for (INT b3 = MAX(0, V2 - v2 - b1 - b2 + g2);
//                b3 < MIN(MIN(V2, v1 - g1), V2 - b1 - b2) + 1; b3++)
//           {
//             for (INT b = MAX(0, v1 - v2 + m1 + g2 - g1); b < MIN(v1 - v2 + m1 + g2 - g1, 0) + 1;
//                  b++)
//             {
//               double sign2 = pow(-1.0, b1 - b2 + b);
//               double bino0 = logMultinomial(g1 - b1, g2 - b2, v1 - g1 - b3, v2 - g2 - V2 + b1 + b2 + b3);
//               double bino1 = logMultinomial(b1, b2, b3, V2 - b1 - b2 - b3);
//               sum += sign2 * exp(bino0 + bino1);
//             }     // b
//           }       // b3
//         }         // b2
//       }           // b1
//     }             // g2
//   }               // g1
//
//   double delts = pow(sdelta, e2);
//   double deltc = pow(sdelta, e1);
//   double res   = sum * sign0 * sign1 * fac0 * delts * deltc;
//   return res;
// }

//==============================================================================
//==============================================================================
//==============================================================================

/** Full expression of the Moshinsky coefficients for the \f$z\f$-part (TESTS ONLY).
 */

Multi<arma::vec> Basis::getFullMoshinskyz(void)
{
  DBG_ENTER;

  //dependencies
  calcTab();

  Multi<arma::vec> fullMoshinskyz;
  // Coefficients for the z-part
  /* mnz() Generated by Maxima
  arma::vec mnz = arma::zeros(n_zGlobalMax * 2);

  for (INT n_z = 0; n_z < n_zGlobalMax * 2; n_z++)
  {
    mnz(n_z) = sqrt(pow(2.0, n_z) / fact[n_z]);
  }
  */
  // size:
  // moshinskyz(2NZ, 2NZ)( 4NZ)

  for (INT n_za = 0; n_za < n_zGlobalMax * 2 + 1; n_za++)
  {
    for (INT n_zb = 0; n_zb < n_zGlobalMax * 2 + 1; n_zb++)
    {
      fullMoshinskyz(n_za, n_zb) = arma::zeros(n_zGlobalMax * 4 + 1);
      arma::vec &mosvec      = fullMoshinskyz(n_za, n_zb);

      for (INT n_zc = 0; n_zc < n_za + n_zb + 1; n_zc++)
      {
        INT    n_zd    = n_za + n_zb - n_zc;
        INT    n_zeMin = MAX(n_zb - n_zc, 0);
        INT    n_zeMax = MIN(n_zb, n_zd) + 1;
        double factor  = pow(2.0, (n_zc + n_zd) / -2.0) /
                         (mnz[n_za] * mnz[n_zb] * mnz[n_zc] * mnz[n_zd]);
        double sum = 0;


        for (INT n_ze = n_zeMin; n_ze < n_zeMax; n_ze++)
        {
          sum += pow(mnz[n_zc - n_zb + n_ze] * mnz[n_zb - n_ze] * mnz[n_zd - n_ze] * mnz[n_ze], 2) *
                 pow(-1.0, n_ze);
        }//n_ze

        mosvec(n_zc) = sum * factor;
      }// n_zc
    }// n_zb
  }// n_za

  DBG_RETURN(fullMoshinskyz);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the log of a multinomial expression.
 */

double Basis::logMultinomial(INT a, INT b, INT c, INT d) const
{
  double loga = 0.0;
  double logb = 0.0;
  double logc = 0.0;
  double logd = 0.0;

  for (INT i = 2; i <= a; i++) loga += log((double)(i));

  for (INT i = 2; i <= b; i++) logb += log((double)(i));

  for (INT i = 2; i <= c; i++) logc += log((double)(i));

  for (INT i = 2; i <= d; i++) logd += log((double)(i));

  double result = -loga - logb - logc - logd;

  for (INT i = 2; i <= a + b + c + d; i++) result += log((double)(i));

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Full expression of the Moshinsky coefficients for the \f$r\f$-part (TESTS ONLY).
 */

double Basis::getFullMoshinskyr(INT m1, INT n1, INT m2, INT n2, INT M1, INT N1, INT M2,
                                INT N2) const
{
  if (M1 + M2 != m1 + m2) return 0.0;

  if (2 * N1 + abs(M1) + 2 * N2 + abs(M2) != 2 * n1 + abs(m1) + 2 * n2 + abs(m2)) return 0.0;

  INT    v1     = n1 + (abs(m1) - m1) / 2;
  INT    v2     = n2 + (abs(m2) - m2) / 2;
  INT    V2     = N2 + (abs(M2) - M2) / 2;
  INT    e1     = 2 * n1 + abs(m1);
  INT    e2     = 2 * n2 + abs(m2);
  double sign0  = pow(-1.0, N1 + N2 - n1 - n2);
  double sign1  = pow(-1.0, V2 + M2);
  double fac0a  = fact[n1 + abs(m1)] * fact[n1];
  double fac0b  = fact[n2 + abs(m2)] * fact[n2];
  double fac0c  = fact[N1 + abs(M1)] * fact[N1];
  double fac0d  = fact[N2 + abs(M2)] * fact[N2];
  double fac0   = sqrt(fac0a * fac0b / (fac0c * fac0d));
  double sum    = 0.0;
  double sdelta = sqrt(2.0) / 2.0;

  if ((M1 >= 0) && (M2 >= 0))
  {
    for (INT g1 = 0; g1 < v1 + 1; g1++)
    {
      for (INT g2 = MAX(0, v2 - v1 + g1 - m1); g2 < v2 + 1; g2++)
      {
        for (INT b1 = 0; b1 < MIN(V2, g1) + 1; b1++)
        {
          for (INT b2 = 0; b2 < MIN(V2, g2) + 1; b2++)
          {
            for (INT b3 = MAX(0, V2 - v2 - b1 - b2 + g2);
                 b3 < MIN(MIN(V2, v1 - g1), V2 - b1 - b2) + 1; b3++)
            {
              for (INT b = MAX(0, v1 - v2 + m1 + g2 - g1 - M1);
                   b < MIN(v1 - v2 + m1 + g2 - g1, M2) + 1; b++)
              {
                double sign2 = pow(-1.0, b1 - b2 + b);
                double bino0 =
                  logMultinomial(g1 - b1, g2 - b2, v1 - g1 - b3, v2 - g2 - V2 + b1 + b2 + b3);
                double bino1 = logMultinomial(b1, b2, b3, V2 - b1 - b2 - b3);
                double bino2 = logBinomial(M1, v1 - v2 + m1 + g2 - g1 - b);
                double bino3 = logBinomial(M2, b);
                sum += sign2 * exp(bino0 + bino1 + bino2 + bino3);
              }     // b
            }       // b3
          }         // b2
        }           // b1
      }             // g2
    }               // g1
  }
  else if ((M1 < 0) && (M2 >= 0))
  {
    for (INT g1 = 0; g1 < v1 + 1; g1++)
    {
      for (INT g2 = MAX(0, v2 - v1 + g1 - m1); g2 < v2 + 1; g2++)
      {
        for (INT b1 = 0; b1 < MIN(V2, g1) + 1; b1++)
        {
          for (INT b2 = 0; b2 < MIN(V2, g2) + 1; b2++)
          {
            for (INT b3 = MAX(0, V2 - v2 - b1 - b2 + g2);
                 b3 < MIN(MIN(V2, v1 - g1), V2 - b1 - b2) + 1; b3++)
            {
              for (INT b = 0; b < MIN(v1 - v2 + m1 + g2 - g1, M2) + 1; b++)
              {
                double sign2 = pow(-1.0, b1 - b2 + b);
                double bino0 =
                  logMultinomial(g1 - b1, g2 - b2, v1 - g1 - b3, v2 - g2 - V2 + b1 + b2 + b3);
                double bino1 = logMultinomial(b1, b2, b3, V2 - b1 - b2 - b3);
                double bino2 = logBinomial(M1, v1 - v2 + m1 + g2 - g1 - b);
                double bino3 = logBinomial(M2, b);
                double signb = binomialSign(M1, v1 - v2 + m1 + g2 - g1 - b);
                sum += sign2 * exp(bino0 + bino1 + bino2 + bino3) * signb;
              }     // b
            }       // b3
          }         // b2
        }           // b1
      }             // g2
    }               // g1
  }
  else if ((M1 >= 0) && (M2 < 0))
  {
    for (INT g1 = 0; g1 < v1 + 1; g1++)
    {
      for (INT g2 = MAX(0, v2 - v1 + g1 - m1); g2 < v2 + 1; g2++)
      {
        for (INT b1 = 0; b1 < MIN(V2, g1) + 1; b1++)
        {
          for (INT b2 = 0; b2 < MIN(V2, g2) + 1; b2++)
          {
            for (INT b3 = MAX(0, V2 - v2 - b1 - b2 + g2);
                 b3 < MIN(MIN(V2, v1 - g1), V2 - b1 - b2) + 1; b3++)
            {
              for (INT b = MAX(0, v1 - v2 + m1 + g2 - g1 - M1); b < v1 - v2 + m1 + g2 - g1 + 1; b++)
              {
                double sign2 = pow(-1.0, b1 - b2 + b);
                double bino0 =
                  logMultinomial(g1 - b1, g2 - b2, v1 - g1 - b3, v2 - g2 - V2 + b1 + b2 + b3);
                double bino1 = logMultinomial(b1, b2, b3, V2 - b1 - b2 - b3);
                double bino2 = logBinomial(M1, v1 - v2 + m1 + g2 - g1 - b);
                double bino3 = logBinomial(M2, b);
                double signb = binomialSign(M2, b);
                sum += sign2 * exp(bino0 + bino1 + bino2 + bino3) * signb;
              }     // b
            }       // b3
          }         // b2
        }           // b1
      }             // g2
    }               // g1
  }
  else if ((M1 < 0) && (M2 < 0))
  {
    for (INT g1 = 0; g1 < v1 + 1; g1++)
    {
      for (INT g2 = MAX(0, v2 - v1 + g1 - m1); g2 < v2 + 1; g2++)
      {
        for (INT b1 = 0; b1 < MIN(V2, g1) + 1; b1++)
        {
          for (INT b2 = 0; b2 < MIN(V2, g2) + 1; b2++)
          {
            for (INT b3 = MAX(0, V2 - v2 - b1 - b2 + g2);
                 b3 < MIN(MIN(V2, v1 - g1), V2 - b1 - b2) + 1; b3++)
            {
              for (INT b = 0; b < v1 - v2 + m1 + g2 - g1 + 1; b++)
              {
                double sign2 = pow(-1.0, b1 - b2 + b);
                double bino0 =
                  logMultinomial(g1 - b1, g2 - b2, v1 - g1 - b3, v2 - g2 - V2 + b1 + b2 + b3);
                double bino1 = logMultinomial(b1, b2, b3, V2 - b1 - b2 - b3);
                double bino2 = logBinomial(M1, v1 - v2 + m1 + g2 - g1 - b);
                double bino3 = logBinomial(M2, b);
                double signb = binomialSign(M2, b);
                double signc = binomialSign(M1, v1 - v2 + m1 + g2 - g1 - b);
                sum += sign2 * exp(bino0 + bino1 + bino2 + bino3) * signb * signc;
              }     // b
            }       // b3
          }         // b2
        }           // b1
      }             // g2
    }               // g1
  }

  double delts = pow(sdelta, e2);
  double deltc = pow(sdelta, e1);
  double res   = sum * sign0 * sign1 * fac0 * delts * deltc;

  return res;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the sign of a binomial expression (TESTS ONLY).
 */

double Basis::binomialSign(INT n, INT m) const
{
  if (n >= 0)
  {
    if ((m >= 0) && (m <= n))
      return 1;
    else
      return 0;
  }
  else
  {
    if (m >= 0)
      return pow(-1.0, m % 2);
    else if (m <= n)
      return pow(-1.0, (n - m) % 2);
    else
      return 0;
  }

  return 0;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the log of a binomial expression (TESTS ONLY).
 */

double Basis::logBinomial(INT n, INT m) const
{
  double result = 0.0;

  if (n >= 0)
  {
    if ((m >= 0) && (m <= n))
    {
      for (INT i = 1; i <= m; i++) result += log((double)(n - m + i)) - log((double)(i));
    }
    else
      return result;
  }
  else
  {
    if (m >= 0)
    {
      for (INT i = 1; i <= m; i++) result += log((double)(i - 1 - n)) - log((double)(i));
    }
    else if (m <= n)
    {
      for (INT i = 1; i <= -1 - n; i++) result += log((double)(n - m + i)) - log((double)(i));
    }
    else
      return result;
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$T_{ab}\f$ arma::mat.
 *
 *  For a definition of the \f$T_{ab}\f$ arma::mat, cf. \ref sec02_01.
 *  Special case of trz12() in hfbcs.f.
 */

void Basis::calcTab(void)
{
  DBG_ENTER;

  if (!Tab.empty())
  {
    DBG_LEAVE;
  }

  Qnumbers QNn_z("n_z", n_zGlobalMax + 2);
  Tab = arma::zeros(QNn_z.nb, QNn_z.nb);

  if (dMax == 1)
  {
    Tab.eye();
    DBG_LEAVE;
  }

  UINT   n_zMax1 = QNn_z.nb;
  UINT   n_zMax2 = QNn_z.nb;
  double b1      = b_z;
  double b2      = b_z;
  double a1      = 1 / pow(b1, 2);
  double a2      = 1 / pow(b2, 2);
  Tab(0, 0) = sqrt(2 * sqrt(a1 * a2) / (a1 + a2)) * exp(-0.5 * pow(d_0, 2) * (a1 * a2) / (a1 + a2));

  for (INT n = 1; n < n_zMax2; n++)
  {
    Tab(0, n) = d_0 * a1 * sqrt(2 * a2) * Tab(0, n - 1) / ((a1 + a2) * sqrt(n));

    if (n > 1)
    {
      Tab(0, n) += (a2 - a1) * sqrt(n - 1) * Tab(0, n - 2) / ((a1 + a2) * sqrt(n));
    }
  }

  for (INT n1 = 1; n1 < n_zMax1; n1++)
  {
    for (INT n2 = 0; n2 < n_zMax2; n2++)
    {
      Tab(n1, n2) = -1 * d_0 * a2 * sqrt(2 * a1) * Tab(n1 - 1, n2) / ((a1 + a2) * sqrt(n1));

      if (n2 > 0)
      {
        Tab(n1, n2) += 2 * sqrt(a1 * a2 * double(n2)) * Tab(n1 - 1, n2 - 1) / ((a1 + a2) * sqrt(n1));
      }

      if (n1 > 1)
      {
        Tab(n1, n2) += (a1 - a2) * sqrt(n1 - 1) * Tab(n1 - 2, n2) / ((a1 + a2) * sqrt(n1));
      }
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the transformation matrix to another Basis instance.
 */

arma::mat Basis::getTransformationMatrix(Basis &otherBasis)
{
  DBG_ENTER;

  Basis &b2 = otherBasis;

  // identical bases ?

  bool identical = true;

  if (mMax != otherBasis.mMax) identical = false;

  for (INT m = 0; m < mMax; m++)
  {
    if (omegaIndexHO(m).n_rows != otherBasis.omegaIndexHO(m).n_rows) identical = false;
  }

  if (fabs(b_r - otherBasis.b_r) > 1e-4) identical = false;

  if (fabs(b_z - otherBasis.b_z) > 1e-4) identical = false;

  if (fabs(d_0 - otherBasis.d_0) > 1e-4) identical = false;

  if (identical)
  {
    DBG_RETURN(arma::eye(HOqn.nb, HOqn.nb));
  }

  // bases are different, let's compute the transformation

  b2.calcTab();
  Multi<arma::mat> oMat          = Basis::getOverlapMatrixHO(b2);

  // for (auto &m : oMat)
  // {
  //   INFO(Tools::matPrint(m.second));
  // }

  arma::mat           overlapOHR = getOverlapMatrixHOR(b2);
  arma::mat           overlapOHZ = getOverlapMatrixHOZ(b2);
  // re-ordering overlapOHZ (Berger2ct style)
  UVEC nzdTodnz1(overlapOHZ.n_rows);
  UVEC nzdTodnz2(overlapOHZ.n_cols);
  INT        i = 0;

  for (INT d = 0; d < dMax; d++)
    for (UINT n_z = 0; n_z < overlapOHZ.n_rows / dMax; n_z++) nzdTodnz1(i++) = n_z * dMax + d;

  i = 0;

  for (INT d = 0; d < b2.dMax; d++)
    for (UINT n_z = 0; n_z < overlapOHZ.n_cols / b2.dMax; n_z++) nzdTodnz2(i++) = n_z * b2.dMax + d;

  arma::mat overlapOHZb = overlapOHZ(nzdTodnz1, nzdTodnz2).t();
  arma::mat transMat    = arma::zeros(b2.HOqn.nb, HOqn.nb);
  UINT      nbBlocks1   = HOqn.calcBlocks({0, 4});        // get blocks of same (m, s) QN
  UINT      nbBlocks2   = b2.HOqn.calcBlocks({0, 4});     // get blocks of same (m, s) QN

  for (UINT msb1 = 0; msb1 < nbBlocks1; msb1++)
  {
    Qnumbers &qnb1 = HOqn.blocks[msb1];
    INT       m1   = qnb1(0, 0);
    INT       s1   = qnb1(4, 0);

    if (m1 - s1 < 0) continue;     // no omega < 0 blocks

    for (UINT msb2 = 0; msb2 < nbBlocks2; msb2++)
    {
      Qnumbers &qnb2 = b2.HOqn.blocks[msb2];
      INT       m2   = qnb2(0, 0);
      INT       s2   = qnb2(4, 0);

      if (m2 - s2 < 0) continue;     // no omega < 0 blocks

      if ((m1 != m2) || (s1 != s2)) continue;     // same (m, s)

      if ((dMax == 2) && (b2.dMax == 2))
      {
        for (INT n2 = 0; n2 < b2.nMax(m2); n2++)
        {
          INT imn2    = b2.QNmn.find({m2, n2});
          auto n_zMax2 = b2.n_zMax(m2, n2);

          // 1. solve
          //
          // |   1    Tz22 |
          // |      \      | E = V E
          // | Tz22    1   |
          arma::mat szt  = arma::eye(n_zMax2 * 2, n_zMax2 * 2);
          arma::mat tz12 = b2.Tab;
          szt.submat(0, n_zMax2, n_zMax2 - 1, n_zMax2 * 2 - 1) =
            tz12.submat(0, 0, n_zMax2 - 1, n_zMax2 - 1);
          szt = arma::symmatu(szt);

          arma::mat E;
          arma::vec V;

          if (!Tools::checkSymmetry(szt, "szt in " + std::string(__PRETTY_FUNCTION__)))
            szt = arma::symmatu(szt);

          Tools::eig_sym(V, E, szt);

          UVEC sortIndex = arma::sort_index(V, "descend");
          V                    = V.elem(sortIndex);
          E                    = E.cols(sortIndex);

          // 2. construct F = { E.col(c) / v(c) if v(c) > 5E-4 }
          arma::mat F  = arma::zeros(arma::size(E));
          UINT      nb = 0;

          for (UINT i = 0; i < E.n_cols; i++)
          {
            if (V(i) > 5e-4) F.col(nb++) = E.col(i) / V(i);
          }

          F = F.head_cols(nb);
          E = E.head_cols(nb);
          // 3. construct A = E * F^T
          arma::mat A = E * F.t();

          // 4. construct Tz^21 = A^T * Tz^21
          auto nr = A.n_rows;
          auto nc = overlapOHZb.n_cols;

          arma::mat Tz21 = arma::zeros(nr, nc);

          auto middler  = overlapOHZb.n_rows / 2;
          auto middler2 = nr / 2;

          Tz21.submat(0, 0, middler2 - 1, nc - 1) = overlapOHZb.submat(0, 0, middler2 - 1, nc - 1);
          Tz21.submat(middler2, 0, middler2 * 2 - 1, nc - 1) =
            overlapOHZb.submat(middler, 0, middler + middler2 - 1, nc - 1);

          Tz21 = A.t() * Tz21;

          for (INT n1 = 0; n1 < nMax(m1); n1++)
          {
            INT imn1 = QNmn.find({m1, n1});

            for (INT n_z1 = 0; n_z1 < n_zMax(m1, n1); n_z1++)
            {
              for (INT n_z2 = 0; n_z2 < b2.n_zMax(m2, n2); n_z2++)
              {
                for (INT d1 = 0; d1 < dMax; d1++)
                {
                  INT i1 = HOqn.find({m1, n1, n_z1, d1, s1});

                  for (INT d2 = 0; d2 < b2.dMax; d2++)
                  {
                    INT i2 = b2.HOqn.find({m2, n2, n_z2, d2, s2});

                    if ((i1 < 0) || (i2 < 0)) continue;

                    transMat(i2, i1) =
                      overlapOHR(imn1, imn2) * Tz21(n_z2 + nr / 2 * d2, n_z1 + nc / 2 * d1);
                  } // d2
                } // d1
              } // n_z2
            } // n_z1
          } // n2
        } // n1
      }
      else
      {
        transMat.submat(qnb2.filter, qnb1.filter) = oMat(m1, s1).t();
      }
    }     // msb2
  }       // msb1

  // transMat generated ! \o/ Thanks Marc and Jean-François !
  DBG_RETURN(transMat);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the corrected overlap matrix between HO1/2ct states.
 */

Multi<arma::mat> Basis::getOverlapMatrixHO(Basis &otherBasis)
{
  DBG_ENTER;

  Multi<arma::mat>  result;
  arma::mat         overlapOHR = getOverlapMatrixHOR(otherBasis);
  arma::mat         overlapOHZ = getOverlapMatrixHOZ(otherBasis);
  Multi<UVEC> msBra      = calcMSBlocks((*this));
  Multi<UVEC> msKet      = calcMSBlocks(otherBasis);

  for (INT m = 0; m < MIN(mMax, otherBasis.mMax); m++)
  {
    for (INT s = 0; s < 2; s++)
    {
      if (m - s < 0) continue;

      // if s > m, omega < 0 so we don't care
      if (s > m) continue;

      result(m, s) = arma::zeros(msBra(m, s).n_elem, msKet(m, s).n_elem);

      for (UINT iB = 0; iB < msBra(m, s).n_elem; iB++)
      {
        for (UINT iK = 0; iK < msKet(m, s).n_elem; iK++)
        {
          result(m, s)(iB, iK) =
            getOverlapHO(otherBasis, msBra(m, s)(iB), msKet(m, s)(iK), overlapOHR, overlapOHZ);
        }
      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get the basis overlap matrix for \f$z\f$ parts.
 *
 *  Cf. trz12() in hfbcs.f.
 */

arma::mat Basis::getOverlapMatrixHOZ(Basis &otherBasis)
{
  DBG_ENTER;

  INT       n_zMax1 = MAX(n_zGlobalMax, otherBasis.n_zGlobalMax);
  INT       n_zMax2 = MAX(n_zGlobalMax, otherBasis.n_zGlobalMax);
  arma::mat tabqq   = arma::zeros(dMax * n_zMax1, otherBasis.dMax * n_zMax2);

  double b_z1 = b_z;
  double b_z2 = otherBasis.b_z;
  double a1   = 1.0 / pow(b_z1, 2);
  double a2   = 1.0 / pow(b_z2, 2);

  for (INT d1 = 0; d1 < dMax; d1++)
  {
    double zd1 = (0.5 - double(d1)) * d_0;

    for (INT d2 = 0; d2 < otherBasis.dMax; d2++)
    {
      double zd2 = (0.5 - double(d2)) * otherBasis.d_0;
      double d0  = zd1 - zd2;
      tabqq(0 * dMax + d1, 0 * otherBasis.dMax + d2) =
        sqrt(2 * sqrt(a1 * a2) / (a1 + a2)) * exp(-0.5 * pow(d0, 2) * (a1 * a2) / (a1 + a2));

      for (INT n = 1; n < n_zMax2; n++)
      {
        tabqq(0 * dMax + d1, n * otherBasis.dMax + d2) =
          d0 * a1 * sqrt(2 * a2) * tabqq(0 * dMax + d1, (n - 1) * otherBasis.dMax + d2) /
          ((a1 + a2) * sqrt(n));

        if (n > 1)
        {
          tabqq(0 * dMax + d1, n * otherBasis.dMax + d2) +=
            (a2 - a1) * sqrt(n - 1) * tabqq(0 * dMax + d1, (n - 2) * otherBasis.dMax + d2) /
            ((a1 + a2) * sqrt(n));
        }
      }

      for (INT n1 = 1; n1 < n_zMax1; n1++)
      {
        for (INT n2 = 0; n2 < n_zMax2; n2++)
        {
          tabqq(n1 * dMax + d1, n2 * otherBasis.dMax + d2) =
            -1 * d0 * a2 * sqrt(2 * a1) * tabqq((n1 - 1) * dMax + d1, n2 * otherBasis.dMax + d2) /
            ((a1 + a2) * sqrt(n1));

          if (n2 > 0)
          {
            tabqq(n1 * dMax + d1, n2 * otherBasis.dMax + d2) +=
              2 * sqrt(a1 * a2 * double(n2)) *
              tabqq((n1 - 1) * dMax + d1, (n2 - 1) * otherBasis.dMax + d2) /
              ((a1 + a2) * sqrt(n1));
          }

          if (n1 > 1)
          {
            tabqq(n1 * dMax + d1, n2 * otherBasis.dMax + d2) +=
              (a1 - a2) * sqrt(n1 - 1) * tabqq((n1 - 2) * dMax + d1, n2 * otherBasis.dMax + d2) /
              ((a1 + a2) * sqrt(n1));
          }
        }
      }
    }
  }

  DBG_RETURN(tabqq);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Get the basis overlap matrix for \f$r_\perp\f$ parts.
 *
 *  Cf. zchro() in hfbcs.f.
 */

arma::mat Basis::getOverlapMatrixHOR(Basis &otherBasis)
{
  DBG_ENTER;

  arma::mat rabqq = arma::zeros(QNmn.nb, otherBasis.QNmn.nb);
  double    bet1  = 1.0 / (b_r * b_r);
  double    bet2  = 1.0 / (otherBasis.b_r * otherBasis.b_r);
  double    a     = 2.0 * sqrt(bet1 * bet2);
  double    blam  = a / (bet1 + bet2);
  a               = (bet1 - bet2) / a;
  double blmu     = blam * a;
  double bmu2     = a * a;
  double bla2     = blam * blam;

  for (INT m = 0; m < MIN(mMax, otherBasis.mMax); m++)
  {
    arma::mat tx = arma::zeros(MAX(nMax(m), otherBasis.nMax(m)), MAX(nMax(m), otherBasis.nMax(m)));
    a            = pow(blam, m + 1);

    for (INT n1 = 0; n1 < MIN(nMax(m), otherBasis.nMax(m)); n1++)
    {
      INT    i  = n1 + m;
      double aw = a;

      for (INT n2 = n1; n2 < MAX(nMax(m), otherBasis.nMax(m)); n2++)
      {
        double b = 1.0;
        INT    j = n2 - n1;

        if (n1 > 0)
        {
          double bw = 1.0;

          for (INT k = 1; k <= n1; k++)
          {
            bw = -(bw * bmu2 / double(k * (j + k))) * double((n1 + 1 - k) * (i + 1 - k));
            b += bw;
          }
        }

        tx(n1, n2) = b * aw;

        aw *= (blmu / double(j + 1)) * sqrt((n2 + 1.0) * (n2 + m + 1.0));
      }

      a *= bla2;
    }

    a = 1.0;

    for (auto n1 = 0; n1 < nMax(m); n1++)
    {
      double b = a;
      UINT   i = MIN(MIN(nMax(m), otherBasis.nMax(m)), n1);

      if (i > 0)
      {
        for (INT n2 = 0; n2 < i; n2++)
        {
          UINT ia       = QNmn.find({m, n1});
          UINT ib       = otherBasis.QNmn.find({m, n2});
          rabqq(ia, ib) = tx(n2, n1) * b;
          b             = -b;
        }
      }

      if (n1 < otherBasis.nMax(m))
      {
        for (auto n2 = n1; n2 < otherBasis.nMax(m); n2++)
        {
          UINT ia       = QNmn.find({m, n1});
          UINT ib       = otherBasis.QNmn.find({m, n2});
          rabqq(ia, ib) = tx(n1, n2);
        }
      }

      a = -a;
    }

    /* for comparison with berger2ct zchro().
        for (INT n1 = 0; n1 < nMax(m); n1++)
        {
          for (INT n2 = 0; n2 < otherBasis.nMax(m); n2++)
          {
            INT ia = QNmn.find({m, n1});
            INT ib = otherBasis.QNmn.find({m, n2});
            printf("%2d %2d %2d %15.12f\n", m+1, n1+1, n2+1, rabqq(ia, ib));
          }
        }
    */
  }

  DBG_RETURN(rabqq);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the overlap between two HO states.
 */

double Basis::getOverlapHO(Basis &otherBasis, UINT iB, UINT iK, const arma::mat &overlapOHR,
                           const arma::mat &overlapOHZ)
{
  //  DBG_ENTER;

  if (HOqn(4, iB) != otherBasis.HOqn(4, iK))
  {
    //    DBG_LEAVE;
    return 0.;
  }

  INT mB   = HOqn(0, iB);
  INT nB   = HOqn(1, iB);
  INT n_zB = HOqn(2, iB);
  INT dB   = HOqn(3, iB);
  INT mK   = otherBasis.HOqn(0, iK);
  INT nK   = otherBasis.HOqn(1, iK);
  INT n_zK = otherBasis.HOqn(2, iK);
  INT dK   = otherBasis.HOqn(3, iK);
  UINT iRr = 0;

  for (UINT blR = 0; blR < mB; blR++)
  {
    iRr += nMax(blR);
  }

  iRr += nB;
  UINT iZr = n_zB * dMax + dB;
  UINT iRc = 0;

  for (UINT blR = 0; blR < mK; blR++)
  {
    iRc += otherBasis.nMax(blR);
  }

  iRc += nK;
  UINT   iZc     = n_zK * otherBasis.dMax + dK;
  double overlap = overlapOHR(iRr, iRc) * overlapOHZ(iZr, iZc);

  //  DBG_LEAVE;
  return overlap;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** PrINT a matrix using Berger's order.
 */

void Basis::debugMat(const arma::mat &mat) const
{
  DBG_ENTER;

  UVEC index = arma::zeros<UVEC>(HOqn.nb);

  for (INT om = 0; om < mMax; om++)
  {
    INT i = 0;

    for (INT s = 0; s < sMax; s++)
    {
      INT m = om + s;

      if (m < 0) continue;

      if (m >= mMax) continue;

      for (INT n = 0; n < nMax(m); n++)
      {
        for (INT d = 0; d < dMax; d++)
        {
          for (INT n_z = 0; n_z < n_zMax(m, n); n_z++)
          {
            index(i) = HOqn.find({m, n, n_z, d, s});
            i++;
          }
        }
      }
    }

    Tools::info(PF("om: %2d", om), mat.submat(index.head(i), index.head(i)), true);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a hermite polynomial for a set of values.
 *
 * \param n_z The \f$n_z\f$ quantity.
 * \param zeta The \f$\zeta\f$ values.
 *
 * The hermite polynomials are calculated using the recursive property
 *
 *  \f{eqnarray*}{
 *  H_0(\zeta)&=&1\\
 *  H_1(\zeta)&=&2\zeta\\
 *  H_{n}(\zeta)&=&2\zeta H_{n-1}(\zeta)-2(n-1)H_{n-2}(\zeta)
 *  \f}
 */

arma::vec Basis::hermite(INT n_z, const arma::vec &zeta)
{
  if (n_z == 0) return arma::ones(zeta.n_elem);

  if (n_z == 1) return 2 * zeta;

#ifdef USE_POLY_CACHING
  static std::map<std::pair<int, std::size_t>, arma::vec>    hermiteBuffer;
  std::pair<int, std::size_t>                                tempPair(n_z, Tools::vec2hash(zeta));
  std::map<std::pair<int, std::size_t>, arma::vec>::iterator it = hermiteBuffer.find(tempPair);

  if (it != hermiteBuffer.end())
  {
    return (*it).second;     // it->second ?
  }

#endif
  arma::vec val = 2 * zeta % hermite(n_z - 1, zeta) - 2 * double(n_z - 1) * hermite(n_z - 2, zeta);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<int, std::size_t>, arma::vec> tempPairPair(tempPair, val);
  hermiteBuffer.insert(tempPairPair);
#endif

  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a generalized laguerre polynomial for a set of values.
 *
 * \param m The \f$m\f$ quantity.
 * \param n The \f$n\f$ quantity.
 * \param eta The \f$\eta\f$ quantities.
 *
 * The generalized laguerre polynomials are calculated using the recursive property
 *
 * \f{eqnarray*}{
 * L^{m}_0(\eta)&=&1\\
 * L^{m}_1(\eta)&=&1+m-\eta\\
 * L^{m}_{n}(\eta)&=&\left(2+\frac{m-1-\eta}{n}\right)L^{m}_{n-1}(\eta)-\left(1+\frac{m-1}{n}\right)L^{m}_{n-2}(\eta)
 * \f}
 */

arma::vec Basis::laguerre(INT m, INT n, const arma::vec &eta)
{
  if (n == 0) return arma::ones(eta.n_elem);

  if (n == 1) return double(1 + m) - eta;

#ifdef USE_POLY_CACHING
  static std::map<std::pair<std::pair<int, int>, std::size_t>, arma::vec> laguerreBuffer;
  std::pair<int, int>                                                     tempPair(m, n);
  std::pair<std::pair<int, int>, std::size_t> tempPairPair(tempPair, Tools::vec2hash(eta));
  std::map<std::pair<std::pair<int, int>, std::size_t>, arma::vec>::iterator it =
    laguerreBuffer.find(tempPairPair);

  if (it != laguerreBuffer.end())
  {
    return it->second;
  }

#endif
  arma::vec val = (2 + (double(m - 1) - eta) / (double)n) % laguerre(m, n - 1, eta) -
                  (1 + double(m - 1) / (double)n) * laguerre(m, n - 2, eta);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<std::pair<int, int>, std::size_t>, arma::vec> tempPairPairPair(
    tempPairPair, val);
  laguerreBuffer.insert(tempPairPairPair);
#endif
  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a legendre polynomial for a set of values.
 *
 * \param l The \f$l\f$ quantity.
 * \param x The \f$x\f$ values.
 *
 * The legendre polynomials are calculated using the recursive property
 *
 *  \f{eqnarray*}
 *  P_0(x)&=&1\\
 *  P_1(x)&=&x\\
 *  P_{l}(x)&=&\frac{1}{l} \left( (2l-1) x P_{l-1}(x)-(l-1)P_{n-2}(x) \right)\f}
 */

arma::vec Basis::legendre(INT l, const arma::vec &x)
{
  if (l == 0) return arma::ones(x.n_elem);

  if (l == 1) return x;

#ifdef USE_POLY_CACHING
  static std::map<std::pair<int, std::size_t>, arma::vec>    legendreBuffer;
  std::pair<int, std::size_t>                                tempPair(l, Tools::vec2hash(x));
  std::map<std::pair<int, std::size_t>, arma::vec>::iterator it = legendreBuffer.find(tempPair);

  if (it != legendreBuffer.end())
  {
    return it->second;
  }

#endif
  arma::vec val =
    (double(2 * l - 1) * x % legendre(l - 1, x) - double(l - 1) * legendre(l - 2, x)) / double(l);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<int, std::size_t>, arma::vec> tempPairPair(tempPair, val);
  legendreBuffer.insert(tempPairPair);
#endif
  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a hermite polynomial value.
 *
 * \param n_z The \f$n_z\f$ quantity.
 * \param zeta The \f$\zeta\f$ quantity.
 *
 * The hermite polynomials are calculated using the recursive property
 *
 *  \f{eqnarray*}
 *  H_0(\zeta)&=&1\\
 *  H_1(\zeta)&=&2\zeta\\
 *  H_{n}(\zeta)&=&2\zeta H_{n-1}(\zeta)-2(n-1)H_{n-2}(\zeta)\f}
 */

double Basis::hermite(INT n_z, double zeta)
{
  if (n_z == 0) return 1;

  if (n_z == 1) return 2 * zeta;

#ifdef USE_POLY_CACHING
  static std::map<std::pair<int, double>, double>    hermiteBuffer2;
  std::pair<int, double>                             tempPair(n_z, zeta);
  std::map<std::pair<int, double>, double>::iterator it =
    hermiteBuffer2.find(tempPair);     // already calculated ?

  if (it != hermiteBuffer2.end()) return it->second;     // yes !

#endif
  double val = 2 * zeta * hermite(n_z - 1, zeta) - 2 * double(n_z - 1) * hermite(n_z - 2, zeta);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<int, double>, double> tempPairPair(tempPair, val);
  hermiteBuffer2.insert(tempPairPair);     // store it
#endif
  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a generalized laguerre polynomial value.
 *
 * \param m The \f$m\f$ quantity.
 * \param n The \f$n\f$ quantity.
 * \param eta The \f$\eta\f$ quantity.
 *
 * The generalized laguerre polynomials are calculated using the recursive property
 *
 * \f{eqnarray*}
 * L^{m}_0(\eta)&=&1\\
 * L^{m}_1(\eta)&=&1+m-\eta\\
 * L^{m}_{n}(\eta)&=&\left(2+\frac{m-1-\eta}{n}\right)L^{m}_{n-1}(\eta)-\left(1+\frac{m-1}{n}\right)L^{m}_{n-2}(\eta)
 * \f}
 */

double Basis::laguerre(INT m, INT n, double eta)
{
  if (n == 0)
  {
    return 1.0;
  }

  if (n == 1)
  {
    return double(1 + m) - eta;
  }

#ifdef USE_POLY_CACHING
  static std::map<std::pair<std::pair<int, int>, double>, double>    laguerreBuffer2;
  std::pair<int, int>                                                tempPair(m, n);
  std::pair<std::pair<int, int>, double>                             tempPairPair(tempPair, eta);
  std::map<std::pair<std::pair<int, int>, double>, double>::iterator it =
    laguerreBuffer2.find(tempPairPair);

  if (it != laguerreBuffer2.end())
  {
    return it->second;
  }

#endif
  double val = (2.0 + (double(m - 1) - eta) / double(n)) * laguerre(m, n - 1, eta) -
               (1.0 + double(m - 1) / double(n)) * laguerre(m, n - 2, eta);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<std::pair<int, int>, double>, double> tempPairPairPair(tempPairPair,
      val);
  laguerreBuffer2.insert(tempPairPairPair);
#endif
  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Evaluate a legendre polynomial value.
 *
 * \param l The \f$l\f$ quantity.
 * \param x The \f$x\f$ values.
 *
 * The legendre polynomials are calculated using the recursive property
 *
 *  \f{eqnarray*}
 *  P_0(x)&=&1\\
 *  P_1(x)&=&x\\
 *  P_{l}(x)&=&\frac{1}{l} \left( (2l-1) x P_{l-1}(x)-(l-1)P_{n-2}(x) \right)\f}
 *
 */

double Basis::legendre(INT l, double x)
{
  if (l == 0) return 1;

  if (l == 1) return x;

#ifdef USE_POLY_CACHING
  static std::map<std::pair<int, double>, double>    legendreBuffer2;
  std::pair<int, double>                             tempPair(l, x);
  std::map<std::pair<int, double>, double>::iterator it =
    legendreBuffer2.find(tempPair);     // already calculated ?

  if (it != legendreBuffer2.end()) return it->second;     // yes !

#endif
  double val =
    (double(2 * l - 1) * x * legendre(l - 1, x) - double(l - 1) * legendre(l - 2, x)) / double(l);
#ifdef USE_POLY_CACHING
  const std::pair<std::pair<int, double>, double> tempPairPair(tempPair, val);
  legendreBuffer2.insert(tempPairPair);     // store it
#endif
  return val;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Talman coefficients for the \f$r\f$-part.
 *
*/

void Basis::calcTalmanr(void)
{
  DBG_ENTER;

  if (!talmanr.empty())
  {
    DBG_LEAVE;
  }

  INT M = mMax;
  Nmaxr = 0;

  for (INT ma = 0 ; ma < M + 1 ; ma++)
  {
    INT Ma = MAX(ma - 1, 0);

    for (INT mb = 0 ; mb < M + 1 ; mb++)
    {
      INT Mb = MAX(mb - 1, 0);

      for (INT na = 0 ; na < nMax(Ma) + 1 ; na++) //ici
      {
        for (INT nb = 0 ; nb < nMax(Mb) + 1 ; nb++) //ici
        {
          INT maxLocal = (2 * (na + nb) + abs(ma) + abs(mb) - abs(ma - mb)) / 2; // TODO: check validity ?

          if (maxLocal > Nmaxr)
          {
            Nmaxr = maxLocal;
          }
        }
      }
    }
  }

  Nmaxr += 1; //For it to be a strict max.

  for (INT ma = 0; ma < M + 1; ma++) //ici
  {
    INT Ma = MAX(ma - 1, 0); //SO purposes for nMax definition.

    for (INT na = 0; na < nMax(Ma) + 1; na++) //ici
    {
      INT xa = 2 * na + abs(ma);

      for (INT mb = -M ; mb < M + 1; mb++) //ici
      {
        INT Mb = mb;//SO purposes for nMax definition.

        if (mb >= 0)
        {
          Mb = MAX(mb - 1, 0);
        }
        if (mb < 0)
        {
          Mb = mb + 1;
        }

        if (mb < 0)
        {
          Mb = mb + 1;
        }

        INT mc = mb - ma;

        for (INT nb = 0; nb < nMax(abs(Mb)) + 1; nb++) //ici
        {
          INT xb = 2 * nb + abs(mb);

          talmanr(ma, na, mb, nb) = arma::zeros(Nmaxr);
          talmanr(-ma, na, -mb, nb) = arma::zeros(Nmaxr);

          arma::vec &talvec = talmanr(ma, na, mb, nb);
          arma::vec &talvec_min = talmanr(-ma, na, -mb, nb);

          for (INT nc = 0; nc < Nmaxr; nc++)
          {
            double fac = sqrt(fact[na] * fact[na + abs(ma)]
                              * fact[nb] * fact[nb + abs(mb)]
                              * fact[nc] * fact[nc + abs(mc)]);

            if ((na + nb - nc) % 2 != 0) fac *= -1.0;

            INT xc = 2 * nc + abs(mc);
            INT xpa = (xc + xb - xa) / 2;
            INT xpb = (xc + xa - xb) / 2;
            INT xpc = (xa + xb - xc) / 2;
            INT mdMin = -1 * MIN(MIN(xpa + ma, xpb + mb), xpc + ma + mb);
            INT mdMax = MIN(MIN(xpa - ma, xpb - mb), xpc - (ma + mb)) + 1;
            double sum = 0.0;

            for (INT m = mdMin; m < mdMax; m += 2)
            {
              double prod0 = fact[(xpa + (m + ma)) / 2];
              double prod1 = fact[(xpa - (m + ma)) / 2];
              double prod2 = fact[(xpb + (m + mb)) / 2];
              double prod3 = fact[(xpb - (m + mb)) / 2];
              double prod4 = fact[(xpc + (m + ma + mb)) / 2];
              double prod5 = fact[(xpc - (m + ma + mb)) / 2];
              sum += 1.0 / (prod0 * prod1 * prod2 * prod3 * prod4 * prod5);
            }

            talvec(nc) = fac * sum;
            talvec_min(nc) = talvec(nc);
          } // nc
        } // nb
      } // mb
    } // na
  } // ma

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Moshinsky coefficients for the \f$r\f$-part.
 *  Those are the reduced Moshinsky, that is to say :
 *  M^(0,0)(0,na+nb+|ma|)_(ma,na,-ma,nb).
 *  It is the only one needed for what follows.
 */

void Basis::calcMoshinskyr(void)
{
  DBG_ENTER;

  // dependencies
  calcTalmanr(); // TODO: remove this dependency

  if (!moshinskyr.empty())
  {
    DBG_LEAVE;
  }

  calcTalmanr();

  INT M = mMax;
  INT nl = 0;

  for (INT ma = M * -2 + 2 ; ma < 2 * M - 1; ma++)
  {
    for (INT nb = 0; nb < Nmaxr; nb++)
    {
      moshinskyr(ma, nb) = arma::zeros(Nmaxr);

      for (INT na = 0; na < Nmaxr; na++)
      {
        nl = na + nb + abs(ma); // only non-zero components
        moshinskyr(ma, nb)(na) = getFullMoshinskyr(ma, nb, -ma, na, 0, 0, 0, nl);
      } // na
    } // nb
  } // ma

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Moshinsky coefficients for the \f$z\f$-part.
 *
 * Appendix D, PhD J.-F. Berger.
 * J.-P. Ebran's derivations, p. 70, Eq. (IV-B.3).
 */

void Basis::calcMoshinskyz(void)
{
  DBG_ENTER;
  // dependencies
  calcTab();

  if (!moshinskyz.empty())
  {
    DBG_LEAVE;
  }

  moshinskyz = arma::zeros(n_zGlobalMax * 2 + 2, n_zGlobalMax * 2 + 2);

  INT n_zd;
  double factor;

  for (INT n_za = 0; n_za < n_zGlobalMax * 2 + 1  ; n_za++)
  {
    for (INT n_zb = 0; n_zb < n_zGlobalMax * 2 + 1 ; n_zb++)
    {
      n_zd = n_za + n_zb;
      factor = pow(2.0, n_zd / -2.0) / (mnz[n_za] * mnz[n_zb] * mnz[n_zd]);
      moshinskyz(n_za, n_zb) = factor * pow(mnz[n_za] * mnz[n_zb], 2) * pow(-1.0, n_zb);
    } // n_zb
  } // n_za

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the full overlap matrix between HO1/2ct states.
 */


arma::mat Basis::getFullOverlap(Basis &otherBasis)
{
  DBG_ENTER;

  arma::mat result = arma::zeros(HOqn.nb, otherBasis.HOqn.nb);

  arma::mat overlapOHR = getOverlapMatrixHOR(otherBasis);
  arma::mat overlapOHZ = getOverlapMatrixHOZ(otherBasis);
  auto msBra = calcMSBlocks((*this));
  auto msKet = calcMSBlocks(otherBasis);

  for (INT m = 0; m < MIN(mMax, otherBasis.mMax); m++)
    for (INT s = 0; s < sMax; s++)
    {
      // if s > m, omega < 0 so we don't care
      if (s > m) continue;

      for (UINT iB = 0; iB < msBra(m, s).n_elem; iB++)
        for (UINT iK = 0; iK < msKet(m, s).n_elem; iK++)
        {
          UINT ia = msBra(m, s)(iB);
          UINT ib = msKet(m, s)(iK);
          result(ia, ib) = getOverlapHO(
                             otherBasis,
                             msBra(m, s)(iB),
                             msKet(m, s)(iK),
                             overlapOHR, overlapOHZ
                           );
        }
    }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculates 4 HO-z overlap <z_alpha,z_beta|z_gamma,z_delta>. Used for testing only.
 */

void Basis::calcIz(void)
{
  DBG_ENTER;

  //dependencies
  calcTalmanz();
  calcMoshinskyz();

  //usefull quantities.
  INT DMax = (dMax - 1) * 3 + 1;
  double fac_bz = 1 / (sqrt(2 * b_z * sqrt(PI)));
  double fac_nabla = 1 / (b_z * sqrt(2));

  //Merges Moshisnky-z and z wave function.
  Multi<arma::mat> addMosz; //addMosz(d_alpha+2*d_gamma,d_delta+2*d_beta)(n_za,n_zb)

  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    double dv_delta = (double(d_delta) - 0.5) * d_0;

    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      double dv_beta = (double(d_beta) - 0.5) * d_0;

      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        double dv_alpha = (double(d_alpha) - 0.5) * d_0;

        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          double dv_gamma = (double(d_gamma) - 0.5) * d_0;
          double k_ag = 0.5 * (dv_alpha + dv_gamma);
          double k_bd = 0.5 * (dv_beta  + dv_delta);
          double d_agbd = (k_ag - k_bd) / sqrt(2);
          double fac_exp = exp(-0.5 * pow(d_agbd / b_z, 2));
          arma::mat &addmosz = addMosz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta);
          addmosz = arma::zeros(n_zGlobalMax * 2 + 2, n_zGlobalMax * 2 + 2);

          for (INT n_za = 0 ; n_za < n_zGlobalMax * 2 + 2 ; n_za++)
          {
            for (INT n_zb = 0 ; n_zb < n_zGlobalMax * 2 + 2 ; n_zb++)
            {
              addmosz(n_za, n_zb) = fac_bz * fac_exp * moshinskyz(n_za, n_zb) * zPartScalar(d_agbd, n_za + n_zb);
            }//n_zb
          }//n_za
        }//d_gamma
      }//d_alpha
    }//d_beta
  }//d_delta


  //Merges addMoshisnky-z and talmanz.
  //Multi<arma::vec,4> Jz;
  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          arma::mat &addmosz = addMosz(d_alpha + 2 * d_gamma, d_delta + 2 * d_beta);

          for (INT n_zbeta = 0 ; n_zbeta < n_zGlobalMax + 1 ; n_zbeta++)
          {
            for (INT n_zdelta = 0 ; n_zdelta < n_zGlobalMax + 1 ; n_zdelta++)
            {
              Jz(n_zdelta, n_zbeta, d_delta + 2 * d_beta, d_alpha + 2 * d_gamma) = addmosz * talmanz(n_zdelta, d_delta, n_zbeta, d_beta);
            }//n_zdelta
          }//n_zbeta
        }//d_gamma
      }//d_alpha
    }//d_beta
  }//d_delta


  //Iz init.
  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT n_zalpha = 0 ; n_zalpha < n_zGlobalMax + 1 ; n_zalpha++)
      {
        for (INT n_zgamma = 0 ; n_zgamma < n_zGlobalMax + 1 ; n_zgamma++)
        {
          Iz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma) = arma::zeros(n_zGlobalMax + 1, n_zGlobalMax + 1, DMax);
        }//n_zgamma
      }//n_zalpha
    }//d_gamma
  }//d_alpha

  //NablaIz init.
  for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
  {
    for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
    {
      for (INT n_zbeta = 0 ; n_zbeta < n_zGlobalMax ; n_zbeta++)
      {
        for (INT n_zdelta = 0 ; n_zdelta < n_zGlobalMax ; n_zdelta++)
        {
          NablaIz(n_zdelta, n_zbeta, d_alpha + 2 * d_gamma) = arma::zeros(n_zGlobalMax, n_zGlobalMax, DMax);
        }//n_zdelta
      }//n_zbeta
    }//d_gamma
  }//d_alpha

  //Merges J-z and talmanz.
  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          for (INT n_zbeta = 0 ; n_zbeta < n_zGlobalMax + 1 ; n_zbeta++)
          {
            for (INT n_zdelta = 0 ; n_zdelta < n_zGlobalMax + 1 ; n_zdelta++)
            {
              arma::vec &jvec = Jz(n_zdelta, n_zbeta, d_delta + 2 * d_beta, d_alpha + 2 * d_gamma);

              for (INT n_zalpha = 0 ; n_zalpha < n_zGlobalMax + 1 ; n_zalpha++)
              {
                for (INT n_zgamma = 0 ; n_zgamma < n_zGlobalMax + 1 ; n_zgamma++)
                {
                  Iz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = arma::accu(talmanz(n_zalpha, d_alpha, n_zgamma, d_gamma) % jvec);
                }//n_zgamma
              }//n_zalpha
            }//n_zdelta
          }//n_zbeta
        }//d_gamma
      }//d_alpha
    }//d_beta
  }//d_delta

  //calc NablaIz (derivatives on z_beta).
  for (INT d_delta = 0 ; d_delta < dMax ; d_delta++)
  {
    for (INT d_beta = 0 ; d_beta < dMax ; d_beta++)
    {
      for (INT d_alpha = 0 ; d_alpha < dMax ; d_alpha++)
      {
        for (INT d_gamma = 0 ; d_gamma < dMax ; d_gamma++)
        {
          for (INT n_zalpha = 0 ; n_zalpha < n_zGlobalMax ; n_zalpha++)
          {
            for (INT n_zgamma = 0 ; n_zgamma < n_zGlobalMax ; n_zgamma++)
            {
              for (INT n_zdelta = 0 ; n_zdelta < n_zGlobalMax ; n_zdelta++)
              {
                //if n_zbeta = 0;
                NablaIz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, 0, d_delta + 2 * d_beta) = - fac_nabla * Iz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, 1, d_delta + 2 * d_beta);

                for (INT n_zbeta = 1 ; n_zbeta < n_zGlobalMax ; n_zbeta++)
                {
                  NablaIz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, n_zbeta, d_delta + 2 * d_beta) = fac_nabla * (sqrt(n_zbeta) * Iz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, n_zbeta - 1, d_delta + 2 * d_beta)
                      - sqrt(n_zbeta + 1) * Iz(n_zalpha, n_zgamma, d_alpha + 2 * d_gamma)(n_zdelta, n_zbeta + 1, d_delta + 2 * d_beta));
                }//n_zbeta
              }//n_zdelta
            }//n_zgamma
          }//n_zalpha
        }//d_gamma
      }//d_alpha
    }//d_beta
  }//d_delta


  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate 4 HO-r overlap <r_alpha,r_beta|r_gamma,r_delta>. Used for testing only.
 */

void Basis::calcIr(void)
{
  DBG_ENTER;

  //dependencies
  calcTalmanr();
  calcMoshinskyr();

  //usefull quantities.
  double fac_br = 1 / (2 * std::pow(b_r, 2) * PI);
  double fac_nabla  = 1 / (b_r * sqrt(2));
  double fac_nabla2 = fac_nabla * fac_nabla;

  INT maxn = Nmaxr;
  arma::mat mosR = arma::zeros(maxn, maxn);

  for (INT n_a = 0 ; n_a < maxn ; n_a++)
  {
    for (INT n_b = 0 ; n_b < maxn ; n_b++)
    {
      mosR(n_a, n_b) = moshinskyr(0, n_a)(n_b);
    }//n_b
  }//n_a

  //Calculates Jr.
  Multi<arma::vec> Jr;

  //If m = 0.
  for (INT m = 0 ; m < mMax + 1 ; m++)
  {
    INT M = MAX(m - 1, 0);

    for (INT n_alpha = 0 ; n_alpha < nMax(M) + 1 ; n_alpha++)
    {
      for (INT n_gamma = 0 ; n_gamma < nMax(M) + 1 ; n_gamma++)
      {
        Jr(m, n_alpha, n_gamma) = fac_br * mosR * talmanr(m, n_alpha, m, n_gamma);
      }//n_gamma
    }//n_alpha
  }//m

  //calculates Ir.
  for (INT m = 0 ; m < mMax + 1 ; m++)
  {
    INT M = MAX(m - 1, 0);

    for (INT mp = 0 ; mp < mMax + 1 ; mp++)
    {
      INT Mp = MAX(mp - 1, 0);

      for (INT n_alpha = 0 ; n_alpha < nMax(M) + 1 ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < nMax(M) + 1 ; n_gamma++)
        {
          Ir(m, n_alpha, n_gamma, mp) = arma::zeros(nMax(Mp) + 1, nMax(Mp) + 1);

          for (INT n_delta = 0 ; n_delta < nMax(Mp) + 1 ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < nMax(Mp) + 1 ; n_beta++)
            {
              Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta) = arma::accu(Jr(m, n_alpha, n_gamma) % talmanr(mp, n_delta, mp, n_beta));
            }//n_beta
          }//n_delta
        }//n_gamma
      }//n_alpha
    }//mp
  }//m

  //Calculates Ir_plus (derivatives on the last indice).
  for (INT m = 0 ; m < mMax ; m++)
  {
    for (INT mp = 1 ; mp < mMax ; mp++)
    {
      for (INT n_alpha = 0 ; n_alpha < nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < nMax(m) ; n_gamma++)
        {
          Ir_plus(m, n_alpha, n_gamma, mp) = arma::zeros(nMax(mp), nMax(mp - 1));

          for (INT n_delta = 0 ; n_delta < nMax(mp) ; n_delta++)
          {
            //If n_beta = 0.
            Ir_plus(m, n_alpha, n_gamma, mp)(n_delta, 0) = fac_nabla * sqrt(mp) * Ir(m, n_alpha, n_gamma, mp)(n_delta, 0);

            for (INT n_beta = 1 ; n_beta < nMax(mp - 1) ; n_beta++)
            {
              Ir_plus(m, n_alpha, n_gamma, mp)(n_delta, n_beta) = fac_nabla *
                (sqrt(mp + n_beta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta) +
                 sqrt(n_beta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta - 1));
            }//n_beta
          }//n_delta
        }//n_gamma
      }//n_alpha
    }//mp
  }//m

  //Calculates Ir_minus (derivatives on the last indice).
  for (INT m = 0 ; m < mMax ; m++)
  {
    for (INT mp = 0 ; mp < mMax - 1 ; mp++)
    {
      for (INT n_alpha = 0 ; n_alpha < nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < nMax(m) ; n_gamma++)
        {
          Ir_minus(m, n_alpha, n_gamma, mp) = arma::zeros(nMax(mp), nMax(mp + 1));

          for (INT n_delta = 0 ; n_delta < nMax(mp) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < nMax(mp + 1) ; n_beta++)
            {
              Ir_minus(m, n_alpha, n_gamma, mp)(n_delta, n_beta) = -fac_nabla *
                (sqrt(mp + 1 + n_beta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta) +
                 sqrt(n_beta + 1) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta + 1));
            }//n_beta
          }//n_delta
        }//n_gamma
      }//n_alpha
    }//mp
  }//m

  //Calculates Ir_mm (derivative + on star beta, derivative - on no-star delta).
  for (INT m = 0 ; m < mMax ; m++)
  {
    for (INT mp = 0 ; mp < mMax - 1 ; mp++)
    {
      for (INT n_alpha = 0 ; n_alpha < nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < nMax(m) ; n_gamma++)
        {
          Ir_mm(m, n_alpha, n_gamma, mp) = arma::zeros(nMax(mp + 1), nMax(mp + 1));

          for (INT n_delta = 0 ; n_delta < nMax(mp + 1) ; n_delta++)
          {
            for (INT n_beta = 0 ; n_beta < nMax(mp + 1) ; n_beta++)
            {
              Ir_mm(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + 1 + n_beta) * sqrt(mp + 1 + n_delta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta);
              Ir_mm(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + 1 + n_beta) * sqrt(n_delta + 1) * Ir(m, n_alpha, n_gamma, mp)(n_delta + 1, n_beta);
              Ir_mm(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + 1 + n_delta) * sqrt(n_beta + 1) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta + 1);
              Ir_mm(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(n_delta + 1) * sqrt(n_beta + 1) * Ir(m, n_alpha, n_gamma, mp)(n_delta + 1, n_beta + 1);
            }//n_beta
          }//n_delta
        }//n_gamma
      }//n_alpha
    }//mp
  }//m

  //Calculates Ir_pp (derivative - on star beta, derivative + on no-star delta).
  for (INT m = 0 ; m < mMax ; m++)
  {
    for (INT mp = 2 ; mp < mMax + 1 ; mp++)
    {
      for (INT n_alpha = 0 ; n_alpha < nMax(m) ; n_alpha++)
      {
        for (INT n_gamma = 0 ; n_gamma < nMax(m) ; n_gamma++)
        {
          Ir_pp(m, n_alpha, n_gamma, mp) = arma::zeros(nMax(mp - 1), nMax(mp - 1));
          //if n_beta = 0 et n_delta=0.
          Ir_pp(m, n_alpha, n_gamma, mp)(0, 0) += -fac_nabla2 * sqrt(mp) * sqrt(mp) * Ir(m, n_alpha, n_gamma, mp)(0, 0);

          //if n_beta = 0 et n_delta>0.
          for (INT n_delta = 1 ; n_delta < nMax(mp - 1) ; n_delta++)
          {
            Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, 0) += -fac_nabla2 * sqrt(mp) * sqrt(mp + n_delta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, 0);
            Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, 0) += -fac_nabla2 * sqrt(mp) * sqrt(n_delta) * Ir(m, n_alpha, n_gamma, mp)(n_delta - 1, 0);
          }//n_delta

          //if n_beta > 0 et n_delta =0.
          for (INT n_beta = 1 ; n_beta < nMax(mp - 1) ; n_beta++)
          {
            Ir_pp(m, n_alpha, n_gamma, mp)(0, n_beta) += -fac_nabla2 * sqrt(mp + n_beta) * sqrt(mp) * Ir(m, n_alpha, n_gamma, mp)(0, n_beta);
            Ir_pp(m, n_alpha, n_gamma, mp)(0, n_beta) += -fac_nabla2 * sqrt(mp) * sqrt(n_beta) * Ir(m, n_alpha, n_gamma, mp)(0, n_beta - 1);
          }//n_beta

          for (INT n_delta = 1 ; n_delta < nMax(mp - 1) ; n_delta++)
          {
            for (INT n_beta = 1 ; n_beta < nMax(mp - 1) ; n_beta++)
            {
              Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + n_beta) * sqrt(mp + n_delta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta);
              Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + n_beta) * sqrt(n_delta) * Ir(m, n_alpha, n_gamma, mp)(n_delta - 1, n_beta);
              Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(mp + n_delta) * sqrt(n_beta) * Ir(m, n_alpha, n_gamma, mp)(n_delta, n_beta - 1);
              Ir_pp(m, n_alpha, n_gamma, mp)(n_delta, n_beta) += -fac_nabla2 * sqrt(n_delta) * sqrt(n_beta) * Ir(m, n_alpha, n_gamma, mp)(n_delta - 1, n_beta - 1);
            }//n_beta
          }//n_delta
        }//n_gamma
      }//n_alpha
    }//mp
  }//m

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Clebsch-Gordan coefficients \f$\lbrace j_1 m_1 j_2 m_2 | j m\rbrace\f$.
 */

double Basis::CGcoeff(INT j1, INT m1, INT j2, INT m2, INT j, INT m)
{

  if (m != m1 + m2) return 0.0;

  INT jmin = abs(j1 - j2);
  INT jmax = j1 + j2;

  if ((j < jmin)||(jmax < j)) return 0.0;

  double s0 = (double(2 * j + 1) * std::tgamma(double(j + j1 - j2 + 1)) * std::tgamma(double(j - j1 + j2 + 1)) * std::tgamma(double(j1 + j2 - j + 1))) / std::tgamma(double(j + j1 + j2 + 2));
  s0 = sqrt(s0);

  double s = std::tgamma(double(j + m + 1)) * std::tgamma(double(j - m + 1));
  double s1 = std::tgamma(double(j1 + m1 + 1)) * std::tgamma(double(j1 - m1 + 1));
  double s2 = std::tgamma(double(j2 + m2 + 1)) * std::tgamma(double(j2 - m2 + 1));
  s = sqrt(s * s1 * s2);

  INT kMax = MIN(MIN(j1 + j2 - j, j1 - m1), j2 + m2);

  double CG = 0.0;

  for (INT k = 0; k <= kMax; k++)
  {
    double k1 = std::tgamma(double(j1 + j2 - j - k + 1));
    double k2 = std::tgamma(double(j1 - m1 - k + 1));
    double k3 = std::tgamma(double(j2 + m2 - k + 1));
    double k4 = std::tgamma(double(j - j2 + m1 + k + 1));
    double k5 = std::tgamma(double(j - j1 - m2 + k + 1));

    double temp = pow(-1, double(k)) / (std::tgamma(double(k + 1)) * k1 * k2 * k3 * k4 * k5);

    CG += temp;
  }

  return s0 * s * CG;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Moshinsky coefficients for the \f$r\f$-part. Those are the reduced Moshinsky, that is to say : M^(0,0)(0,na+nb+|ma|)_(ma,na,-ma,nb). It is the only one needed for what follows.
 */

// Moshinsky coefficient used for the spin-orbit interaction
void Basis::calcMoshinskyDG(void)
{
  DBG_ENTER;

  if (!moshinskyDG.empty()) DBG_LEAVE;

  //Dependencies
  calcTalmanr(); // TODO: remove this dependency

  INT M = mMax;
  NmaxDG = 0;
  for (INT alpha = -1; alpha < 2; alpha++)
  {
    for (INT m_nu = - 2 * M; m_nu < 2 * M + 1; m_nu++)
    {
      INT m_mu = - m_nu - alpha;
      for (INT n_nu = 0; n_nu < Nmaxr; n_nu++)
      {
        moshinskyDG(alpha, m_nu, n_nu) = arma::zeros(Nmaxr);

        for (INT n_mu = 0; n_mu < Nmaxr; n_mu++)
        {
          INT n_sigma = n_mu + n_nu + (abs(m_mu) + abs(m_nu) - abs(alpha)) / 2; // TODO: check validity ?

          if(NmaxDG < n_sigma)
          {
            NmaxDG = n_sigma;
          }
          moshinskyDG(alpha, m_nu, n_nu)(n_mu) = getFullMoshinskyr(m_mu, n_mu, m_nu, n_nu, 0, 0, -alpha, n_sigma);
        } //n_mu
      } //n_nu
    } //m_nu
  } //alpha

  NmaxDG += 1; // to get a strict maximum

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the Moshinsky coefficients for the \f$r\f$-part. Those are the reduced Moshinsky, that is to say : M^(0,0)(0,na+nb+|ma|)_(ma,na,-ma,nb). It is the only one needed for what follows.
 */

// Moshinsky coefficient used for the tensor interaction
void Basis::calcMoshinskyTDG(void)
{
  DBG_ENTER;

  if (!moshinskyTDG.empty()) DBG_LEAVE;

  //Dependencies
  calcTalmanr(); // TODO: remove this dependency

  INT M = mMax;
  NmaxTDG = 0; // upper bound of n_sigma
  for (INT k = -1; k < 3; k++)
  {
    for (INT m_nu = - M; m_nu < 2 * M + 1; m_nu++)
    {
      INT m_mu = - m_nu + k;
      for (INT n_nu = 0; n_nu < Nmaxr; n_nu++)
      {
        moshinskyTDG(k, m_nu, n_nu) = arma::zeros(Nmaxr);

        for (INT n_mu = 0; n_mu < Nmaxr; n_mu++)
        {
          INT n_sigma = n_mu + n_nu + (abs(m_mu) + abs(m_nu) - abs(k)) / 2; // TODO: check validity ?

          if (NmaxTDG < n_sigma)
          {
            NmaxTDG = n_sigma;
          }
          moshinskyTDG(k, m_nu, n_nu)(n_mu) = getFullMoshinskyr(m_mu, n_mu, m_nu, n_nu, 0, 0, k, n_sigma);
        } //n_nu
      } //n_mu
    } //m_mu
  } //k

  NmaxTDG += 1; // to get a strict maximum

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice string representation of the basis parameters.
 */

const std::string Basis::getNiceInfo(void)
{
  DBG_ENTER;

  std::string result;

  if (d_0 > 0.0) result += Tools::treeStr({
                                            {"br", PF("%.3f", b_r)},
                                            {"bz", PF("%.3f", b_z)},
                                            {"d0", PF("%.3f", d_0)},
                                          }, true);
  else           result += Tools::treeStr({
                                            {"br", PF("%.3f", b_r)},
                                            {"bz", PF("%.3f", b_z)},
                                          }, true);

  DBG_RETURN(result);
}
