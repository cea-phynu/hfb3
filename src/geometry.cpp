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

#include "geometry.h"
#include "tools.h"
#include "discrete.h"
#include "multipole_operators.h"

/** \file
 *  \brief Methods of the Geometry class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** The dummy constructor.
 */

Geometry::Geometry(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 * \param _mesh a reference to a State object.
 */

Geometry::Geometry(State &state) :
  mesh(Mesh::regular(0, 0, -20, 20, 0, 20, 401, 1, 401))
{
  DBG_ENTER;

  if (!state.checkSolution())
  {
    DBG_LEAVE;
  }

  Discrete discrete(&state.basis, mesh);
  dens = discrete.getLocalXZ(state.rho(NEUTRON) + state.rho(PROTON ), true);

  (*this) = Geometry(mesh, dens, state.sys);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor.
 *
 * \param _mesh A Mesh object, defining the integration domain.
 * \param _dens The local density to be integrated.
 * \param _izNeck z-Index of the neck, if a neck has been found.
 */

Geometry::Geometry(Mesh _mesh, arma::mat _dens, System _system, INT _izNeck) : system(_system), izNeck(_izNeck), mesh(_mesh), dens(_dens)
{
  DBG_ENTER;

  neckPos  = -999.0;
  neckDens = -999.0;
  intLeft  = -999.0;
  intRight = -999.0;
  intTotal = -999.0;
  q_neck   = 999.0;

  calcAll();

  if (izNeck == -1)
  {
#ifdef USE_NECK_SCHUNCK
    calc_neck_abcissa_schunck();
#else
    calc_neck_abcissa();
#endif
  }

  // TODO make a more generic formula to find z = 0 index
  if ((izNeck < 0) || (izNeck > mesh.az.nb - 1)) izNeck = -1;

  calc_fragment_properties();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the neck abcissa.
 *
 * zNeck is determined such that the integral of C - rho(z) in the neck area is
 * splitted in two equal parts. The constant reads
 * C = rhoMin + 0.25 * (rhoMax - rhoMin),
 * with: rhoMin = minimum of rho(z) in the neck
 *       rhoMax = maximum of rho(z).
 */

void Geometry::calc_neck_abcissa(void)
{
  DBG_ENTER;
  INT ixZero = -1;

  for (INT ix = 0; ix < mesh.ax.nb; ix++)
  {
    if (fabs(mesh.ax.p(ix)) < EPSILON)
    {
      ixZero = ix;
      break;
    }
  }

  ASSERT(ixZero != -1, "Mesh does no contain z = 0.0");

  // Compute the density integrated over x
  arma::mat z = Tools::matFromRow(mesh.az.p.t(), mesh.ax.nb);
  arma::mat r = Tools::matFromCol(mesh.ax.p, mesh.az.nb);
  arma::mat wx = Tools::matFromCol(mesh.ax.we, mesh.az.nb);
  arma::mat velem = r * 2.0 * PI;
  arma::mat func = dens % wx % velem;
  arma::rowvec zDens = arma::sum(func);

  // Maximum of the density integrated over x
  double dMax = -1.0;

  for (INT iz = 0; iz < mesh.az.nb; iz++)
  {
    double d = zDens(iz);

    if (d > dMax) dMax = d;
  }

  // INFO("maximum density: %e", dMax);
  ASSERT(dMax >= 0.0, "Negative maximum density");

  // Left half-density point.
  INT izHalfDensL = -1;

  for (INT iz = 0; iz < mesh.az.nb; iz++)
  {
    double d = zDens(iz);

    if (d > dMax * 0.1)
    {
      izHalfDensL = iz;
      break;
    }
  }

  // INFO("left partial-density: %d", izHalfDensL);
  ASSERT(izHalfDensL >= 0, "Did not find the left partial-density point");

  // Right half-density point.
  INT izHalfDensR = -1;

  for (INT iz = INT(mesh.az.nb) - 1; iz >= 0; iz--)
  {
    double d = zDens(iz);

    if (d > dMax * 0.1)
    {
      izHalfDensR = iz;
      break;
    }
  }

  // INFO("right partial-density: %d", izHalfDensR);
  ASSERT(izHalfDensR >= 0, "Did not find the right partial-density point");

  // "Middle" of the density.
  INT izMiddle = (izHalfDensR + izHalfDensL) / 2;
  // Left maximum density.
  INT izMaxL = -1;
  double dMaxL = 0.0;

  for (INT iz = 0; iz <= izMiddle; iz++)
  {
    double d = zDens(iz);

    if (d > dMaxL)
    {
      dMaxL = d;
      izMaxL = iz;
    }
  }

  // INFO("left maximum density: %e", dMaxL);
  ASSERT(izMaxL >= 0, "Did not find the left maximum density");

  // Right maximum density.
  INT izMaxR = -1;
  double dMaxR = 0.0;

  for (INT iz = INT(mesh.az.nb) - 1; iz >= izMiddle; iz--)
  {
    double d = zDens(iz);

    if (d > dMaxR)
    {
      dMaxR = d;
      izMaxR = iz;
    }
  }

  // INFO("right maximum density: %e", dMaxR);
  ASSERT(izMaxR >= 0, "Did not find the right maximum density");

  // Minimum density between the two peaks
  double dMin = 9999.0;
  INT izMin = -1;

  for (INT iz = izMaxL; iz <= izMaxR; iz++)
  {
    double d = zDens(iz);

    if (d < dMin)
    {
      dMin = d;
      izMin = iz;
    }
  }

  // INFO("center minimum density: %e %d", dMin, izMin);
  ASSERT(izMin >= 0, "Did not find the minimum density between peaks");

  izNeck = izMin;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the fragment properties.
 */

void Geometry::calc_fragment_properties(void)
{
  DBG_ENTER;
  INT ixZero = -1;

  for (INT ix = 0; ix < INT(mesh.ax.nb); ix++)
  {
    if (fabs(mesh.ax.p(ix)) < EPSILON)
    {
      ixZero = ix;
      break;
    }
  }

  if (ixZero == -1)
  {
    Tools::warning("mesh does no contain x = 0.0 in Geometry::calc_neck_abcissa()");
    DBG_LEAVE;
    return;
  }

  if (izNeck == -1)
  {
    // no neck has been found
    neckPos  = 999.0;
    neckDens = 999.0;
    intLeft  = 999.0;
    intRight = 999.0;
    intTotal = integrate(0, INT_DENS_T);
    q_neck   = 999.0;
  }
  else
  {
    neckPos  = mesh.az.p(izNeck);
    neckDens = dens(ixZero, izNeck);
    intLeft  = integrate(izNeck, INT_DENS_L);
    intRight = integrate(izNeck, INT_DENS_R);
    intTotal = integrate(izNeck, INT_DENS_T);
    q_neck   = integrate(izNeck, INT_Q_NECK);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the neck abcissa (N. Schunck's method).
 *
 * The neck's abcissa minimizes <qNeck>.
 */

void Geometry::calc_neck_abcissa_schunck(void)
{
  DBG_ENTER;
  INT ixZero = -1;

  for (INT ix = 0; ix < mesh.ax.nb; ix++)
  {
    if (fabs(mesh.ax.p(ix)) < EPSILON)
    {
      ixZero = ix;
      break;
    }
  }

  if (ixZero == -1)
  {
    std::stringstream ss;
    ss << "mesh does no contain z = 0.0 in Geometry::calc_neck_abcissa_schunck()";
    Tools::warning(ss.str());
    DBG_LEAVE;
    return;
  }

  // Compute the density integrated over x
  arma::mat r = Tools::matFromCol(mesh.ax.p, mesh.az.nb);
  arma::mat wx = Tools::matFromCol(mesh.ax.we, mesh.az.nb);
  arma::mat wz = Tools::matFromRow(mesh.az.we.t(), mesh.az.nb);
  arma::mat velem = r * 2.0 * PI;
  arma::mat func = dens % wx % velem;
  double aNeck = 1.0;
  arma::colvec zDens = arma::sum(func).t();
  // Compute the qNeck for any zNeck
  arma::mat neckFunc(mesh.az.nb, mesh.az.nb);

  for (INT izn = 0; izn < mesh.az.nb; izn++)
  {
    for (INT iz = izn; iz < mesh.az.nb; iz++)
    {
      double zDiff = mesh.az.p(iz) - mesh.az.p(izn);
      double expVal = exp( - 1.0 * zDiff * zDiff / aNeck );
      neckFunc(izn, iz) = expVal;
      neckFunc(iz, izn) = expVal;
    }
  }

  arma::colvec qNeckOfZ = neckFunc % wz * zDens;
  // Maximum of the density integrated over x
  double dMax = -1.0;

  for (INT iz = 0; iz < mesh.az.nb; iz++)
  {
    double d = zDens(iz);

    if (d > dMax) dMax = d;
  }

  if (dMax < 0.0)
  {
    std::stringstream ss;
    ss << "negative maximum density in Geometry::calc_neck_abcissa_schunck()";
    Tools::warning(ss.str());
    DBG_LEAVE;
    return;
  }

  // Left half-density point.
  INT izHalfDensL = -1;

  for (INT iz = 0; iz < mesh.az.nb; iz++)
  {
    double d = zDens(iz);

    if (d > dMax / 2)
    {
      izHalfDensL = iz;
      break;
    }
  }

  if (izHalfDensL < 0)
  {
    std::stringstream ss;
    ss << "did not find the left half-density point in Geometry::calc_neck_abcissa_schunck()";
    Tools::warning(ss.str());
    DBG_LEAVE;
    return;
  }

  // Right half-density point.
  INT izHalfDensR = -1;

  for (INT iz = INT(mesh.az.nb) - 1; iz >= 0; iz--)
  {
    double d = zDens(iz);

    if (d > dMax / 2)
    {
      izHalfDensR = iz;
      break;
    }
  }

  if (izHalfDensR < 0)
  {
    std::stringstream ss;
    ss << "did not find the right half-density point in Geometry::calc_neck_abcissa_schunck()";
    Tools::warning(ss.str());
    DBG_LEAVE;
    return;
  }

  // Search for qNeck local and global minima in the interval [izHalfDensL, izHalfDensR]
  INT neckCount = 0;
  double qNeckMin = 9999.0;

  for (INT iz = izHalfDensL + 1; iz <= izHalfDensR - 1; iz++)
  {
    bool isLocalMin = qNeckOfZ(iz) - qNeckOfZ(iz - 1) < 0. && qNeckOfZ(iz + 1) - qNeckOfZ(iz) > 0.;

    if (isLocalMin)
    {
      neckCount++;

      if (qNeckOfZ(iz) < qNeckMin)
      {
        qNeckMin = qNeckOfZ(iz);
        izNeck = iz;
      }
    }
  }

  // Choose neck value
  if (neckCount == 0)
  {
    izNeck = (izHalfDensR + izHalfDensL) / 2;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the multipole moments.
 */

void Geometry::calcQlm(void)
{
  DBG_ENTER;

  qlm     = arma::zeros(7);
  qlm(0)  = integrate(0, INT_Q00);
  qlm(1)  = integrate(0, INT_Q10);
  qlm(2)  = integrate(0, INT_Q20);
  qlm(3)  = integrate(0, INT_Q30);
  qlm(4)  = integrate(0, INT_Q40);
  qlm(5)  = integrate(0, INT_Q50);
  qlm(6)  = integrate(0, INT_Q60);

  calcBeta();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate all quantities.
 */

void Geometry::calcAll(void)
{
  DBG_ENTER;

  calcQlm();

  State state;

  intTotal = integrate(0, INT_DENS_T);
  int_r   = integrate(0, INT_R);
  int_r2  = integrate(0, INT_R2);
  int_rp2 = integrate(0, INT_RP2);
  rms     = sqrt(int_r2 / intTotal);
  radius  = int_r / intTotal;

  /*
  * J. L. Friar and J. W. Negele
  * Adv. Nucl. Phys. 8 219 (1975)
  */
  chargeRMS = sqrt(int_r2 / intTotal + PROTON_RADIUS_SQ + DARWINFOLDY + NEUTRON_RADIUS_SQ * (double)(system.nNeut) / (double)(system.nProt));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

double Geometry::getOldChargeRMS(double Z, double N)
{
  DBG_ENTER;

  double hbarOmega = 1.85 + 35.5 * pow(Z + N, -1.0 / 3.0);
  double b = 41.47 / hbarOmega / (double)(Z + N);
  double gb = 0.7144;

  int_r2   = integrate(0, INT_R2);
  intTotal = integrate(0, INT_DENS_T);

  double chRMS = sqrt(int_r2 / intTotal + 1.5 * (gb * gb - b) - 0.1161 * N / Z);

  DBG_RETURN(chRMS);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the cylindrical integral of the given density.
 *
 * \param izCut The integral boundary.
 * \param what The integral type.
 */

double Geometry::integrate(INT izCut, INT what)
{
  arma::mat z = Tools::matFromRow(mesh.az.p.t(), mesh.ax.nb);
  arma::mat z2 = z % z;
  arma::mat z3 = z2 % z;
  arma::mat z4 = z3 % z;
  arma::mat z5 = z4 % z;
  arma::mat z6 = z5 % z;
  arma::mat rp = Tools::matFromCol(mesh.ax.p, mesh.az.nb);
  arma::mat rp2 = rp % rp;
  arma::mat rp4 = rp2 % rp2;
  arma::mat rp6 = rp4 % rp2;
  arma::mat r = arma::sqrt(rp2 + z2);
  arma::mat wz = Tools::matFromRow(mesh.az.we.t(), mesh.ax.nb);
  arma::mat wx = Tools::matFromCol(mesh.ax.we, mesh.az.nb);
  arma::mat ones = arma::ones(mesh.ax.nb, mesh.az.nb);
  arma::mat velem = rp * 2.0 * PI;
  double aNeck = 1.0; // [fm]
  arma::mat vneck = arma::exp(-1.0 * ((z - neckPos) / aNeck) % ((z - neckPos) / aNeck));
  arma::mat func;
  double result = 0.0;

  switch (what)
  {
    case INT_DENS_T:
      func = dens % wx % wz % velem;
      result = arma::accu(func);
      break;

    case INT_DENS_L:
      func = dens % wx % wz % velem;
      result = arma::accu(func.cols(0, izCut));
      break;

    case INT_DENS_R:
      func = dens % wx % wz % velem;
      result = arma::accu(func.cols(izCut + 1, mesh.az.nb - 1));
      break;

    case INT_Q00:
      func = 1. / 2.  * sqrt( 1. / PI) * dens % wx % wz % velem;
      // func = dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q10:
      func = 1. / 2.  * sqrt( 3. / PI) * z % dens % wx % wz % velem;
      func = z % dens % wx % wz % velem;
      func = z % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q20:
      if      (general.compatibility == General::COMPAT_NONE   ) func = 1. / 4.  * sqrt( 5. / PI) * (2 * z2 - rp2) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_BERGER ) func = (2 * z2 - rp2) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_ROBLEDO) func = (    z2 - 1. / 2. * rp2) % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q30:
      if      (general.compatibility == General::COMPAT_NONE   ) func = 1. / 4.  * sqrt( 7. / PI) * (2 * z3 - 3 * z % rp2) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_BERGER ) func = (z3 - 3.0 / 2.0 * z % rp2) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_ROBLEDO) func = (z3 - 3.0 / 2.0 * z % rp2) % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q40:
      if      (general.compatibility == General::COMPAT_NONE   ) func = 1. / 16. * sqrt( 9. / PI) * (8 * z4 - 24 * z2 % rp2 + 3 * rp4) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_BERGER ) func = (z4 - 24.0 / 8.0 * z2 % rp2 + 3.0 / 8.0 * rp4) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_ROBLEDO) func = (z4 - 24.0 / 8.0 * z2 % rp2 + 3.0 / 8.0 * rp4) % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q50:
      if      (general.compatibility == General::COMPAT_NONE   ) func = 1. / 16. * sqrt(11. / PI) * (8 * z5 - 40 * z3 % rp2 + 15 * z % rp4) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_BERGER ) func = (z5 - 40.0 / 8.0 * z3 % rp2 + 15.0 / 8.0 * z % rp4) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_ROBLEDO) func = (z5 - 40.0 / 8.0 * z3 % rp2 + 15.0 / 8.0 * z % rp4) % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_Q60:
      if      (general.compatibility == General::COMPAT_NONE   ) func = 1. / 32. * sqrt(13. / PI) * (16 * z6 - 120 * z4 % rp2 + 90 * z2 % rp4 - 5 * rp6) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_BERGER ) func = (z6 - 120.0 / 16.0 * z4 % rp2 + 90.0 / 16.0 * z2 % rp4 - 5.0 / 16.0 * rp6) % dens % wx % wz % velem;
      else if (general.compatibility == General::COMPAT_ROBLEDO) func = (z6 - 120.0 / 16.0 * z4 % rp2 + 90.0 / 16.0 * z2 % rp4 - 5.0 / 16.0 * rp6) % dens % wx % wz % velem;

      result = arma::accu(func);
      break;

    case INT_RP2:
      func = rp2 % dens % wx % wz % velem;
      result = arma::accu(func);
      break;

    case INT_R:
      func = r % dens % wx % wz % velem;
      result = arma::accu(func);
      break;

    case INT_R2:
      func = r % r % dens % wx % wz % velem;
      result = arma::accu(func);
      break;


    case INT_Q_NECK:
      func = vneck % dens % wx % wz % velem;
      result = arma::accu(func);
      break;

    default:
      std::stringstream ss;
      ss << "Unknown 2nd argument in Geometry::integrate(INT, INT): what = " << what;
      Tools::warning(ss.str());
      DBG_RETURN(0.0);
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the multipole expansion of the density.
 */

arma::cube Geometry::getDensityMultipoleExpansion(State &state, const arma::vec &rVals)
{
  DBG_ENTER;

  INT quadOrder = 40;

  Axis quadAxis(Axis::GAUSS_LEGENDRE, quadOrder);

  arma::vec p = quadAxis.p;
  arma::vec w = quadAxis.w;

  INT lambdaMax = 9;

  // shperical harmonics
  arma::mat ylm = arma::zeros(quadOrder, lambdaMax);
  arma::mat mp = arma::zeros(quadOrder, lambdaMax);

  for (INT lambda = 0; lambda < lambdaMax; lambda++)
    mp.col(lambda) = arma::pow(-p, lambda);

  ylm.col(0) = 1. / 2.  * sqrt( 1. / PI) * (mp.col(0));
  ylm.col(1) = 1. / 2.  * sqrt( 3. / PI) * (     mp.col(1)     );
  ylm.col(2) = 1. / 4.  * sqrt( 5. / PI) * (3. * mp.col(2) - 1.);
  ylm.col(3) = 1. / 4.  * sqrt( 7. / PI) * (5. * mp.col(3) - 3. * mp.col(1));
  ylm.col(4) = 1. / 16. * sqrt( 9. / PI) * (35. * mp.col(4) - 30. * mp.col(2) + 3.);
  ylm.col(5) = 1. / 16. * sqrt(11. / PI) * (63. * mp.col(5) - 70. * mp.col(3) + 15. * mp.col(1));
  ylm.col(6) = 1. / 32. * sqrt(13. / PI) * (231. * mp.col(6) - 315. * mp.col(4) + 105. * mp.col(2) - 5.);
  ylm.col(7) = 1. / 32. * sqrt(15. / PI) * (429. * mp.col(7) - 693. * mp.col(5) + 315. * mp.col(3) - 35. * mp.col(1));
  ylm.col(8) = 1. / 256.* sqrt(17. / PI) * (6435. * mp.col(8) - 12012. * mp.col(6) + 6930. * mp.col(4) - 1260. * mp.col(2) + 35.);

  ASSERT(p.n_elem == quadOrder, "Wrong number of nodes");
  ASSERT(w.n_elem == quadOrder, "Wrong number of weights");

  arma::cube result = arma::zeros(rVals.n_elem, lambdaMax, 2);
  INT ir = 0;

  for (auto r: rVals)
  {
    arma::vec rPerp = r * arma::sqrt(1 - p % p);
    arma::vec z     = -r * p;

    arma::mat dens = arma::zeros(quadOrder, 2);
    for (INT i = 0; i < quadOrder; i++)
    {
      arma::vec res = Discrete::getLocalDensity(state, rPerp(i), 0.0, z(i));
      dens.row(i) = res.t();
    }


    for (INT lambda = 0; lambda < lambdaMax; lambda+=2)
    {
      result(ir, lambda, NEUTRON) = 2.0 * PI * arma::accu(dens.col(NEUTRON) % w % ylm.col(lambda));
      result(ir, lambda, PROTON ) = 2.0 * PI * arma::accu(dens.col(PROTON ) % w % ylm.col(lambda));
    }

    ir++;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the \f$\beta\f$ values.
 */

void Geometry::calcBeta(void)
{
  DBG_ENTER;

  beta = 0.0;

  double Z = system.nProt;
  double N = system.nNeut;
  double A = Z + N;

  beta = MultipoleOperators::getBetaFromQ20(A, A, qlm(2));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Geometry::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  std::vector<std::pair<std::string, std::string>> list;

  if (isShort)
  {
    if (qlm.empty())
    {
      DBG_RETURN(PF_RED("empty"));
    }

    for (UINT l = 0; l < 4; l++)
    {
      list.push_back({PF("Q%1d0", l), PF("%.6e", qlm(l))});
    }

    list.push_back({"beta", PF("%.6e", beta)});
  }
  else
  {
    list.push_back({"Geometry", ""});

    if (qlm.empty())
    {
      list.push_back({"qlm", PF_RED("empty")});
    }
    else
    {
      for (UINT l = 0; l < qlm.n_rows; l++)
      {
        list.push_back({PF("<Q%1d0> ", l),
                        "(" + PF_YELLOW("total") + ":" + PF_BLUE("%13.6e", qlm(l)) + ")"});
      }

      list.push_back({"beta  ", "(" + PF_YELLOW("total") +
                      ":" + PF_BLUE("%13.6e", beta) + ")"});
    }
  }

  result += Tools::treeStr(list, isShort);

  DBG_RETURN(result);
}
