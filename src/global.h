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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <map>
#include <set>
#include <zlib.h>
#include <armadillo>
#include <cstdarg>

/** \file
 *  \brief Global definitions.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/// Disable bound check in arma::* objects.
// #define ARMA_NO_DEBUG

//=====================
//===== Typedefs ======
//=====================

/// Positive integer (mainly used for indexing).
typedef long long unsigned int UINT;

/// Integer (mainly used for quantum number values).
typedef int INT;

/// Vector container for real numbers.
typedef arma::Col<double> VEC;

/// Vector container for positive integers.
typedef arma::Col<UINT> UVEC;

/// Vector container for integers.
typedef arma::Col<INT> IVEC;

/// Matrix container for real numbers.
typedef arma::Mat<double> MAT;

/// Matrix container for positive integers.
typedef arma::Mat<UINT> UMAT;

/// Matrix container for integers.
typedef arma::Mat<INT> IMAT;

/// Cube container for real numbers.
typedef arma::Cube<double> CUBE;

/// Cube container for positive integers.
typedef arma::Cube<UINT> UCUBE;

/// Cube container for integers.
typedef arma::Cube<INT> ICUBE;

//=====================
//===== Constants =====
//=====================

/// \f$\pi\f$ numerical constant.
#define PI 3.141592653589793238462643

/** \f$\hbar.c\f$ numerical constant.
 *
 * Using 2019 SI base units:
 *
 * h         [J.s] := 6.62607015e-34
 * c         [m/s] := 2.99792458e+08
 * e         [A.s] := 1.602176634e-19
 * hbarc    [eV.m]  = h[J.s]/(2pi).c[m/s]/e[A.s]
 *                  = 1.9732698045930246e-07
 * hbarc  [MeV.fm]  = hbarc [eV.m] * 1e-06 * 1e+15
 *                  = 1.9732698045930246e+02
 * alpha        [] := 7.2973525664e-03
 */
#define HBARC 1.9732698045930246e+02

/// Default version string (should be set at compile time).
#ifndef CFG_GIT_VERSION
#define CFG_GIT_VERSION "undefined"
#endif

/// \f$\alpha\f$ numerical constant (fine-structure constant).
#define ALPHA 7.2973525664e-03

/// \f$(\hbar.c)^2\f$ numerical constant.
#define HBARC2 HBARC*HBARC

/// The mass of a proton.
#define MASS_P 938.27 // MeV/c^2

/// The mass of a neutron.
#define MASS_N 939.57 // MeV/c^2

/// The mass of a nucleon.
// #define NUCLEON_MASS (MASS_P + MASS_N) / 2.0 // MeV/c^2
#define NUCLEON_MASS (HBARC2 / 41.47) // MeV/c^2

/// The \f$R_0\f$ quantity used in the formula \f$R = R_0 A^{1/3}\f$.
#define R_0 1.2

/// Proton sq. charge radius
/*
*   G. G. Simon, C. Schmidt, F. Borkowski, and V. H. Walter
*   Nucl. Phys. A 333, 381 ~1980!;
*/
#define PROTON_RADIUS_SQ 0.769 // fm^2

/// Neutron sq. charge radius
/*
*    Atac, H., Constantinou, M., Meziani, ZE. et al.
*    Nat Commun 12, 1759 (2021).
*/
#define NEUTRON_RADIUS_SQ -0.1152 // ± 0.0017 fm^2

/// Darwin-Foldy contribution
/*
*    J. L. Friar, J. Martorell, and D. W. L. Sprung
*    Phys. Rev. A 56, 4579
*/
#define DARWINFOLDY 0.033 // fm^2

/// Authors of the program.
#define CFG_AUTHORS "N. Dubray, J.-P. Ebran, P. Carpentier, M. Frosini, A. Zdeb, N. Pillet, J. Newsome, M. Verrière, G. Accorto, D. Regnier"

/// Constant used for equality test between doubles.
#define EPSILON 1e-9

/// Neutron id
#define NEUTRON 0

/// Proton id
#define PROTON 1

/// Total id
#define TOTAL 2

//==================
//===== Macros =====
//==================

/// Common ranges of quantum numbers
#define TRANGE  {0, 2}                    // iso
#define SRANGE  {0, 2}                    // spin
#define SSRANGE {0, 4}                    // spin spin
#define MRANGE  {-basis.mMax, basis.mMax} // m
#define MPRANGE {0, basis.mMax}           // m>0
#define NRANGE  {0, basis.nGlobalMax}     // n
#define NZRANGE {0, basis.n_zGlobalMax}   // nz
#define DRANGE  {0, 2}                    // d
#define DDRANGE {0, 4}                    // dd

//==================
//===== Macros =====
//==================

/// Define MIN macro.
#define MIN(X,Y) (X<Y?X:Y)

/// Define MAX macro.
#define MAX(X,Y) (X>Y?X:Y)

//===================
//===== Options =====
//===================

/// Print Berger's basis deformation parameters.
// #define PRINT_BERGER_BASIS_PARAMETERS

/// Use a cache mechanism for orthogonal polynomials.
#define USE_POLY_CACHING

/// New / old way of calculating local densities.
#define NEW_DENSIT

/// Ignore the constraint on \<Q10\> during Woods-Saxon optimizations.
#define WS_IGNORE_Q10_CONSTRAINT

/// Use N. Schunck's method to calculate the neck's abcissa.
// #define USE_NECK_SCHUNCK

/// Include Plotly.js source in binary.
#define INCLUDE_PLOTLYJS_IN_BINARY

/// Make the generated VTK files valid XML files.
//#define VTK_EXPORT_USE_BASE64

/// Clamp chemical potential values.
// #define CLAMP_CHEMPOTS

/// Short info() by default (used by SWIG).
#define USE_SHORT_INFO false

/// Check bounds of FMULTI objects.
// #define CHECKBOUNDS_FMULTI

/// Use ANSI "256 colors".
#define ANSI_256_COLORS

/// Use ANSI trueColors (if supported by terminal).
// #define ANSI_TRUECOLORS

//================
//===== Misc =====
//================

///@{
/** Format strings for table construction.
 */
#define TABLE_TD          std::string("##TD")
#define TABLE_TR          std::string("##TR")
#define TABLE_CENTER      std::string("##CE")
#define TABLE_LEFT        std::string("##LE")
#define TABLE_RIGHT       std::string("##RI")
#define TABLE_RED         std::string("##RE")
#define TABLE_GREEN       std::string("##GR")
#define TABLE_BLUE        std::string("##BL")
#define TABLE_YELLOW      std::string("##YE")
#define TABLE_NORM        std::string("##NO")
///@}

#ifdef SKILL0
/// SKILL level 0 (cf. Makefile for meaning)
#define SKILL "SKILL0 (I'm too young to die)"
#endif

#ifdef SKILL1
/// SKILL level 1 (cf. Makefile for meaning)
#define SKILL "SKILL1 (Hey, not too rough)"
#endif

#ifdef SKILL2
/// SKILL level 2 (cf. Makefile for meaning)
#define SKILL "SKILL2 (Hurt me plenty)"
#endif

#ifdef SKILL3
/// SKILL level 3 (cf. Makefile for meaning)
#define SKILL "SKILL3 (Ultra-Violence)"
#endif

#ifdef SKILL4
/// SKILL level 4 (cf. Makefile for meaning)
#define SKILL "SKILL4 (Nightmare!)"
#endif

// Default SKILL string value (should be set above).
#ifndef SKILL
/// No SKILL defined !?
#define SKILL "NO SKILL DEFINED !?"
#endif

//==============================================================================
//==============================================================================
//==============================================================================

#include "logger.h"
#include "general.h"

//==============================================================================
//==============================================================================
//==============================================================================

/// Categories of messages that must be displayed on the standard output.
extern std::set<INT> msgToOut;

/// Possible verbosity level values.
enum
{
  MSG_ERROR,
  MSG_INFO,
  MSG_DATA,
  MSG_TIMELINE,
  MSG_WARNING,
  MSG_DEBUG,
  MSG_PID,
  MSG_MAIN,
  MSG_TIME
};

/// Force the loading of an invalid DataTree object.
extern bool forceInvalidDataTree;

/// Save the DataTree to a file.
extern bool saveToFile;

/// Use colors in output.
extern bool useColors;

/// Generate bokeh plots.
extern bool useBokeh;

/// A Logger instance.
extern Logger logger;

/// A General instance.
extern General general;

/// How to react in case of an error (throw exception or exit).
extern bool exitOnError;

/// Style of the printed tables.
extern int tableStyle;

/// Use UTF-8 characters in the output.
extern bool useUtf8;

// /// Convenient sqrt version with explicit conversion.
// extern double sqrt(const INT);
// extern double sqrt(const int);
// extern double pow(const double, const INT);
// extern double pow(const double, const int);

/** Equality operators for armadillo objects.
 * @{
 */

extern bool operator==(const UVEC &, const UVEC &);
extern bool operator==(const IVEC &, const IVEC &);
extern bool operator==(const arma::vec &, const arma::vec &);

extern bool operator==(const UMAT &, const UMAT &);
extern bool operator==(const IMAT &, const IMAT &);
extern bool operator==(const arma::mat &, const arma::mat &);

extern bool operator==(const UCUBE &, const UCUBE &);
extern bool operator==(const ICUBE &, const ICUBE &);
extern bool operator==(const arma::cube &, const arma::cube &);

extern bool operator!=(const UVEC &, const UVEC &);
extern bool operator!=(const IVEC &, const IVEC &);
extern bool operator!=(const arma::vec &, const arma::vec &);

extern bool operator!=(const UMAT &, const UMAT &);
extern bool operator!=(const IMAT &, const IMAT &);
extern bool operator!=(const arma::mat &, const arma::mat &);

extern bool operator!=(const UCUBE &, const UCUBE &);
extern bool operator!=(const ICUBE &, const ICUBE &);
extern bool operator!=(const arma::cube &, const arma::cube &);

/**
 * @}
 */

#endif // GLOBAL_H
