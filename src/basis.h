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

#ifndef BASIS_H
#define BASIS_H

/** \file
 *  \brief Headers for the Basis class.
 */

#include "global.h"
#include "generic.h"
#include "datatree.h"
#include "qnumbers.h"
#include "multi.h"

/** \brief Cylindrical, one- or two-center basis.
 *
 * This class implements a cylindrical one- or two-center basis.
 */

class Basis : public Generic
{
public :

  // Constructors
  Basis(void);
  explicit Basis(const std::string &filename);                         // #TEST#
  explicit Basis(const DataTree &dataTree,                             // #TEST#
                 const std::string &prefix = "");

  explicit Basis(double _d_0,                                          // #TEST#
                 double _b_r,
                 double _b_z,
                 INT _nOscil,
                 INT _n_zMaxImposed,
                 double _g_q);

  static Basis fromBerger2ct(INT _nOscil,                              // #TEST#
                         double _g_q,
                         double _hw,
                         double _qzr,
                         double _d_0,
                         INT _n_zMaxImposed);

  // Basis truncation / optimization
  void cylTruncate(void);                                              // #TEST#

  // Polynomial evaluation
  static double hermite(INT n_z, double z);                            // #TEST#
  static double laguerre(INT m, INT n, double r);                      // #TEST#
  static double legendre(INT l, double x);                             // #TEST#
  static arma::vec hermite(INT n_z, const arma::vec &zeta);            // #TEST#
  static arma::vec laguerre(INT m,                                     // #TEST#
                            INT n,
                            const arma::vec &eta);
  static arma::vec legendre(INT l, const arma::vec &x);                // #TEST#

  // Wave function evaluation
  double zPartScalar(double z, INT n_z, INT d = -1);                   // #TEST#
  double rPartScalar(double r, INT m, INT n);                          // #TEST#
  arma::vec zPart(const arma::vec &z, INT n_z, INT d = -1);            // #TEST#
  arma::vec rPart(const arma::vec &r, INT m, INT n);                   // #TEST#

  void debugMat(const arma::mat &mat) const;

  // Normalized wave function evaluation
  arma::vec rPartNorm(const arma::vec &eta, INT m, INT n);
  arma::vec zPartNorm(const arma::vec &zeta, INT n_z);
  arma::vec rPartNormReduced(const arma::vec &eta,
                             INT m,
                             INT n);
  arma::vec zPartNormReduced(const arma::vec &zeta, INT n_z);

  // 1st derivative of normalized wave function evaluation
  arma::vec rPartNormd(const arma::vec &eta,
                       INT m,
                       INT n,
                       INT ipm);
  arma::vec zPartNormd(const arma::vec &zeta, INT n_z);
  arma::vec rPartL0(const arma::vec &eta, INT m, INT n);
  arma::vec rPartLavecm(const arma::vec &eta, INT m, INT n);

  // Misc
  DataTree getDataTree(const std::string &prefix = "") const;          // #TEST#

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string niceStr(void) const;                               // #TEST#

  bool operator==(const Basis &other) const;                           // #TEST#
  bool operator!=(const Basis &other) const                            // #TEST#
  {
    return !((*this) == other);
  };

  // Basis orthogonalization
  void calcWDN(void);                                                  // #TEST#
  void calcRecouv(void);
  void calcHharmo2ct(void);                                            // #TEST#

  // Basis transformation
  arma::mat getTransformationMatrix(Basis &other);                     // #TEST#
  arma::mat getOverlapMatrixHOZ(Basis &other);                         // #TEST#
  arma::mat getOverlapMatrixHOR(Basis &other);                         // #TEST#
  double getBasisDistance(const Basis &otherBasis);                    // #TEST#

  double getOverlapHO(Basis &other,                                    // #TEST#
                      UINT iB,
                      UINT iK,
                      const arma::mat &overlapMatrixHOR,
                      const arma::mat &overlapMatrixHOZ);
  Multi<arma::mat> getOverlapMatrixHO(Basis &other);                // #TEST#

  // Basis overlap calculation
  void calcTab(void);                                                  // #TEST#
  double tabzd(INT nza, INT nzb, INT da, INT db) const;                // #TEST#

  void calcTalmanr(void);                                              // #TEST#
  void calcTalmanz(void);                                              // #TEST#
  void calcMoshinskyr(void);                                           // #TEST#
  void calcMoshinskyz(void);                                           // #TEST#

  arma::mat getFullOverlap(Basis &otherBasis);

  // The following method should ONLY BE USED FOR TESTING
  double getFullMoshinskyr(INT m1,                                     // #TEST#
                           INT n1,
                           INT m2,
                           INT n2,
                           INT M1,
                           INT N1,
                           INT M2,
                           INT N2) const;

  // double getMoshinskyr(INT m1, INT n1, INT m2, INT n2, INT N2) const;

  Multi<arma::vec> getFullMoshinskyz(void);
  double logMultinomial(INT a, INT b, INT c, INT d) const;

  // Only used for tests
  double binomialSign(INT n, INT m) const;
  double logBinomial(INT n, INT m) const;

  // Used only for Spin-Obrit testing. (4HOF INTegrals).
  void calcIz(void);
  void calcIr(void);

  // DG stuff
  void calcMoshinskyDG(void);                                          // #TEST#
  void calcMoshinskyTDG(void);                                         // #TEST#

  // Clebsch-Gordan coefficients
  double CGcoeff(INT j, INT m, INT j1, INT m1, INT j2, INT m2);

  const std::string getNiceInfo(void);                                 // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// The \f$d_0\f$ parameter.
  double d_0 = -1.0;

  /// The maximum value for \f$d\f$.
  INT dMax = 0;

  /// The maximum value for \f$s\f$.
  INT sMax = 0;

  /// The number of oscillator shells.
  INT nOscil = 11;

  /// The global maximum value of \f$n_z + 1\f$.
  INT n_zGlobalMax = 0;

  /// The global maximum value of \f$n + 1\f$.
  INT nGlobalMax = 0;

  /// The maximum value for \f$m\f$.
  INT mMax = 0;

  /// The maximum value for talman-r and moshinsky-r.
  INT Nmaxr = 0;

  /// The maximum values for \f$n\f$.
  IVEC nMax;

  /// The maximum values for \f$n_z\f$.
  IMAT n_zMax;

  /// The maximum imposed value for \f$n_z\f$.
  INT n_zMaxImposed = 24;

  ///  The "deformation" of the basis.
  double g_q = 1.0;

  /// The \f$b_\perp\f$ parameter.
  double b_r = 1.5;

  /// The \f$b_z\f$ parameter.
  double b_z = 1.5;

  /// The n_zMax values in berger2ct.
  IVEC n_zMaxBerger;

  /// The 'alf' variable in berger2ct.
  double alfBerger = 0.0;

  /// The 'bet' variable in berger2ct.
  double betBerger = 0.0;

  /// The 'q' variable in berger2ct.
  double qBerger = 0.0;

  /// The 'hwr' variable in berger2ct.
  double hwrBerger = 0.0;

  /// The 'hwz' variable in berger2ct.
  double hwzBerger = 0.0;

  /// The 'hw' variable in berger2ct.
  double hwBerger = 0.0;

  /// The indices of the omega blocks in the HO basis.
  Multi<UVEC> omegaIndexHO;

  /// The sizes of the m-blocks.
  UVEC mSize;

  /// The \f$W\f$ arma::mat.
  Multi<arma::mat> W;

  /// The \f$N\f$ arma::mat.
  Multi<arma::mat> N;

  /// The full overlap arma::mat.
  arma::mat S;

  /// The transformation from basis OR to basis HO.
  arma::mat ORtoHO;

  /// The \f$H\f$ arma::mat.
  Multi<arma::mat> H;

  /// The \f$T_{ab}\f$ arma::mat.
  arma::mat Tab;

  /// The Talman coefficients (\f$z\f$-part).
  Multi<arma::vec> talmanz;

  /// The Talman coefficients (\f$r\f$-part).
  Multi<arma::vec> talmanr;

  /// The Moshinsky coefficients (\f$r\f$-part).
  Multi<arma::vec> moshinskyr;

  /// The Moshinsky coefficients (\f$z\f$-part).
  arma::mat moshinskyz;

  /// The set of quantum numbers for the HO basis.
  Qnumbers HOqn;

  /// The set of quantum numbers for the Hartree-Fock (HF) basis.
  Qnumbers HFqn;

  /// The set of quantum numbers {m, n}.
  Qnumbers QNmn;

  /// The set of quantum numbers {m, n, n_z, d}.
  Qnumbers QNmnn_zd;

  /// The set of quantum numbers {n_z}.
  Qnumbers QNn_zMax;

  /// The set of quantum numbers {n_z, d}.
  Qnumbers QNn_zMaxd;

  /// The set of quantum numbers for the Orthogonal (OR) basis.
  Qnumbers ORqn;

  /// The transformation from basis HO to basis OR.
  arma::mat HOtoOR;

  /// The indices of the omega blocks in the OR basis (Cylindrical basis only).
  Multi<UVEC> omegaIndexOR;

  /// The incomplete z-overlaps as defined in Eq. (E-16), p. E-2, PhD J.-F. Berger.
  arma::mat Recouvz;

  /// The incomplete overlaps as defined in  Eq. (E-17), p. E-3, PhD J.-F. Berger.
  arma::mat Recouv;

  /// 4 HO-r overlap. For testing only.
  Multi<arma::mat> Ir;

  /// 4 HO-r overlap (nabla + on the last indice). For testing only.
  Multi<arma::mat> Ir_plus;

  /// 4 HO-r overlap (nabla - on the last indice). For testing only.
  Multi<arma::mat> Ir_minus;

  /// 4 HO-r overlap. For testing only.
  Multi<arma::mat> Ir_mm;

  /// 4 HO-r overlap. For testing only.
  Multi<arma::mat> Ir_pp;

  //4 HO-z overlap. For testing only.
  Multi<arma::cube> Iz;

  //4 HO-z overlap (nabla 0 on the last indice). For testing only.
  Multi<arma::cube> NablaIz;

  //For testing only.
  Multi<arma::vec> Jz;

  //Moshinsky-r for the SO INTeraction
  Multi<arma::vec> moshinskyDG;

  //Moshinsky-r for the TS INTeraction
  Multi<arma::vec> moshinskyTDG;

  //Global max for n_sigma (SO INTeraction)
  INT NmaxDG = 0;

  //Global max for n_sigma (TS INTeraction)
  INT NmaxTDG = 0;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  // Auxiliary function for cylTruncate().
  INT calcNumberOfStates(INT n,
                         double p,
                         double q,
                         INT nx,
                         INT ny,
                         INT nz) const;

  // Auxiliary function for getOverlapMatrixHO().
  Multi<UVEC> calcMSBlocks(Basis &b);
};

#endif // BASIS_H
