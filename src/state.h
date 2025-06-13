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

#ifndef STATE_H
#define STATE_H

/** \file
 *  \brief Headers for the State class.
 */

#include "global.h"
#include "generic.h"
#include "basis.h"
#include "system.h"
#include "constraint.h"
#include "states.h"

class DataTree;

/** \brief Cylindrical one- or two-center state.
 *
 *  This class represents a cylindrical one- or two-center state (Bogoliubov- or EFA-state).
 */

class State : public Generic
{
public :
  State(void);                                                         // #TEST#
  State(const std::string &filename);                                  // #TEST#
  State(const DataTree &dataTree);                                     // #TEST#

  DataTree getDataTree(void);                                          // #TEST#
  void convertFrom(Basis &);                                           // #TEST#
  State convertTo(Basis &);                                            // #TEST#
  double getDensityDistance(State                                      // #TEST#
                            &otherState,
                            bool pairing = false);
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string getNiceInfo(const std::string &what = "");         // #TEST#
  bool empty(void) const;                                              // #TEST#
  const arma::vec getOverlap(State &otherState);   // #TEST#
  void calcInertia(const IVEC &collectiveCoordinates = {});            // #TEST#
  bool checkSolution(void) const;                                      // #TEST#
  void calcCanonical(void);                                            // #TEST#
  const arma::vec getCanonicalV2(const arma::mat &_rho);
  void calcRhoKappaFromUV(void);                                       // #TEST#
  void calcUVFromRhoKappa(const DataTree &dataTree);                   // #TEST#
  const std::string getOmegaContributionsInfo(void) const;             // #TEST#
  void printCanonicalStatesInfo(const arma::mat &_rho, INT iso);       // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// A state is given in a specified Basis
  Basis basis;

  /// A state is associated to a physical System
  System sys;

  /// A state is computed for a set of contraints
  std::map<std::string, Constraint> constraints;

  /// Blocked quasi-particle states [NEUTRON, PROTON].
  Multi<INT> blockedQP;

  // TODO: store omega-blocks only
  /// The rho matrices [NEUTRON, PROTON].
  Multi<arma::mat> rho;

  /// The kappa matrices [NEUTRON, PROTON].
  Multi<arma::mat> kappa;

  /// An arma::mat object representing \f$U\f$ in the HO*OR representation [NEUTRON, PROTON].
  Multi<arma::mat> U;

  /// An arma::mat object representing \f$V\f$ in the HO*OR representation [NEUTRON, PROTON].
  Multi<arma::mat> V;

  /// Occupation numbers \f$v^2\f$ [NEUTRON, PROTON].
  Multi<arma::vec> vecOc;

  /// Individual HF energies [NEUTRON, PROTON].
  Multi<arma::vec> eneQP;

  /// HF state omega and index [NEUTRON, PROTON].
  Multi<IMAT> oaiHF;

  /// Chemical potentials [NEUTRON, PROTON].
  arma::vec chemPot = arma::zeros(2);

  /// Total binding energy.
  double totalEnergy = 1e99;

  /// Is this state converged ?
  bool converged = false;

  /// Number of iterations performed.
  INT nbIter = 0;

  /// Length of the calculation [s].
  double calculationLength = 0.0;

  //===== collective quantities =====

  /// collective coordinates
  IVEC collectiveCoordinates;

  /// metric (GCM and ATDHF)
  arma::mat metric;

  /// masses (GCM)
  arma::mat massGCM;

  /// ZPE energies (GCM)
  arma::mat zpeGCM;

  /// masses (ATDHF)
  arma::mat massATDHF;

  /// ZPE energies (ATDHF)
  arma::mat zpeATDHF;

  //===== Canonical objects (in OR representation) =====

  /// u_c quantities
  Multi<arma::vec> u_c;

  /// v_c quantities
  Multi<arma::vec> v_c;

  /// D transformation matrix (OR -> canonical)
  Multi<arma::mat> D_c;

  std::map<std::string, double> energyContributions;

  //===== Individual states =====

  /// Particle states (iso)
  Multi<States> partStates;

  /// Quasi-particle states (iso)
  Multi<States> qpStates;
};

#endif // STATE_H
