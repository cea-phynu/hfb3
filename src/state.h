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

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define STATE_VALID_KEYS \
{ "state/U"                          , "U matrices"                         , ""     , "MM" }, \
{ "state/V"                          , "V matrices"                         , ""     , "MM" }, \
{ "state/rho"                        , "Rho matrices"                       , ""     , "MM" }, \
{ "state/kappa"                      , "Kappa matrices"                     , ""     , "MM" }, \
{ "state/chemicalPotential"          , "Chemical potentials"                , ""     , "V"  }, \
{ "state/totalEnergy"                , "Total energy"                       , ""     , "D"  }, \
{ "state/converged"                  , "Converged state ?"                  , "False", "B"  }, \
{ "state/nbIter"                     , "Number of iterations"               , "0"    , "I"  }, \
{ "state/blockedQPsIndex"            , "Index of the blocked QP states"     , ""     , "IV" }, \
{ "state/blockedQPsOmega"            , "Omega of the blocked QP states"     , ""     , "IV" }, \
{ "state/blockedQPsIsospin"          , "Isospin of the blocked QP states"   , ""     , "IV" }, \
{ "state/qpStatesEnergy"             , "Energy of the quasi-particle states", ""     , "MV" }, \
{ "state/qpStatesOccupation"         , "Occupation of individual states"    , ""     , "MV" }, \
{ "state/qpStatesIndex"              , "Index of the quasi-particle states" , ""     , "MIV"}, \
{ "state/calculationLength"          , "Length of the calculation [s]"      , "0.0"  , "D"  }, \
{ "basis/useStateBasis"              , "If possible, use the basis of the starting state", "True" , "B"  }

//==============================================================================
//==============================================================================
//==============================================================================

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
  void calcInertia(const std::string &interactionName = "", const std::list<std::pair<INT, INT> > &_collectiveCoordinates = {});
  bool checkSolution(void) const;                                      // #TEST#
  void calcCanonical(void);                                            // #TEST#
  const arma::vec getCanonicalV2(const arma::mat &_rho);
  void calcRhoKappaFromUV(void);                                       // #TEST#
  void calcUVFromRhoKappa(void);                                       // #TEST#
  const std::string getOmegaContributionsInfo(void) const;             // #TEST#
  void printCanonicalStatesInfo(const arma::mat &_rho, INT iso);       // #TEST#

  //============================================================================
  //============================================================================
  //============================================================================

  /// A state is given in a specified Basis
  Basis basis;

  /// A state is associated to a physical System
  System sys;

  /// A state is computed for a set of contraints
  std::map<std::string, Constraint> constraints;

  /// Blocked quasi-particle states
  std::set<StateId> blockedQPStates;

  // TODO: store omega-blocks only ?
  /// The rho matrices [NEUTRON, PROTON]
  Multi<arma::mat> rho;

  /// The kappa matrices [NEUTRON, PROTON]
  Multi<arma::mat> kappa;

  /// An arma::mat object representing \f$U\f$ in the HO*OR representation [NEUTRON, PROTON]
  Multi<arma::mat> U;

  /// An arma::mat object representing \f$V\f$ in the HO*OR representation [NEUTRON, PROTON]
  Multi<arma::mat> V;

  /// HF state omega and index [NEUTRON, PROTON]
  Multi<IMAT> oaiHF;

  /// Chemical potentials [NEUTRON, PROTON]
  arma::vec chemPot = arma::zeros(2);

  /// Total binding energy.
  double totalEnergy = 1e99;

  /// Is this state converged ?
  bool converged = false;

  /// Number of iterations performed
  INT nbIter = 0;

  /// Length of the calculation [s]
  double calculationLength = 0.0;

  //===== collective quantities =====

  /// collective coordinates for the inertia calculation
  // std::list<std::pair<INT, INT> > collectiveCoordinates = {{2, 0}, {3, 0}, {4, 0}};
  std::list<std::pair<INT, INT> > collectiveCoordinates = {{2, 0}, {3, 0}, {4, 0}, {2, 1}, {2, 2}};

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
