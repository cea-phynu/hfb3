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

#include "state.h"
#include "global.h"
#include "tools.h"
#include "discrete.h"
#include "multipole_operators.h"
#include "interaction.h"
#include "solver_hfb_broyden.h"

/** \file
 *  \brief Methods of the State class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > State::validKeys =
  {
    { "state/rho"                        , "Rho matrices"                   , ""     , "MM" },
    { "state/kappa"                      , "Kappa matrices"                 , ""     , "MM" },
    { "state/individualStates/energy"    , "Energy of individual states"    , ""     , "MV" },
    { "state/individualStates/occupation", "Occupation of individual states", ""     , "MV" },
    { "state/chemicalPotential"          , "Chemical potentials"            , ""     , "V"  },
    { "state/totalEnergy"                , "Total energy"                   , ""     , "D"  },
    { "state/converged"                  , "Converged state ?"              , "False", "B"  },
    { "state/nbIter"                     , "Number of iterations"           , "0"    , "I"  },
    { "state/calculationLength"          , "Length of the calculation [s]"  , "0.0"  , "D"  },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor.
 */

State::State(void)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Construct a State from a filename.
 *
 *  \param filename A filename.
 */

State::State(const std::string &filename) : State(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Construct a State from a DataTree object.
 *
 *  \param dataTree A DataTree object.
 */

State::State(const DataTree &dataTree) :
  basis(Basis(dataTree)),
  sys(System(dataTree))
{
  DBG_ENTER;

  rho(NEUTRON) = arma::mat();
  rho(PROTON ) = arma::mat();
  kappa(NEUTRON) = arma::mat();
  kappa(PROTON ) = arma::mat();
  U(NEUTRON) = arma::mat();
  U(PROTON ) = arma::mat();
  V(NEUTRON) = arma::mat();
  V(PROTON ) = arma::mat();
  blockedQP(NEUTRON) = -1;
  blockedQP(PROTON ) = -1;

  // Try to construct rho and kappa from the dataTree instance
  Multi<arma::mat> multiRho;
  Multi<arma::mat> multiKappa;
  dataTree.get(multiRho,               "state/rho", true);
  dataTree.get(multiKappa,             "state/kappa", true);

  if ((multiRho.size() > 2) && (multiKappa.size() > 2))
  {
    for (INT iso: {NEUTRON, PROTON})
    {
      rho(iso) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
      kappa(iso) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

      for (INT omega = 0; omega < basis.mMax; omega++)
      {
        // INFO("%d %d", iso, omega);
        // Tools::info("multiRho", multiRho(iso, omega));
        // Tools::info("multiKappa", multiKappa(iso, omega));

        rho(  iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) = multiRho(iso, omega);
        kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) = multiKappa(iso, omega);
      }
    }
  }
  else if ((multiRho.size() > 0) && (multiKappa.size() > 0))
  {
    // try to load using previous storage format. TODO: to be removed
    dataTree.get(rho,             "state/rho", true);
    dataTree.get(kappa,           "state/kappa", true);
  }

  dataTree.get(eneQP,             "state/individualStates/energy", true);
  dataTree.get(chemPot,           "state/chemicalPotential", true);

  dataTree.get(totalEnergy,       "state/totalEnergy", true);
  dataTree.get(converged,         "state/converged", true);

  dataTree.get(nbIter,            "state/nbIter", true);
  dataTree.get(calculationLength, "state/calculationLength", true);

  // Construct constraints
  constraints = Constraint::fromDataTree(dataTree);

  basis.calcWDN();

  // If a basis is provided, use it (unless specified otherwise)
  if (dataTree.contains("state/basis/b_r"))
  {
    Basis stateBasis = Basis(dataTree, "state/");

    bool useStateBasis = true;
    dataTree.get(useStateBasis, "basis/useStateBasis", true);

    if (stateBasis != basis)
    {

      if (useStateBasis)
      {
        basis = stateBasis;
        basis.calcWDN();

        Tools::debug("The basis specified in 'state/basis/' will be used.");
        Tools::debug("To use the basis specified in 'basis/', set `basis/useStateBasis` to False.");
      }
      else
      {
        Tools::info("Converting the state from the basis in `state/basis/` to the basis in `basis/`.");

        // Convert the state FROM the stateBasis to the current basis
        convertFrom(stateBasis);
      }
    }
  }

  // If missing rho and kappa, try to recover them from U and V
  if (rho(NEUTRON).empty() || rho(PROTON).empty() ||
      kappa(NEUTRON).empty() || kappa(PROTON).empty())
  {
    calcRhoKappaFromUV();
  }

  // // If missing U and V, try to recover them from rho and kappa
  // if (U(NEUTRON).empty() || U(PROTON).empty() ||
  //     V(NEUTRON).empty() || V(PROTON).empty())
  // {
  //   calcUVFromRhoKappa(dataTree);
  // }

  // Tools::mesg("BogCnv", getOmegaContributionsInfo());

  // Set collectiveCoordinates (used for inertia calculation)
  for (auto &c : constraints)
    if (c.second.useForInertias)
      Tools::growIVec(collectiveCoordinates, c.second.lm);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Update a DataTree instance.
 */

DataTree State::getDataTree(void)
{
  DBG_ENTER;

  DataTree state;

  // State datatree must contain the associated basis
  state.merge(basis.getDataTree("state/"));

  // State datatree must contain the associated system
  state.merge(sys.getDataTree());

  // Save rho and kappa matrices
  Multi<arma::mat> multiRho;
  Multi<arma::mat> multiKappa;
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < basis.mMax; omega++)
    {

      arma::mat rhoHO   = rho(  iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
      arma::mat kappaHO = kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));

      multiRho(iso, omega) = rhoHO;
      multiKappa(iso, omega) = kappaHO;

      // INFO("%d %d", iso, omega);
      // Tools::info("multiRho", multiRho(iso, omega));
      // Tools::info("multiKappa", multiKappa(iso, omega));
    }
  }

  state.set("state/rho",                         multiRho);
  state.set("state/kappa",                       multiKappa);

  state.set("state/individualStates/energy",     eneQP);
  state.set("state/individualStates/occupation", vecOc);
  state.set("state/chemicalPotential",           chemPot);
  state.set("state/totalEnergy",                 totalEnergy);
  state.set("state/converged",                   converged);
  state.set("state/nbIter",                      nbIter);
  state.set("state/calculationLength",           calculationLength);

  // set the energy constraints
  for (auto &c : constraints)
  {
    state.merge(c.second.getDataTree());
  }

  DBG_RETURN(state);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check for emptiness.
 */

bool State::empty(void) const
{
  for (INT iso: {NEUTRON, PROTON})
  {
    if (!rho(iso).empty()) return false;
    if (!kappa(iso).empty()) return false;
    if (!U(iso).empty()) return false;
    if (!V(iso).empty()) return false;
  }

  return true;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert from a different basis.
 *
 *  \param initialBasis The initial Basis object.
 */

void State::convertFrom(Basis &initialBasis)
{
  DBG_ENTER;

  if (basis == initialBasis)
  {
    // Tools::mesg("State.", PF_YELLOW("Identical bases => no conversion"));

    DBG_LEAVE;
  }

  if (empty())
  {
    // Tools::mesg("State.", PF_YELLOW("empty state => no conversion"));

    DBG_LEAVE;
  }

  arma::mat transMat = initialBasis.getTransformationMatrix(basis);

  // basis.calcWDN();
  // arma::mat Mitotal = basis.ORtoHO;
  // initialBasis.calcWDN();
  // arma::mat Mtotal  = initialBasis.HOtoOR;
  // U(NEUTRON) = transMat * U(NEUTRON) * Mtotal * transMat.t() * Mitotal;
  // U(PROTON ) = transMat * U(PROTON ) * Mtotal * transMat.t() * Mitotal;
  // V(NEUTRON) = transMat * V(NEUTRON) * Mtotal * transMat.t() * Mitotal;
  // V(PROTON ) = transMat * V(PROTON ) * Mtotal * transMat.t() * Mitotal;

  rho(  NEUTRON) = transMat * rho(  NEUTRON) * transMat.t();
  rho(  PROTON ) = transMat * rho(  PROTON ) * transMat.t();
  kappa(NEUTRON) = transMat * kappa(NEUTRON) * transMat.t();
  kappa(PROTON ) = transMat * kappa(PROTON ) * transMat.t();

  // basis conversion
  Tools::mesg("State.", "State basis conversion:");
  Tools::mesg("State.", "from: " + initialBasis.info(true));
  Tools::mesg("State.", "to:   " + basis.info(true));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Convert to a different basis.
 *
 *  \param initialBasis The target Basis object.
 */

State State::convertTo(Basis &targetBasis)
{
  DBG_ENTER;

  State state = (*this);
  state.basis = targetBasis;
  state.convertFrom(basis);

  DBG_RETURN(state);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the overlap with another state.
 */

const arma::vec State::getOverlap(State &otherState)
{
  DBG_ENTER;

  //Dependencies
  basis.calcWDN();
  otherState.basis.calcWDN();

  // in case of missing U and V matrices, perform one HFB Broyden iteration.
  if (U(NEUTRON).empty() || U(PROTON).empty() || V(NEUTRON).empty() || V(PROTON).empty())
  {
    SolverHFBBroyden solver(getDataTree());
    solver.nextIter();
    U = solver.state.U;
    V = solver.state.V;
  }

  // in case of missing U and V matrices, perform one HFB Broyden iteration.
  if (otherState.U(NEUTRON).empty() ||
      otherState.U(PROTON ).empty() ||
      otherState.V(NEUTRON).empty() ||
      otherState.V(PROTON ).empty())
  {
    SolverHFBBroyden solver(otherState.getDataTree());
    solver.nextIter();
    otherState.U = solver.state.U;
    otherState.V = solver.state.V;
  }

  //Useful quantities
  Basis &basis1 = otherState.basis;
  arma::mat R = basis.HOtoOR * basis.getFullOverlap(basis1) * basis1.HOtoOR.t(); //Overlap in ORxOR basis.
  arma::mat Mt_0 = basis.ORtoHO.t();
  arma::mat Mt_1 = basis1.ORtoHO.t();
  arma::vec overlap(3);
  overlap(0) = 1;
  overlap(1) = 1;

  for (INT omega = 0 ; omega < std::min(basis.mMax, basis1.mMax) ; omega++)
  {
    arma::mat R_omega = R.submat(basis.omegaIndexOR(omega), basis1.omegaIndexOR(omega));
    arma::mat Mt_omega_0 = Mt_0.submat(basis.omegaIndexOR(omega), basis.omegaIndexHO(omega));
    arma::mat Mt_omega_1 = Mt_1.submat(basis1.omegaIndexOR(omega), basis1.omegaIndexHO(omega));

    Multi<arma::mat> U_0;
    Multi<arma::mat> V_0;
    Multi<arma::mat> U_1;
    Multi<arma::mat> V_1;

    U_0(NEUTRON) = Mt_omega_0 * U(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
    U_0(PROTON ) = Mt_omega_0 * U(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
    V_0(NEUTRON) = Mt_omega_0 * V(NEUTRON).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
    V_0(PROTON ) = Mt_omega_0 * V(PROTON ).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));

    U_1(NEUTRON) = Mt_omega_1 * otherState.U(NEUTRON).submat(basis1.omegaIndexHO(omega), basis1.omegaIndexOR(omega));
    U_1(PROTON ) = Mt_omega_1 * otherState.U(PROTON ).submat(basis1.omegaIndexHO(omega), basis1.omegaIndexOR(omega));
    V_1(NEUTRON) = Mt_omega_1 * otherState.V(NEUTRON).submat(basis1.omegaIndexHO(omega), basis1.omegaIndexOR(omega));
    V_1(PROTON ) = Mt_omega_1 * otherState.V(PROTON ).submat(basis1.omegaIndexHO(omega), basis1.omegaIndexOR(omega));

    for (INT iso: {NEUTRON, PROTON})
    {
      arma::mat U_1_inv = arma::inv(U_1(iso));

      double det_plus = fabs(arma::det(U_0(iso).t() + V_0(iso).t() * R_omega * V_1(iso) * U_1_inv * R_omega.t()));
      double det_minus = fabs(arma::det(U_0(iso).t() - V_0(iso).t() * R_omega * V_1(iso) * U_1_inv * R_omega.t()));

      overlap(iso) = fabs(arma::det(U_1(iso))) * MAX(det_plus, det_minus) * overlap(iso);
    } //iso
  } //omega

  overlap(TOTAL) = overlap(NEUTRON) * overlap(PROTON);

  DBG_RETURN(overlap);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate U and V matrices from rho and kappa matrices IF NEEDED.
 */

void State::calcUVFromRhoKappa(const DataTree &dataTree)
{
  DBG_ENTER;

  if (rho(NEUTRON).empty()) DBG_LEAVE;
  if (rho(PROTON ).empty()) DBG_LEAVE;
  if (kappa(NEUTRON).empty()) DBG_LEAVE;
  if (kappa(PROTON ).empty()) DBG_LEAVE;

  if (!U(NEUTRON).empty()) DBG_LEAVE;
  if (!U(PROTON ).empty()) DBG_LEAVE;
  if (!V(NEUTRON).empty()) DBG_LEAVE;
  if (!V(PROTON ).empty()) DBG_LEAVE;

  Tools::mesg("State.", "Updating U and V matrices from rho and kappa matrices");

  // perform a single SolverHFBBroyden iteration to generate U and V matrices
  SolverHFBBroyden solverHFBBroyden(dataTree, *this);
  solverHFBBroyden.nextIter();

  U = solverHFBBroyden.state.U;
  V = solverHFBBroyden.state.V;

  //============================================================================
  //============================================================================
  //============================================================================

  // TODO: Other way to do the same thing using Cholesky decomposition ?

  // // J.-P. Ebran's derivations, p. 8, Eq. (II-B.18).
  // arma::mat Mtot = state.basis.HOtoOR.t();
  //
  // Multi<arma::mat> rho;
  // rho(NEUTRON) = state.rho(NEUTRON);
  // rho(PROTON ) = state.rho(PROTON );
  //
  // Multi<arma::mat> kappa;
  // kappa(NEUTRON) = state.kappa(NEUTRON);
  // kappa(PROTON ) = state.kappa(PROTON );
  //
  // Multi<arma::mat> Utot;
  // Multi<arma::mat> Vtot;
  //
  // Utot(PROTON ) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  // Utot(NEUTRON) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  // Vtot(NEUTRON) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  // Vtot(PROTON ) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  //
  // for (INT iso: {NEUTRON, PROTON})
  // {
  //   for (INT omega = 0; omega < state.basis.mMax; omega++)
  //   {
  //     arma::mat M = Mtot.submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexOR(omega));
  //
  //     arma::mat rhoHO   = rho(  iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega));
  //     arma::mat kappaHO = kappa(iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega));
  //
  //     arma::mat rhoOR   = M.t() * rhoHO   * M;
  //     arma::mat kappaOR = M.t() * kappaHO * M;
  //
  //     UINT ORsize_omega = state.basis.omegaIndexOR(omega).n_rows;
  //     arma::mat RmatOR = arma::zeros(ORsize_omega * 2, ORsize_omega * 2);
  //
  //     // J.-P. Ebran's derivations, p. 8, Eq. (II-C.14).
  //     RmatOR.submat(           0,            0,     ORsize_omega - 1,     ORsize_omega - 1) =  rhoOR;
  //     RmatOR.submat(ORsize_omega, ORsize_omega, 2 * ORsize_omega - 1, 2 * ORsize_omega - 1) = arma::eye(ORsize_omega, ORsize_omega) - rhoOR;
  //     RmatOR.submat(           0, ORsize_omega,     ORsize_omega - 1, 2 * ORsize_omega - 1) = -kappaOR;
  //     RmatOR.submat(ORsize_omega,            0, 2 * ORsize_omega - 1,     ORsize_omega - 1) = -kappaOR.t();
  //
  //     arma::vec E;
  //     arma::mat B;
  //
  //     if (!Tools::checkSymmetry(RmatOR, "RmatOR in " + std::string(__PRETTY_FUNCTION__)))
  //       RmatOR = arma::symmatu(RmatOR);
  //
  //     Tools::eig_sym(E, B, RmatOR);
  //
  //     // Extraction of U and V matrices from B.
  //     // J.-P. Ebran's derivations, p. 8, Eq. (II-C.4).
  //     arma::mat U = B.submat(0, 0, ORsize_omega - 1, ORsize_omega - 1);
  //     arma::mat V = B.submat(0, ORsize_omega, ORsize_omega - 1, 2 * ORsize_omega - 1);
  //
  //     Utot(iso).submat(state.basis.omegaIndexOR(omega), state.basis.omegaIndexOR(omega)) = U;
  //     Vtot(iso).submat(state.basis.omegaIndexOR(omega), state.basis.omegaIndexOR(omega)) = V;
  //   }
  // }
  //
  // arma::mat M = state.basis.HOtoOR.t();
  //
  // // in HO*OR representation
  // state.U(NEUTRON) = M * Utot(NEUTRON);
  // state.U(PROTON ) = M * Utot(PROTON );
  // state.V(NEUTRON) = M * Vtot(NEUTRON);
  // state.V(PROTON ) = M * Vtot(PROTON );

  //============================================================================
  //============================================================================
  //============================================================================

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate rho and kappa matrices from U and V matrices.
 */

void State::calcRhoKappaFromUV(void)
{
  DBG_ENTER;

  if (U(NEUTRON).empty()) DBG_LEAVE;
  if (U(PROTON ).empty()) DBG_LEAVE;
  if (V(NEUTRON).empty()) DBG_LEAVE;
  if (V(PROTON ).empty()) DBG_LEAVE;

  Tools::mesg("State.", "Updating rho and kappa matrices from U and V matrices");

  rho(NEUTRON)   = V(NEUTRON) * V(NEUTRON).t();
  rho(PROTON )   = V(PROTON ) * V(PROTON ).t();
  kappa(NEUTRON) = V(NEUTRON) * U(NEUTRON).t();
  kappa(PROTON ) = V(PROTON ) * U(PROTON ).t();

  rho(NEUTRON)   = arma::symmatu(rho(NEUTRON));
  rho(PROTON )   = arma::symmatu(rho(PROTON ));
  kappa(NEUTRON) = arma::symmatu(kappa(NEUTRON));
  kappa(PROTON ) = arma::symmatu(kappa(PROTON ));


  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the density distance from another state.
 */

double State::getDensityDistance(State &otherState, bool pairing)
{
  DBG_ENTER;

  Discrete d1(&basis);
  Discrete d2(&otherState.basis);

  arma::mat v1;
  arma::mat v2;

  if (pairing)
  {
    v1 =                      kappa(NEUTRON) +                      kappa(PROTON );
    v2 = otherState.kappa(NEUTRON) + otherState.kappa(PROTON );
  }
  else
  {
    v1 =                      rho(NEUTRON) +                      rho(PROTON );
    v2 = otherState.rho(NEUTRON) + otherState.rho(PROTON );
  }

  arma::mat dens1 = d1.getLocalXZ(v1, true);
  arma::mat dens2 = d2.getLocalXZ(v2, true);
  Mesh &m = d1.mesh;
  arma::mat r = Tools::matFromCol(m.ax.p, m.az.nb);
  arma::mat velem = r * 2.0 * PI;
  arma::mat wz = Tools::matFromRow(m.az.we.t(), m.ax.nb);
  arma::mat wx = Tools::matFromCol(m.ax.we, m.az.nb);
  arma::mat diff = arma::abs(dens1 - dens2) % wx % wz % velem;
  double state = arma::accu(diff);

  DBG_RETURN(state);
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string State::info(bool isShort) const
{
  DBG_ENTER;

  std::string state = "";

  if (isShort)
  {
    state += Tools::treeStr(
    {
      {"sys.", sys.info(true)},
      {"cst.", Constraint::niceStr(constraints, true)},
      {"cvg.", converged ? "true" : "false"},
      {"ene.", PF("%.3e", totalEnergy)},
    }, true);
  }
  else
  {
    state += Tools::treeStr(
    {
      {"State", ""},
      {"System", sys.info(false)},
      {"Basis ", basis.info(true)},
      {"Const.", Constraint::niceStr(constraints, true)},
      {"conver", Tools::infoStr(converged)},
      {"rho(n)", Tools::infoStr(rho(NEUTRON))},
      {"rho(p)", Tools::infoStr(rho(PROTON ))},
      {"kap(n)", Tools::infoStr(kappa(NEUTRON))},
      {"kap(p)", Tools::infoStr(kappa(PROTON ))},
      {"U(n)  ", Tools::infoStr(U(NEUTRON))},
      {"U(p)  ", Tools::infoStr(U(PROTON ))},
      {"V(n)  ", Tools::infoStr(V(NEUTRON))},
      {"V(p)  ", Tools::infoStr(V(PROTON ))},
      {"bloc.n", Tools::infoStr(blockedQP(NEUTRON))},
      {"bloc.p", Tools::infoStr(blockedQP(PROTON ))},
      {"eneTot", Tools::infoStr(totalEnergy)},
      {"chemPn", Tools::infoStr(chemPot(NEUTRON))},
      {"chemPp", Tools::infoStr(chemPot(PROTON ))},
      {"ZpeGcm", Tools::infoStr(arma::trace(zpeGCM))},
      {"ZpeAtd", Tools::infoStr(arma::trace(zpeATDHF))},
    }, false);
  }

  DBG_RETURN(state);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return info in a nice form.
 */

const std::string State::getNiceInfo(const std::string &what)
{
  DBG_ENTER;

  std::string state = "";

  std::list<std::string> names;
  std::list<std::string> units;
  std::list<std::list<std::string> > valuesMetric;
  std::list<std::list<std::string> > valuesMassGCM;
  std::list<std::list<std::string> > valuesZpeGCM;
  std::list<std::list<std::string> > valuesMassATDHF;
  std::list<std::list<std::string> > valuesZpeATDHF;

  if ((what == "inertia") && (!metric.empty()))
  {
    INT ic = 0;
    for (auto &c : collectiveCoordinates)
    {
      names.push_back(PF("q%01d0", c));
      units.push_back(PF("[fm%01d]", c));

      std::list<std::string> val;
      val.push_back(PF("q%01d0 [fm%01d]", c, c));
      for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", metric(ic, ic2)));
      valuesMetric.push_back(val);

      val.clear();
      val.push_back(PF("q%01d0 [fm%01d]", c, c));
      for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", massGCM(ic, ic2)));
      valuesMassGCM.push_back(val);

      val.clear();
      val.push_back(PF("q%01d0 [fm%01d]", c, c));
      for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", zpeGCM(ic, ic2)));
      valuesZpeGCM.push_back(val);

      val.clear();
      val.push_back(PF("q%01d0 [fm%01d]", c, c));
      for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", massATDHF(ic, ic2)));
      valuesMassATDHF.push_back(val);

      val.clear();
      val.push_back(PF("q%01d0 [fm%01d]", c, c));
      for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", zpeATDHF(ic, ic2)));
      valuesZpeATDHF.push_back(val);

      ic++;
    }

    state += Tools::valueTable("Metric"         , names, units, valuesMetric   ) + "\n";
    state += Tools::valueTable("Mass GCM"       , names, units, valuesMassGCM  ) + "\n";
    state += Tools::valueTable("ZPE GCM [MeV]"  , names, units, valuesZpeGCM   ) + "\n";
    state += Tools::valueTable("Mass ATDHF"     , names, units, valuesMassATDHF) + "\n";
    state += Tools::valueTable("ZPE ATDHF [MeV]", names, units, valuesZpeATDHF );
  }
  else
  {
    state += info();
  }

  DBG_RETURN(state);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the inertia tensor, the metric and the ZPE energy (ATDHF and GCM prescriptions).
 */

void State::calcInertia(const IVEC &_collectiveCoordinates)
{
  DBG_ENTER;

  ASSERT(!U(NEUTRON).empty(), "Empty U(NEUTRON)");
  ASSERT(!U(PROTON ).empty(), "Empty U(PROTON )");
  ASSERT(!V(NEUTRON).empty(), "Empty V(NEUTRON)");
  ASSERT(!V(PROTON ).empty(), "Empty V(PROTON )");

  if (!_collectiveCoordinates.empty()) collectiveCoordinates = _collectiveCoordinates;

  Multi<arma::mat> Qtot;

  for (auto &lambda : collectiveCoordinates)
  {
    Qtot(NEUTRON, lambda) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);
    Qtot(PROTON , lambda) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);
  }

  MultipoleOperators multipoleOperators(*this);
  multipoleOperators.calcQl0Matrices();

  Multi<arma::mat> Utot, Vtot;

  // go to OR*OR representation
  Utot(NEUTRON) = basis.ORtoHO.t() * U(NEUTRON);
  Utot(PROTON ) = basis.ORtoHO.t() * U(PROTON );
  Vtot(NEUTRON) = basis.ORtoHO.t() * V(NEUTRON);
  Vtot(PROTON ) = basis.ORtoHO.t() * V(PROTON );

  for (INT iso: {NEUTRON, PROTON})
  {
    INT ic = 0;

    for (auto &lambda : collectiveCoordinates)
    {
      Qtot(iso, ic) = Utot(iso).t() * basis.HOtoOR * multipoleOperators.ql0(lambda) * basis.HOtoOR.t() * Vtot(iso)
                    + Vtot(iso).t() * basis.HOtoOR * multipoleOperators.ql0(lambda) * basis.HOtoOR.t() * Utot(iso);
      ic++;
    }
  }

  Multi<arma::mat> M;

  for (INT iso: {NEUTRON, PROTON})
  {
    arma::mat eitot = iso == NEUTRON ? Tools::matFromCol(eneQP(NEUTRON), eneQP(NEUTRON).n_elem) : Tools::matFromCol(eneQP(PROTON), eneQP(PROTON).n_elem);
    arma::mat tempMat = 1.0 / (eitot + eitot.t());

    for (INT k = 0; k < 4; k++)
    {
      M(iso, k) = arma::zeros(collectiveCoordinates.n_elem, collectiveCoordinates.n_elem);

      for (INT i = 0; i < collectiveCoordinates.n_elem; i++)
      {
        for (INT j = 0; j < collectiveCoordinates.n_elem; j++)
        {
          M(iso, k)(i, j) = 2.0 * arma::accu(Qtot(iso, i) % Qtot(iso, j) % arma::pow(tempMat, double(k)));
        }
      }
      // Tools::info(PF("M(%s, %d)", Tools::strIsospin(iso).c_str(), k), M(iso, k), true);
    }

  }

  // for (INT k = 0; k < 4; k++)
  // {
  //   arma::mat toto = M(NEUTRON, k) + M(PROTON, k);
  //   Tools::info(PF("k=%d", k), toto, true);
  // }

  arma::mat M1 = M(NEUTRON, 1) + M(PROTON, 1);
  arma::mat M2 = M(NEUTRON, 2) + M(PROTON, 2);
  arma::mat M3 = M(NEUTRON, 3) + M(PROTON, 3);

  arma::mat M1i = M1.i();

  // metric calculation
  metric = 0.5 * M1i * M2 * M1i;

  //===== GCM prescription =====
  massGCM = 4.0 * metric * M1 * metric;
  zpeGCM = 0.5 * massGCM.i() * metric;

  //===== ATDHF prescription =====
  massATDHF = M1i * M3 * M1i;
  zpeATDHF = 0.5 * massATDHF.i() * metric;

  // Tools::info("M1", M1, true);
  // Tools::info("M1.i", M1i, true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Check if rho and kappa matrices are filled, to identify an instance
 *  constructed from a pure input DataTree.
 */

bool State::checkSolution(void) const
{
  DBG_ENTER;

  bool value = (!rho(NEUTRON).empty()) && (!rho(PROTON ).empty()) && (!kappa(NEUTRON).empty()) && (!kappa(PROTON ).empty());

  if (!value)
  {
    Tools::warning("This State instance does not contain a solution state.");
  }

  DBG_RETURN(value);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print canonical states occupations by omega block.
 */

void State::printCanonicalStatesInfo(const arma::mat &_rho, INT iso)
{
  DBG_ENTER;

  std::string state;

  // Dependencies
  basis.calcWDN();

  // Useful quantities.
  arma::mat Mt = basis.ORtoHO.t();

  for (INT omega = 0 ; omega < basis.mMax ; omega++)
  {
    arma::mat Mt_omega = Mt.submat(basis.omegaIndexOR(omega), basis.omegaIndexHO(omega));

    arma::mat rho_omega = _rho.submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
    arma::mat D_omega;
    arma::vec v_square_omega;

    arma::eig_sym(v_square_omega, D_omega, rho_omega);
    v_square_omega = arma::abs(v_square_omega);

    INT m = basis.ORqn(0, basis.omegaIndexHO(omega)(0));
    INT s = basis.ORqn(4, basis.omegaIndexHO(omega)(0));
    std::string label = PF("%d/2", 2 * m - 2 * s + 1);

    Tools::info(PF("Can. %s ", iso == NEUTRON ? "n" : "p") + label, v_square_omega, true);
  } // omega

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the occupation numbers in the canonical representation.
 */

const arma::vec State::getCanonicalV2(const arma::mat &_rho)
{
  DBG_ENTER;

  // Dependencies
  basis.calcWDN();

  // Useful quantities.
  Multi<arma::mat> U_omega;
  Multi<arma::mat> V_omega;
  arma::mat Mt = basis.ORtoHO.t();

  arma::vec state = arma::zeros(basis.ORqn.nb);

  for (INT omega = 0 ; omega < basis.mMax ; omega++)
  {
    arma::mat Mt_omega = Mt.submat(basis.omegaIndexOR(omega), basis.omegaIndexHO(omega));

      arma::mat rho_omega = _rho.submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
      arma::mat D_omega;
      arma::vec v_square_omega;

      arma::eig_sym(v_square_omega, D_omega, rho_omega);
      v_square_omega = arma::abs(v_square_omega);
      state.elem(basis.omegaIndexOR(omega)) = v_square_omega;
  } // omega

  DBG_RETURN(state);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate u_c, v_c and D_c in the OR basis.
 */
void State::calcCanonical(void)
{
  DBG_ENTER;

  // Dependencies
  basis.calcWDN();

  // Useful quantities.
  arma::mat Mt = basis.ORtoHO.t();

  // Initialization.
  v_c(PROTON ) = arma::zeros(basis.ORqn.nb);
  v_c(NEUTRON) = arma::zeros(basis.ORqn.nb);
  u_c(PROTON ) = arma::zeros(basis.ORqn.nb);
  u_c(NEUTRON) = arma::zeros(basis.ORqn.nb);
  D_c(PROTON ) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);
  D_c(NEUTRON) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);

  for (int omega = 0 ; omega < basis.mMax ; omega++)
  {
    arma::mat Mt_omega = Mt.submat(basis.omegaIndexOR(omega), basis.omegaIndexHO(omega));

    for (int iso: {NEUTRON, PROTON})
    {
      // Compute rho and kappa for this omega block in the OR basis
      arma::mat rho_omega   = Mt_omega * rho(iso).submat(basis.omegaIndexHO(omega),   basis.omegaIndexHO(omega)) * Mt_omega.t();
      arma::mat kappa_omega = Mt_omega * kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * Mt_omega.t();

      // Diagonalize rho
      arma::mat D_omega;
      arma::vec v_square_omega;
      arma::eig_sym(v_square_omega, D_omega, rho_omega);

      // Clamp to avoid small negative eigenvalues
      double epsilon = 1e-28;
      v_square_omega = arma::clamp(v_square_omega, epsilon, 1.-epsilon);

      // Extract the sign of uv.
      arma::vec uv_omega = arma::diagvec(D_omega.t() * kappa_omega * D_omega);
      arma::vec sign = arma::sign(uv_omega);
      sign = sign + arma::ones(sign.n_rows) - sign%sign; // set sign to 1 when the value is 0

      // Store u_c with the convention u_c(i) >= 0 in OR representation.
      u_c(iso).elem(basis.omegaIndexOR(omega)) = arma::sqrt( 1. - v_square_omega );
      v_c(iso).elem(basis.omegaIndexOR(omega)) = sign % arma::sqrt(v_square_omega);
      D_c(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = D_omega;
    } // iso
  } // omega

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the omega contributions to rho matrices.
 */

const std::string State::getOmegaContributionsInfo(void) const
{
  DBG_ENTER;

  Multi<arma::vec> contribs;

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT type = 0; type < 4; type++)
    {
      contribs(iso, type) = arma::zeros(basis.mMax);

      for (INT omega = 0; omega < basis.mMax; omega++)
      {
        switch(type)
        {
          case 0:
            if (!rho(iso).empty())
              contribs(iso, type)(omega) = arma::norm(rho(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)), "inf");
            break;
          case 1:
            if (!kappa(iso).empty())
              contribs(iso, type)(omega) = arma::norm(kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)), "inf");
            break;
          case 2:
            if (!U(iso).empty())
              contribs(iso, type)(omega) = arma::norm(U(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)), "inf");
            break;
          case 3:
            if (!V(iso).empty())
              contribs(iso, type)(omega) = arma::norm(V(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)), "inf");
            break;
        }
      }
    }
  }

  std::string state = "";

  state += Tools::treeStr(
      {
      {"OmegaContributions", ""},
      {"rho  n", Tools::vecToStr(contribs(NEUTRON, 0))},
      {"rho  p", Tools::vecToStr(contribs(PROTON , 0))},
      {"kappan", Tools::vecToStr(contribs(NEUTRON, 1))},
      {"kappap", Tools::vecToStr(contribs(PROTON , 1))},
      {"Un    ", Tools::vecToStr(contribs(NEUTRON, 2))},
      {"Up    ", Tools::vecToStr(contribs(PROTON , 2))},
      {"Vn    ", Tools::vecToStr(contribs(NEUTRON, 3))},
      {"Vp    ", Tools::vecToStr(contribs(PROTON , 3))},
      }, false);

  DBG_RETURN(state);
}
