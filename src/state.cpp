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
#include "states.h"

/** \file
 *  \brief Methods of the State class.
 */

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

  rho(NEUTRON)       = arma::mat();
  rho(PROTON )       = arma::mat();
  kappa(NEUTRON)     = arma::mat();
  kappa(PROTON )     = arma::mat();
  U(NEUTRON)         = arma::mat();
  U(PROTON )         = arma::mat();
  V(NEUTRON)         = arma::mat();
  V(PROTON )         = arma::mat();

  // // Construct U and V from the dataTree instance
  // dataTree.get(U, "state/U", true);
  // dataTree.get(V, "state/V", true);
  // calcRhoKappaFromUV();

  // Try to construct rho and kappa from the dataTree instance
  Multi<arma::mat> multiRho;
  Multi<arma::mat> multiKappa;
  dataTree.get(multiRho  , "state/rho"  , true);
  dataTree.get(multiKappa, "state/kappa", true);

  if ((multiRho.size() > 2) && (multiKappa.size() > 2))
  {
    Basis stateBasis = Basis(dataTree, "state/");

    for (INT iso: {NEUTRON, PROTON})
    {
      rho(iso) = arma::zeros(stateBasis.HOqn.nb, stateBasis.HOqn.nb);
      kappa(iso) = arma::zeros(stateBasis.HOqn.nb, stateBasis.HOqn.nb);

      for (INT omega = 0; omega < stateBasis.mMax; omega++)
      {
        // INFO("%d %d", iso, omega);
        // Tools::info("multiRho", multiRho(iso, omega));
        // Tools::info("multiKappa", multiKappa(iso, omega));

        rho(  iso).submat(stateBasis.omegaIndexHO(omega), stateBasis.omegaIndexHO(omega)) = multiRho(iso, omega);
        kappa(iso).submat(stateBasis.omegaIndexHO(omega), stateBasis.omegaIndexHO(omega)) = multiKappa(iso, omega);
      }
    }
  }
  else if ((multiRho.size() > 0) && (multiKappa.size() > 0))
  {
    // try to load using previous storage format. TODO: to be removed
    dataTree.get(rho,             "state/rho", true);
    dataTree.get(kappa,           "state/kappa", true);
  }

  dataTree.get(totalEnergy,       "state/totalEnergy", true);
  dataTree.get(converged,         "state/converged", true);
  dataTree.get(nbIter,            "state/nbIter", true);
  dataTree.get(calculationLength, "state/calculationLength", true);

  //============================================================================

  // read blocked QP states
  IVEC blockedQPsOmega;
  dataTree.get(blockedQPsOmega,   "state/blockedQPsOmega", true);
  IVEC blockedQPsIndex;
  dataTree.get(blockedQPsIndex,   "state/blockedQPsIndex", true);
  IVEC blockedQPsIsospin;
  dataTree.get(blockedQPsIsospin, "state/blockedQPsIsospin", true);

  ASSERT((blockedQPsOmega.n_elem == blockedQPsIndex.n_elem) && (blockedQPsOmega.n_elem == blockedQPsIsospin.n_elem), "Sizes mismatch for blocked QP states.");
  blockedQPStates.clear();

  for (INT i = 0; i < blockedQPsOmega.n_elem; i++)
  {
    StateId is = {blockedQPsIndex(i), blockedQPsOmega(i), blockedQPsIsospin(i)};
    blockedQPStates.insert(is);
  }

  //============================================================================

  // read chemical potentials
  dataTree.get(chemPot, "state/chemicalPotential", true);

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
        // basis.calcWDN();

        Tools::debug("The basis specified in 'state/basis/' will be used.");
        Tools::debug("To use the basis specified in 'basis/', set `basis/useStateBasis` to False.");
      }
      else
      {
        Tools::debug("Converting the state from the basis in `state/basis/` to the basis in `basis/`.");

        // Convert the state FROM the stateBasis to the current basis
        convertFrom(stateBasis);
      }
    }
  }

  //============================================================================

  basis.calcWDN();

  //============================================================================

  // read QP states
  Multi<VEC> qpStatesEnergy;
  Multi<IVEC> qpStatesIndex;
  Multi<VEC> qpStatesOccupation;
  dataTree.get(qpStatesEnergy,     "state/qpStatesEnergy", true);
  dataTree.get(qpStatesIndex,      "state/qpStatesIndex", true);
  dataTree.get(qpStatesOccupation, "state/qpStatesOccupation", true);

  for (auto iso: {NEUTRON, PROTON})
  {
    qpStates(iso) = States(PF("QP. %s", iso == NEUTRON ? "n" : "p"));
  }

  bool qpOK = true;
  INT nbQP = 0;

  if (qpStatesEnergy.size() > 0)
  {
    nbQP = basis.ORqn.nb;

    if (qpStatesEnergy(    NEUTRON).n_elem != nbQP) { Tools::warning("wrong size of qpStatesEnergy(n)"    ); qpOK = false;}
    if (qpStatesEnergy(    PROTON ).n_elem != nbQP) { Tools::warning("wrong size of qpStatesEnergy(p)"    ); qpOK = false;}
    if (qpStatesIndex(     NEUTRON).n_elem != nbQP) { Tools::warning("wrong size of qpStatesIndex(n)"     ); qpOK = false;}
    if (qpStatesIndex(     PROTON ).n_elem != nbQP) { Tools::warning("wrong size of qpStatesIndex(p)"     ); qpOK = false;}
    if (qpStatesOccupation(NEUTRON).n_elem != nbQP) { Tools::warning("wrong size of qpStatesOccupation(n)"); qpOK = false;}
    if (qpStatesOccupation(PROTON ).n_elem != nbQP) { Tools::warning("wrong size of qpStatesOccupation(p)"); qpOK = false;}
  }
  else
  {
    qpOK = false;
  }

  if (qpOK)
  {
    for (INT i = 0; i < nbQP; i++)
    {
      UINT index = qpStatesIndex(NEUTRON)(i);

      INT m = basis.ORqn(0, index);
      INT s = basis.ORqn(4, index);

      qpStates(NEUTRON).add(qpStatesIndex(NEUTRON)(i),
                            qpStatesEnergy(NEUTRON)(i),
                            qpStatesOccupation(NEUTRON)(i),
                            {basis.blockIdOR(index), m - s, NEUTRON});
      qpStates(PROTON ).add(qpStatesIndex(PROTON)(i),
                            qpStatesEnergy(PROTON)(i),
                            qpStatesOccupation(PROTON)(i),
                            {basis.blockIdOR(index), m - s, PROTON});
    }
    Tools::debug(PF("%d QP states loaded", nbQP));
  }
  else
  {
    Tools::debug("QP states not loaded");
  }

  //============================================================================

  // Construct constraints
  constraints = Constraint::fromDataTree(dataTree);

  // If missing rho and kappa and U and V present, reconstruct from U and V
  if (rho(NEUTRON).empty() || rho(PROTON).empty() ||
      kappa(NEUTRON).empty() || kappa(PROTON).empty())
  {
    calcRhoKappaFromUV();
  }

  // If missing U and V and rho and kappa present, reconstruct from rho and kappa
  if (U(NEUTRON).empty() || U(PROTON).empty() ||
      V(NEUTRON).empty() || V(PROTON).empty())
  {
    calcUVFromRhoKappa();
  }

  // Tools::mesg("BogCnv", getOmegaContributionsInfo());

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Save the state in a DataTree instance.
 */

DataTree State::getDataTree(void)
{
  DBG_ENTER;

  DataTree dt;

  // State datatree must contain the associated basis
  dt.merge(basis.getDataTree("state/"));

  // State datatree must contain the associated system
  dt.merge(sys.getDataTree());

  // Save rho and kappa matrices
  Multi<arma::mat> multiRho;
  Multi<arma::mat> multiKappa;
  for (INT iso: {NEUTRON, PROTON})
  {
    if (rho(iso).empty()) break; // no rho or kappa matrices to save

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

  dt.set("state/rho",                         multiRho);
  dt.set("state/kappa",                       multiKappa);
  // dt.set("state/U",                           U);
  // dt.set("state/V",                           V);

  // TODO: save individual states

  IVEC blockedQPsOmega   = arma::zeros<IVEC>(blockedQPStates.size());
  IVEC blockedQPsIndex   = arma::zeros<IVEC>(blockedQPStates.size());
  IVEC blockedQPsIsospin = arma::zeros<IVEC>(blockedQPStates.size());

  INT i = 0;
  for (auto &bsid: blockedQPStates)
  {
    blockedQPsOmega(i) = bsid.omega;
    blockedQPsIndex(i) = bsid.index;
    blockedQPsIsospin(i) = bsid.isospin;
    i++;
  }
  dt.set("state/blockedQPsOmega", blockedQPsOmega);
  dt.set("state/blockedQPsIndex", blockedQPsIndex);
  dt.set("state/blockedQPsIsospin", blockedQPsIsospin);

  Multi<VEC> qpStatesEnergy;
  qpStatesEnergy(NEUTRON) = qpStates(NEUTRON).getEnergy();
  qpStatesEnergy(PROTON ) = qpStates(PROTON ).getEnergy();

  Multi<IVEC> qpStatesIndex;
  qpStatesIndex(NEUTRON) = qpStates(NEUTRON).getIndex();
  qpStatesIndex(PROTON ) = qpStates(PROTON ).getIndex();

  Multi<VEC> qpStatesOccupation;
  qpStatesOccupation(NEUTRON) = qpStates(NEUTRON).getOccupation();
  qpStatesOccupation(PROTON ) = qpStates(PROTON ).getOccupation();

  dt.set("state/qpStatesEnergy",     qpStatesEnergy);
  dt.set("state/qpStatesIndex",      qpStatesIndex);
  dt.set("state/qpStatesOccupation", qpStatesOccupation);
  dt.set("state/chemicalPotential",  chemPot);
  dt.set("state/totalEnergy",        totalEnergy);
  dt.set("state/converged",          converged);
  dt.set("state/nbIter",             nbIter);
  dt.set("state/calculationLength",  calculationLength);

  // set the energy constraints
  for (auto &c : constraints)
  {
    dt.merge(c.second.getDataTree());
  }

  DBG_RETURN(dt);
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
    // Tools::mesg("State ", PF_YELLOW("Identical bases => no conversion"));

    DBG_LEAVE;
  }

  if (empty())
  {
    // Tools::mesg("State ", PF_YELLOW("empty state => no conversion"));

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
  Tools::mesg("State ", "State basis conversion:");
  Tools::mesg("State ", "from: " + initialBasis.info(true));
  Tools::mesg("State ", "to:   " + basis.info(true));

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

void State::calcUVFromRhoKappa(void)
{
  DBG_ENTER;

  if (rho(NEUTRON  ).empty()) { DBG_LEAVE; };
  if (rho(PROTON   ).empty()) { DBG_LEAVE; };
  if (kappa(NEUTRON).empty()) { DBG_LEAVE; };
  if (kappa(PROTON ).empty()) { DBG_LEAVE; };

  Tools::debug("Updating U and V matrices from rho and kappa matrices");
  // Tools::mesg("State ", info());

  //============================================================================
  //============================================================================
  //============================================================================

  // basis orthonormalization
  basis.calcWDN();

  arma::mat M = basis.HOtoOR.t();
  arma::mat Mi = basis.ORtoHO.t();

  for (INT iso: {NEUTRON, PROTON})
  {
    U(iso) = arma::zeros(basis.HOqn.nb, basis.ORqn.nb);
    V(iso) = arma::zeros(basis.HOqn.nb, basis.ORqn.nb);

    for (INT omega = 0; omega < basis.mMax; omega++)
    {
      arma::mat M_omega  = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
      arma::mat Mi_omega  = Mi.submat(basis.omegaIndexOR(omega), basis.omegaIndexHO(omega));

      UINT dimOR = M_omega.n_cols;

      // rho matrix in OR representation
      arma::mat rhoOR   = Mi_omega * rho(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * Mi_omega.t();

      // kappa matrix in OR representation
      arma::mat kappaOR = Mi_omega * kappa(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * Mi_omega.t();

      // construct R matrix in 2*OR representation
      arma::mat R = arma::zeros(dimOR * 2, dimOR * 2);
      R.submat(0, 0, dimOR - 1, dimOR - 1) = rhoOR;
      R.submat(dimOR, dimOR, 2 * dimOR - 1, 2 * dimOR - 1) = arma::eye(dimOR, dimOR) - rhoOR;
      R.submat(0, dimOR, dimOR - 1, 2 * dimOR - 1) = kappaOR;
      R.submat(dimOR, 0, 2 * dimOR - 1, dimOR - 1) = kappaOR;

      // diagonalize
      arma::vec eige;
      arma::mat eigv;
      Tools::eig_sym(eige, eigv, R);

      // sort by decreasing eigenvalues.
      UVEC Eind = arma::sort_index(eige, "descend");
      eigv = eigv.cols(Eind);
      eige = eige.elem(Eind);

      // Tools::info("eige", eige, true);

      // extraction of U and V matrices
      arma::mat U_omega = eigv.submat(dimOR, 0, 2 * dimOR - 1, dimOR - 1);
      arma::mat V_omega = eigv.submat(0, 0, dimOR - 1, dimOR - 1);

      // back to HO representation
      arma::mat Uho = M_omega * U_omega;
      arma::mat Vho = M_omega * V_omega;

      U(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)) = Uho;
      V(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)) = Vho;
    }
  }

  // Tools::mesg("State ", "Updating U and V matrices from rho and kappa matrices - after:");
  // Tools::mesg("State ", info());

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

  Tools::mesg("State ", "Updating rho and kappa matrices from U and V matrices");
  // Tools::mesg("State ", info());

  rho(NEUTRON)   = V(NEUTRON) * V(NEUTRON).t();
  rho(PROTON )   = V(PROTON ) * V(PROTON ).t();
  kappa(NEUTRON) = V(NEUTRON) * U(NEUTRON).t();
  kappa(PROTON ) = V(PROTON ) * U(PROTON ).t();

  // rho(NEUTRON)   = arma::symmatu(rho(NEUTRON));
  // rho(PROTON )   = arma::symmatu(rho(PROTON ));
  // kappa(NEUTRON) = arma::symmatu(kappa(NEUTRON));
  // kappa(PROTON ) = arma::symmatu(kappa(PROTON ));

  // Tools::mesg("State ", "Updating rho and kappa matrices from U and V matrices - after:");
  // Tools::mesg("State ", info());

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

  Discrete d1(basis);
  Discrete d2(otherState.basis);

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
    ASSERT(qpStates.contains(NEUTRON), "no NEUTRON key in state.qpStates");
    ASSERT(qpStates.contains(PROTON ), "no PROTON  key in state.qpStates");
    ASSERT(rho.contains(NEUTRON), "no NEUTRON key in state.rho");
    ASSERT(rho.contains(PROTON ), "no PROTON  key in state.rho");
    ASSERT(chemPot.n_elem == 2, "wrong dimension for state.chemPot");
    ASSERT(U.contains(NEUTRON), "no NEUTRON key in state.U");
    ASSERT(U.contains(PROTON ), "no PROTON  key in state.U");
    ASSERT(V.contains(NEUTRON), "no NEUTRON key in state.V");
    ASSERT(V.contains(PROTON ), "no PROTON  key in state.V");
    ASSERT(kappa.contains(NEUTRON), "no NEUTRON key in state.kappa");
    ASSERT(kappa.contains(PROTON ), "no PROTON  key in state.kappa");

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
      {"block.", Tools::infoStr(blockedQPStates)},
      {"eneTot", Tools::infoStr(totalEnergy)},
      {"chemPn", Tools::infoStr(chemPot(NEUTRON))},
      {"chemPp", Tools::infoStr(chemPot(PROTON ))},
      {"QPs n.", PF("%d states", qpStates(NEUTRON).n_elem)},
      {"QPs p.", PF("%d states", qpStates(PROTON ).n_elem)},
      // {"ZpeGcm", Tools::infoStr(arma::trace(zpeGCM))},
      // {"ZpeAtd", Tools::infoStr(arma::trace(zpeATDHF))},
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

  std::string result = "";

  std::list<std::string> names;
  std::list<std::string> units;
  std::list<std::list<std::string> > valuesMetric;
  std::list<std::list<std::string> > valuesMassGCM;
  std::list<std::list<std::string> > valuesZpeGCM;
  std::list<std::list<std::string> > valuesMassATDHF;
  std::list<std::list<std::string> > valuesZpeATDHF;

  if (what == "inertia")
  {
    if (!metric.empty())
    {
      INT ic = 0;
      for (auto &p : collectiveCoordinates)
      {
        INT lambda = p.first;
        INT mu = p.second;

        names.push_back(PF("q%01d%01d", lambda, mu));
        units.push_back(PF("[fm%01d]", lambda));

        std::list<std::string> val;
        val.push_back(PF("q%01d%01d [fm%01d]", lambda, mu, lambda));
        for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", metric(ic, ic2)));
        valuesMetric.push_back(val);

        val.clear();
        val.push_back(PF("q%01d%01d [fm%01d]", lambda, mu, lambda));
        for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", massGCM(ic, ic2)));
        valuesMassGCM.push_back(val);

        val.clear();
        val.push_back(PF("q%01d%01d [fm%01d]", lambda, mu, lambda));
        for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", zpeGCM(ic, ic2)));
        valuesZpeGCM.push_back(val);

        val.clear();
        val.push_back(PF("q%01d%01d [fm%01d]", lambda, mu, lambda));
        for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", massATDHF(ic, ic2)));
        valuesMassATDHF.push_back(val);

        val.clear();
        val.push_back(PF("q%01d%01d [fm%01d]", lambda, mu, lambda));
        for (INT ic2 = 0; ic2 < collectiveCoordinates.size(); ic2++) val.push_back(PF("%13.6e", zpeATDHF(ic, ic2)));
        valuesZpeATDHF.push_back(val);

        ic++;
      }

      result += Tools::valueTable("Metric"         , names, units, valuesMetric   ) + "\n";
      result += Tools::valueTable("Mass GCM"       , names, units, valuesMassGCM  ) + "\n";
      result += Tools::valueTable("ZPE GCM [MeV]"  , names, units, valuesZpeGCM   ) + "\n";
      result += Tools::valueTable("Mass ATDHF"     , names, units, valuesMassATDHF) + "\n";
      result += Tools::valueTable("ZPE ATDHF [MeV]", names, units, valuesZpeATDHF );
    }
    else
    {
      result += PF_YELLOW("Inertia have not been calculated.");
    }
  }
  else
  {
    result += info();
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the inertia tensor, the metric and the ZPE energy (ATDHF and GCM prescriptions).
 */

void State::calcInertia(const std::string &interactionName, const std::list<std::pair<INT, INT> > &_collectiveCoordinates)
{
  DBG_ENTER;

  if (qpStates(NEUTRON).n_elem == 0)
  {
    Tools::mesg("State ", "(Re)calculating the QP states and U/V matrices => performing iterations of the Broyden solver");
    DataTree dt = DataTree::getDefault() + getDataTree();
    dt.set("interaction/name", interactionName);

    SolverHFBBroyden solver(dt, *this);
    solver.init();

    Tools::mesg("State ", solver.info());

    while(solver.nextIter());
    *this = solver.state;
  }

  qpStates(NEUTRON).sort("energy");
  qpStates(PROTON ).sort("energy");

  // INFO("before calculating inertia...");
  // INFO(qpStates(NEUTRON).info(-1, true));
  // INFO(qpStates(PROTON ).info(-1, true));

  ASSERT(!U(NEUTRON).empty(), "Empty U(NEUTRON)");
  ASSERT(!U(PROTON ).empty(), "Empty U(PROTON )");
  ASSERT(!V(NEUTRON).empty(), "Empty V(NEUTRON)");
  ASSERT(!V(PROTON ).empty(), "Empty V(PROTON )");
  ASSERT(!qpStates.empty()  , "No QP states");

  if (!_collectiveCoordinates.empty()) collectiveCoordinates = _collectiveCoordinates;

  // ===========================================================================
  // ===========================================================================
  // ===========================================================================

  // Construct and print the list of collective coordinates
  std::string listStr;
  INT i = 0;
  for (auto p: collectiveCoordinates)
  {
    INT lambda = p.first;
    INT mu = p.second;

    listStr += PF_YELLOW("Q(%01d,%01d)", lambda, mu);
    i++;
    if (i != collectiveCoordinates.size()) listStr += ", ";
  }
  Tools::mesg("State ", "Collective coordinates for the inertia calculation: (" + listStr + ")");

  // ===========================================================================
  // ===========================================================================
  // ===========================================================================

  // Initialize Qtot matrices
  Multi<arma::mat> Qtot;

  for (auto iso: {NEUTRON, PROTON})
  {
    for (auto &p: collectiveCoordinates)
    {
      INT lambda = p.first;
      INT mu = p.second;

      Qtot(iso, lambda, mu) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);
    }
  }

  // ===========================================================================

  // Fill Qtot matrices
  MultipoleOperators multipoleOperators(*this);
  multipoleOperators.calcQlmHO();

  for (INT iso: {NEUTRON, PROTON})
  {
    for (auto &p: collectiveCoordinates)
    {
      INT lambda = p.first;
      INT mu = p.second;

      ASSERT(multipoleOperators.qlmHO.contains(lambda, mu),
             PF("HO representation of multipole operator missing for lambda=%d and mu=%d.", lambda, mu));

      arma::mat qlmHO = multipoleOperators.qlmHO(lambda, mu);

      if (mu != 0)
      {
        // TODO: understand why this factor is needed to reproduce Berger2ct values
        qlmHO = sqrt(sqrt(2.0 * PI)) * ( multipoleOperators.qlmHO(lambda,  mu)
                                       + multipoleOperators.qlmHO(lambda, -mu));
      }

      Qtot(iso, lambda, mu) = U(iso).t() * qlmHO * V(iso)
                            + V(iso).t() * qlmHO * U(iso);
    }
  }

  // ===========================================================================

  // INFO(info());
  // Tools::info("U(NEUTRON)", U(NEUTRON), true);
  // Tools::info("Qtot(NEUTRON, 2, 0)", Qtot(NEUTRON, 2, 0));
  // Tools::end();

  // ===========================================================================

  INT nbDirections = collectiveCoordinates.size();

  Multi<arma::mat> M;

  Multi<VEC> eneQP;
  eneQP(NEUTRON) = qpStates(NEUTRON).getEnergy();
  eneQP(PROTON ) = qpStates(PROTON ).getEnergy();

  for (INT iso: {NEUTRON, PROTON})
  {
    arma::mat eitot = Tools::matFromCol(eneQP(iso), eneQP(iso).n_elem);
    arma::mat tempMat = 1.0 / (eitot + eitot.t());

    for (INT order = 0; order < 4; order++)
    {
      M(iso, order) = arma::zeros(nbDirections, nbDirections);

      INT i = 0;
      for (auto pi: collectiveCoordinates)
      {
        INT lambdai = pi.first;
        INT mui = pi.second;

        INT j = 0;
        for (auto pj: collectiveCoordinates)
        {
          INT lambdaj = pj.first;
          INT muj = pj.second;

          M(iso, order)(i, j) = 2.0 * arma::accu(Qtot(iso, lambdai, mui) % Qtot(iso, lambdaj, muj) % arma::pow(tempMat, double(order)));
          j++;
        }
        i++;
      }
    }
  }

  // DEBUG
  // for (INT order = 0; order < 4; order++)
  // {
  //   arma::mat toto = M(NEUTRON, order) + M(PROTON, order);
  //   Tools::info(PF("order=%d", order), toto, true);
  // }
  // DEBUG

  arma::mat M1 = M(NEUTRON, 1) + M(PROTON, 1);
  arma::mat M2 = M(NEUTRON, 2) + M(PROTON, 2);
  arma::mat M3 = M(NEUTRON, 3) + M(PROTON, 3);

  arma::mat M1i;
  try
  {
    M1i = M1.i();
  }
  catch (...)
  {
    Tools::info("M1", M1, true);
    Tools::warning("M1 is a singular matrix");
    DBG_LEAVE;
  }

  // metric calculation
  metric = 0.5 * M1i * M2 * M1i;

  //===== GCM prescription =====
  //    Eq. (95) in N. Schunck and L. Robledo, Rep. Prog. Phys. 79 (2016) 116301
  //    Eq. (45) in R. Navarro Perez et al, Comp. Phys. Comm. 220, (2017) 263
  massGCM = 4.0 * metric * M1 * metric;
  //    Eq. (82) and Eq. (96) in N. Schunck and L. Robledo, Rep. Prog. Phys. 79 (2016) 116301
  // != Eq. (41) and Eq. (46) in R. Navarro Perez et al, Comp. Phys. Comm. 220, (2017) 263

  try
  {
    zpeGCM = 0.5 * massGCM.i() * metric;
  }
  catch (...)
  {
    Tools::warning("massGCM is a singular matrix");
    DBG_LEAVE;
  }

  //===== ATDHF prescription =====
  // != Eq. (110) in N. Schunck and L. Robledo, Rep. Prog. Phys. 79 (2016) 116301
  // != Eq. (48)  in R. Navarro Perez et al, Comp. Phys. Comm. 220, (2017) 263
  massATDHF = M1i * M3 * M1i;

  try
  {
    // May or may not be used. See discussion p. 28-29 in N. Schunck and L. Robledo, Rep. Prog. Phys. 79 (2016) 116301
    zpeATDHF = 0.5 * massATDHF.i() * metric;
  }
  catch (...)
  {
    Tools::warning("massATDHF is a singular matrix");
    DBG_LEAVE;
  }

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
