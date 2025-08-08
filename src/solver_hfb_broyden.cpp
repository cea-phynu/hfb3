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

#include "solver_hfb_broyden.h"
#include "tools.h"
#include "plot.h"
#include "mixing.h"

/** \file
 *  \brief Methods of the SolverHFBBroyden class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > SolverHFBBroyden::validKeys =
  {
    { "solver/broyden/cvgTarget"           , "Convergence target value."                                         , "1e-6" , "D" },
    { "solver/broyden/maxIter"             , "Maximum number of iterations"                                      , "200"  , "I" },
    { "solver/broyden/emptyKappaProtection", "Force non-empty kappa matrices"                                    , "False", "B" },
    { "solver/broyden/cvgTargetLambda"     , "Convergence target value for lambda-iterations"                    , "1e-5" , "D" },
    { "solver/broyden/maxIterLambda"       , "Maximum number of iterations for lambda-iterations"                , "20"   , "I" },
    { "solver/broyden/earlyLambdaMixing"   , "Mix Lagrange multipliers in the linear phase of the Broyden mixing", "True" , "B" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

SolverHFBBroyden::SolverHFBBroyden(const std::string &filename) : SolverHFBBroyden(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

SolverHFBBroyden::SolverHFBBroyden(const DataTree &dataTree) : SolverHFBBroyden(dataTree, State(dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

SolverHFBBroyden::SolverHFBBroyden(const DataTree &dataTree, State _state) : Solver(dataTree, _state),
  mixing(dataTree),
  interaction(dataTree, &state),
  discrete(&state.basis,  Mesh::regular(-10.0, 0.0, -15.0, 10.0, 0.0, 15.0, 101, 1, 151)),
  multipoleOperators(state)
{
  DBG_ENTER;

  dataTree.get(maxIter,              "solver/broyden/maxIter",              true);
  dataTree.get(emptyKappaProtection, "solver/broyden/emptyKappaProtection", true);
  dataTree.get(cvgTarget,            "solver/broyden/cvgTarget",            true);
  dataTree.get(lambdaMax,            "solver/broyden/lambdaMax",            true);
  dataTree.get(lambdaIterMax,        "solver/broyden/lambdaIterMax",        true);

  INT earlyLambdaMixingInt = 1;
  dataTree.get(earlyLambdaMixingInt, "solver/broyden/earlyLambdaMixing",    true);

  earlyLambdaMixing = (earlyLambdaMixingInt != 0);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the HFB loop.
 */

void SolverHFBBroyden::init()
{
  DBG_ENTER;

  startTime = Tools::clock();

  nbIter = 0;

  if (maxIter < 1) maxIter = 1;

  // Chemical potentials (lagrange multipliers for the numbers of particles).
  // state.chemPot = arma::zeros(2);

  // Initialize objects.
  rhoNext(NEUTRON)   = state.rho(NEUTRON);
  rhoNext(PROTON )   = state.rho(PROTON );
  kappaNext(NEUTRON) = state.kappa(NEUTRON);
  kappaNext(PROTON ) = state.kappa(PROTON );

  // Avoid empty kappa matrices.
  if (arma::norm(state.kappa(NEUTRON), "inf") < 1e-14)
  {
    Tools::mesg("SolBro", "Detected an empty kappa matrix [NEUT]");

    if (emptyKappaProtection == 1)
    {
      Tools::mesg("SolBro", "Forcing non-empty kappa matrix [NEUT]");
      state.kappa(NEUTRON) = state.rho(NEUTRON);
    }
  }

  if (arma::norm(state.kappa(PROTON ), "inf") < 1e-14)
  {
    Tools::mesg("SolBro", "Detected an empty kappa matrix [PROT]");

    if (emptyKappaProtection)
    {
      Tools::mesg("SolBro", "Forcing non-empty kappa matrix [PROT]");
      state.kappa(PROTON ) = state.rho(PROTON );
    }
  }

  ASSERT(!state.rho(NEUTRON  ).empty(), "empty rho(NEUTRON).");
  ASSERT(!state.rho(PROTON   ).empty(), "empty rho(PROTON).");
  ASSERT(!state.kappa(NEUTRON).empty(), "empty kappa(NEUTRON).");
  ASSERT(!state.kappa(PROTON ).empty(), "empty kappa(PROTON).");

  ASSERT(state.constraints.size() != 0, "no constraints !?");

  // Calculate the multipole moment operator matrices.
  multipoleOperators.calcQl0Matrices();

  for (auto &c : state.constraints)
  {
    if (c.second.gender == Constraint::SD || c.second.gender == Constraint::MA || c.second.gender == Constraint::QN)
    {
      fragInLoop = true;
      break;
    }
  }

  if (fragInLoop)
  {
    GeometricalOperators geometricalOperatorsInit(state);
    geometricalOperators = geometricalOperatorsInit;
    geometricalOperators.calcQNeckMatrix();
    geometricalOperators.calcQneck(state.rho(NEUTRON), state.rho(PROTON));
  }

  // Calculate the basis orthogonalization.
  state.basis.calcWDN();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Plot some quantities using Bokeh.
 */

void SolverHFBBroyden::bokehPlot(void)
{
  DBG_ENTER;

  if (useBokeh)
  {
    Plot::slot(0);
    Plot::curve("#iter", "Ene HFB [MeV]", nbIter, ene);

    Plot::slot(1);
    Plot::curve("#iter", "Convergence", nbIter, log10(value));

    if (plotDensities)
    {
      Plot::slot(2);
      arma::mat densn = discrete.getLocalXZ(rhoNext(NEUTRON), true);
      Plot::map("HFB local density (neut)", densn, discrete.mesh);

      Plot::slot(3);
      arma::mat densp = discrete.getLocalXZ(rhoNext(PROTON), true);
      Plot::map("HFB local density (prot)", densp, discrete.mesh);
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Prepare and perform the next HFB iteration.
 */

INT SolverHFBBroyden::nextIter()
{
  DBG_ENTER;

  //============================================================================

// #define PRINT_ITER_DURATION
#ifdef PRINT_ITER_DURATION
  double startTime = Tools::clock();
#endif

  if (!singleHFBiter())
  {
    status = Solver::DIVERGED;
    bestEne = 0.0;

    DBG_RETURN(status);
  }

  bokehPlot();

#ifdef PRINT_ITER_DURATION
  double iterLength = Tools::clock() - startTime;
#endif

  //============================================================================

  // print iteration message

  std::string eneStr = ((ene > -9999.999 ) && (ene < 0.0)) ? PF_GREEN("%9.3f", ene) : PF_RED("%9.2e", ene);
  Tools::mesg("SolBro",
              PF("#it: %03d ", nbIter) +
              PF("cvg: %8.2e ", value) +
              PF("ene: ") + eneStr + " " +
#ifdef PRINT_ITER_DURATION
              PF("tot: %6.3fs (fld: %6.3fs)", iterLength, interaction.calcLength) + " " +
#endif
              PF("ln: %7.3f(#%2d) ", state.chemPot(NEUTRON), lambdaIter(NEUTRON)) +
              PF("lp: %7.3f(#%2d)", state.chemPot(PROTON ), lambdaIter(PROTON ))
              // PF_YELLOW(" " + interaction.getWarningStr())
             );
  // Tools::mesg("SolBro", multipoleOperators.getNiceInfo());
  // Tools::mesg("SolBro", state.info());
  // Tools::mesg("SolBro", state.getOmegaContributionsInfo());

  //============================================================================

  nbIter++;
  state.nbIter = nbIter;
  status = Solver::ITERATING;

  if (nbIter >= maxIter)
  {
    Tools::mesg("SolBro", "Maximum number of iterations reached, exiting HFB loop");
    status = Solver::MAXITER;
    DBG_RETURN(status);
  }

  if ((value < cvgTarget) && (nbIter > 1))
  {
    Tools::mesg("SolBro", "Target value reached, exiting HFB loop");
    state.totalEnergy = ene;
    state.converged = true;

    status = Solver::CONVERGED;
    DBG_RETURN(status);
  }

  if (fabs(value) > 1e20)
  {
    Tools::mesg("SolBro", "Convergence is too high, exiting HFB loop");
    status = Solver::DIVERGED;
    DBG_RETURN(status);
  }

  if (fabs(ene) > 1e16)
  {
    Tools::mesg("SolBro", "Energy is too high, exiting HFB loop");
    status = Solver::DIVERGED;
    DBG_RETURN(status);
  }

  state.calculationLength = Tools::clock() - startTime;

  //============================================================================

  // Do the mixing
  mixStates();

  // Re-calculate multipole moments.
  multipoleOperators.calcQlm(state.rho);

  if (fragInLoop)
  {
    if (geometricalOperators.izNeck != -1)
    {
      geometricalOperators.calcQneck(state.rho(NEUTRON), state.rho(PROTON ));
    }
  }

  for (auto &c : state.constraints)
  {
    if (c.second.gender == Constraint::MM)
    {
      c.second.measuredVal = multipoleOperators.qlm(c.second.lm, c.second.iso) * c.second.factor;
    }
    else if (c.second.gender == Constraint::QN)
    {
      if (geometricalOperators.izNeck != -1)
      {
        c.second.measuredVal = geometricalOperators.qneck(c.second.iso) * c.second.factor;
      }
      else
      {
        c.second.measuredVal = c.second.val;
      }
    }
  }

  DBG_RETURN(status);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Mix the previous states using the Broyden Mixing method.
 */

bool SolverHFBBroyden::mixStates(void)
{
  DBG_ENTER;

  // Prepare the Mixing vector.
  UINT nRowsRho   = rhoNext(NEUTRON).n_rows;
  UINT nColsRho   = rhoNext(NEUTRON).n_cols;
  UINT nRowsKappa = kappaNext(NEUTRON).n_rows;
  UINT nColsKappa = kappaNext(NEUTRON).n_cols;

  UINT sizeRho   = nRowsRho * nColsRho;
  UINT sizeKappa = nRowsKappa * nColsKappa;
  UINT sizeConstraints = state.constraints.size();

  UINT index0 = 0;
  UINT index1 = index0 + sizeRho;
  UINT index2 = index1 + sizeRho;
  UINT index3 = index2 + sizeKappa;
  UINT index4 = index3 + sizeKappa;

  UINT totalSize = index4 + sizeConstraints;
  arma::vec vec = arma::zeros(totalSize);

  vec.subvec(index0, index1 - 1) = vectorise(arma::symmatu(rhoNext(NEUTRON)));
  vec.subvec(index1, index2 - 1) = vectorise(arma::symmatu(rhoNext(PROTON )));
  vec.subvec(index2, index3 - 1) = vectorise(arma::symmatu(kappaNext(NEUTRON)));
  vec.subvec(index3, index4 - 1) = vectorise(arma::symmatu(kappaNext(PROTON )));

  UINT i = 0;
  for (auto &c : state.constraints)
  {
    vec(index4 + i) = c.second.lambda;
    i++;
  }

  // Do the mixing !
  if (!mixing.newVec(vec))
  {
    // mixing failed !?
    bestEne = 0.0;

    DBG_RETURN(false);
  }

  // Get the mixed vector.
  vec = mixing.getVec();

  arma::mat oldRhon   = state.rho(NEUTRON);
  arma::mat oldRhop   = state.rho(PROTON );
  arma::mat oldKappan = state.kappa(NEUTRON);
  arma::mat oldKappap = state.kappa(PROTON );

  // Extract the mixed quantities.
  state.rho(NEUTRON) = arma::reshape(vec.subvec(index0, index1 - 1), nRowsRho, nColsRho);
  state.rho(PROTON ) = arma::reshape(vec.subvec(index1, index2 - 1), nRowsRho, nColsRho);
  state.kappa(NEUTRON) = arma::reshape(vec.subvec(index2, index3 - 1), nRowsKappa, nColsKappa);
  state.kappa(PROTON ) = arma::reshape(vec.subvec(index3, index4 - 1), nRowsKappa, nColsKappa);

  /* INFO("variation min rhon: %f"  , (state.rho(NEUTRON)   - oldRhon  ).min()); */
  /* INFO("variation min rhop: %f"  , (state.rho(PROTON )   - oldRhop  ).min()); */
  /* INFO("variation min kappan: %f", (state.kappa(NEUTRON) - oldKappan).min()); */
  /* INFO("variation min kappap: %f", (state.kappa(PROTON ) - oldKappap).min()); */
  /**/
  /* INFO("variation max rhon: %f"  , (state.rho(NEUTRON)   - oldRhon  ).max()); */
  /* INFO("variation max rhop: %f"  , (state.rho(PROTON )   - oldRhop  ).max()); */
  /* INFO("variation max kappan: %f", (state.kappa(NEUTRON) - oldKappan).max()); */
  /* INFO("variation max kappap: %f", (state.kappa(PROTON ) - oldKappap).max()); */

// #define BROYDEEN_SLOWING_FACTOR 0.1
#ifdef BROYDEEN_SLOWING_FACTOR
  double slowingFactor = BROYDEEN_SLOWING_FACTOR;

  state.rho(NEUTRON)   = (1.0 - slowingFactor) * state.rho(NEUTRON)   + slowingFactor * oldRhon;
  state.rho(PROTON )   = (1.0 - slowingFactor) * state.rho(PROTON )   + slowingFactor * oldRhop;
  state.kappa(NEUTRON) = (1.0 - slowingFactor) * state.kappa(NEUTRON) + slowingFactor * oldKappan;
  state.kappa(PROTON ) = (1.0 - slowingFactor) * state.kappa(PROTON ) + slowingFactor * oldKappap;
#endif

  if (mixing.isLinearMode() || earlyLambdaMixing)
  {
    i = 0;
    for (auto &c : state.constraints)
    {
      c.second.lambda = vec(index4 + i); // DEBUG
      // INFO("lambda: %e", c.second.lambda);
      i++;
    }
  }

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Perform the HFB loop.
 */

void SolverHFBBroyden::calc()
{
  DBG_ENTER;

  init();

  //==========================
  //===== HFB LOOP START =====
  //==========================

  while (nextIter() == Solver::ITERATING);

  //========================
  //===== HFB LOOP END =====
  //========================

  if (status == Solver::CONVERGED)
    Tools::mesg("SolBro", PF("iter: %3d, ene: %9.3f MeV", nbIter + 1, bestEne));
  else
    Tools::mesg("SolBro", PF("iter: %3d, NOT CONVERGED", nbIter + 1));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Perform one HFB iteration.
 */

bool SolverHFBBroyden::singleHFBiter(void)
{
  DBG_ENTER;

  //Useful quantities.
  Basis &basis = state.basis;
  arma::mat M = state.basis.HOtoOR.t();

  // Particle numbers.
  arma::vec nPart(2);
  nPart(NEUTRON) = (double)(state.sys.nNeut);
  nPart(PROTON ) = (double)(state.sys.nProt);

  // New rho and kappa matrices initialization.
  rhoNext(NEUTRON)   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  rhoNext(PROTON )   = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappaNext(NEUTRON) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);
  kappaNext(PROTON ) = arma::zeros(basis.HOqn.nb, basis.HOqn.nb);

  // Force some field recalculations.
  interaction.clear();

  // Calculate the energy contributions (and the interaction).
  interaction.calcEnergies();

  // Print the energy contributions.
  // INFO(interaction.getNiceInfo());

  // Print info of a specific field
  // INFO(interaction("central")->info());

  // J.-P. Ebran's derivations, p. 16, Eq. (III-A.9).
  Multi<arma::mat> hTilde;
  Multi<arma::mat> deltaTilde;
  Multi<arma::vec> VV;

  for (INT iso: {NEUTRON, PROTON})
  {
    hTilde(iso) = interaction.getHamiltonianContributions(iso, Field::DIRECT) + interaction.getHamiltonianContributions(iso, Field::EXCHANGE);
    deltaTilde(iso) = interaction.getHamiltonianContributions(iso, Field::PAIRING);

    U(iso) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);
    V(iso) = arma::zeros(basis.ORqn.nb, basis.ORqn.nb);

    qpEne(iso) = arma::zeros(basis.ORqn.nb);
    partEne(iso) = arma::zeros(basis.ORqn.nb);
    VV(iso) = arma::zeros(basis.ORqn.nb);
  }

  // Tools::info("hTilde", hTilde(NEUTRON));

  // for (auto &c : state.constraints)
  // {
  //   Tools::info(c.second.info());
  // }

  //============================================================================
  // lambda-iterations
  // -> adjust the lagrange multipliers for the particle number operators

  // double multipMass;
  // double traceRho;
  // double sumV2;
  // arma::vec v2;

  for (INT iso: {NEUTRON, PROTON})
  {
    // if (state.blockedQP(iso) > -1)
    // {
    //   state.calcCanonicalFromRho(rhoNext);
    //   multipoleOperators.calcQlm(state.rho);
    //   v2 = arma::square(state.v_c(iso));
    //   sumV2 = arma::sum(v2) * 2.0;
    //   traceRho = arma::trace(state.rho(iso)) * 2.0;
    //   multipMass = multipoleOperators.nPart(iso);
    //   INFO("initial sum v2: %9.6f | trace: %9.6f | multip: %9.6f", sumV2, traceRho, multipMass);
    // }

    state.U(iso) = arma::zeros(basis.HOqn.nb, basis.ORqn.nb);
    state.V(iso) = arma::zeros(basis.HOqn.nb, basis.ORqn.nb);

    lambdaIter(iso) = 0;

    arma::vec  hfEnergiesTotal = arma::zeros(state.basis.ORqn.nb);
    IVEC omegaTotal      = arma::zeros<IVEC >(state.basis.ORqn.nb);
    IVEC indexTotal      = arma::zeros<IVEC >(state.basis.ORqn.nb);

    double minLambda = 0;
    double maxLambda = 0;

    arma::vec initialLambda = state.chemPot;

    bool minErrorComputed = false;
    bool maxErrorComputed = false;
    double minError = 0;
    double maxError = 0;

    while (true)
    {
      // INFO("trace rho (before)    : %9.6f", arma::trace(rhoNext(iso)) * 2.0);

      for (INT omega = 0; omega < basis.mMax; omega++)
      {
        //== Compute constraintHam
        arma::mat M_omega  = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
        UINT dimHO = M_omega.n_rows;
        UINT dimOR = M_omega.n_cols;

        arma::mat hTilde_omega     = hTilde(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));
        arma::mat constraintHam = arma::zeros(dimHO, dimHO);

        for (auto &c : state.constraints)
        {
          if ((c.second.iso == iso) || (c.second.iso == TOTAL))
          {
            if (c.second.gender == Constraint::MM)
            {
              constraintHam += multipoleOperators.ql0(c.second.lm).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega))
                               * c.second.lambda * c.second.factor;
            }
            else if (c.second.gender == Constraint::QN)
            {
              if (geometricalOperators.izNeck != -1)
              {
                constraintHam += geometricalOperators.qneck0(0).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega))
                                 * c.second.lambda * c.second.factor;
              }
              else
              {
                continue;
              }
            }
          }
        }

        hTilde_omega -= constraintHam;
        arma::mat h_omega     = M_omega.t() * hTilde_omega * M_omega;

        //== Compute particle states
        arma::vec eige;
        arma::mat eigv;
        arma::eig_sym(eige, eigv, h_omega);

        partEne(iso).elem(state.basis.omegaIndexOR(omega)) = eige;

        arma::mat e_omega = h_omega - arma::eye(dimOR, dimOR) * state.chemPot(iso);

        arma::mat delta_omega = M_omega.t() * deltaTilde(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) * M_omega;

        arma::mat hHFB = arma::zeros(dimOR * 2, dimOR * 2);
        hHFB.submat(0, 0, dimOR - 1, dimOR - 1) = e_omega;
        hHFB.submat(dimOR, dimOR, 2 * dimOR - 1, 2 * dimOR - 1) = -e_omega;
        hHFB.submat(0, dimOR, dimOR - 1, 2 * dimOR - 1) = -delta_omega;
        hHFB.submat(dimOR, 0, 2 * dimOR - 1, dimOR - 1) = -delta_omega;


        //== HFB Hamiltonian diagonalization.
        if (!Tools::checkSymmetry(hHFB, "hHFB in " + std::string(__PRETTY_FUNCTION__)))
          hHFB = arma::symmatu(hHFB);

        // Tools::info("hHFB", hHFB);

        // check for invalid values
        if (hHFB.has_nan() || hHFB.has_inf())
        {
          Tools::warning("In SolverHFBBroyden::singleHFBiter(): the hamiltonian matrix contains NaN or inf => exit");
          DBG_RETURN(false);
        }

        // Diagonalization
        arma::vec qpEne_omega;
        arma::mat B_omega;
        arma::eig_sym(qpEne_omega, B_omega, hHFB);

        // Sort by decreasing eigenvalues.
        UVEC Eind = arma::sort_index(qpEne_omega, "descend");
        B_omega = B_omega.cols(Eind);
        qpEne_omega = qpEne_omega.elem(Eind);

        // Remove E<0 quasi-particles.
        qpEne_omega = qpEne_omega.head(dimOR);

        // Extraction of U and V matrices from B.
        arma::mat U_omega = B_omega.submat(0, 0, dimOR - 1, dimOR - 1);
        arma::mat V_omega = B_omega.submat(dimOR, 0, 2 * dimOR - 1, dimOR - 1); // equivalent to -V

        // J.-P. Ebran's derivations, p. 10, Eq. (II-C.9).
        arma::mat Uho = M_omega * U_omega;
        arma::mat Vho = M_omega * V_omega;

        // reconstruction of the total rho and kappa matrices.
        rhoNext(iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega)) =  Vho * Vho.t();
        kappaNext(iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega)) =  Vho * Uho.t();

        qpEne(iso).elem(state.basis.omegaIndexOR(omega)) = qpEne_omega;
        U(iso).submat(state.basis.omegaIndexOR(omega), state.basis.omegaIndexOR(omega)) = U_omega;
        V(iso).submat(state.basis.omegaIndexOR(omega), state.basis.omegaIndexOR(omega)) = V_omega;


        state.U(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)) = Uho;
        state.V(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega)) = Vho;

        INT localId = -1;

        // blocking for rho and kappa matrices
        if (state.blockedQP(iso) > -1)
        {
          if (arma::any(arma::find(state.basis.omegaIndexOR(omega) == state.blockedQP(iso))))
          {
            localId = 0;
            while(state.basis.omegaIndexOR(omega)(localId) != state.blockedQP(iso)) localId++;
            // INFO("found state ! omega: %d %d %d", omega, localId, state.basis.omegaIndexOR(omega)(localId));
          }

          if (localId != -1)
          {
            // Do the blocking
            rhoNext(iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega))
              += 0.5 * (Uho.col(localId) * Uho.col(localId).t() - Vho.col(localId) * Vho.col(localId).t());
            kappaNext(iso).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega))
              += 0.5 * (-1.0 * Uho.col(localId) * Vho.col(localId).t() - Vho.col(localId) * Uho.col(localId).t());
          }
        } // blocking ?
      } // omega


      // Estimator for v2 ?
      VV(iso) = arma::clamp(arma::sum(arma::square(V(iso))), 0.0, 1.0).t();

      // if (state.blockedQP(iso) > -1)
      // {
      //   VV(iso)(state.blockedQP(iso)) = 0.5;
      // }

      // Store particle states.
      state.partStates(iso) = States(PF("Part. %s", iso == NEUTRON ? "n" : "p"));
      for (INT i = 0; i < partEne(iso).n_elem; i++)
      {
        INT m = basis.ORqn(0, i);
        INT s = basis.ORqn(4, i);
        std::string label = PF("%d/2", 2 * m - 2 * s + 1);

        state.partStates(iso).add(i, partEne(iso)(i), 0.0, label);
      } // i

      // Store quasi-particle states.
      state.qpStates(iso) = States(PF("QP. %s", iso == NEUTRON ? "n" : "p"));
      for (INT i = 0; i < qpEne(iso).n_elem; i++)
      {
        INT m = basis.ORqn(0, i);
        INT s = basis.ORqn(4, i);

        std::string strBlocked = "";
        if (i == state.blockedQP(iso)) strBlocked = "* ";

        std::string label = strBlocked + PF("%d/2", 2 * m - 2 * s + 1);

        state.qpStates(iso).add(i, qpEne(iso)(i), VV(iso)(i), label);
      } // i

      // Sort the quasi-particle states by ascending energy.
      state.qpStates(iso).sort();

      // optionally print qpStates
      // Tools::mesg("SolBro", state.qpStates(iso).info());

      multipoleOperators.calcQlm(rhoNext);

      double v2error = multipoleOperators.nPart(iso) - nPart(iso);

      if ((lambdaIter(iso) >= lambdaIterMax) || (fabs(v2error) < lambdaMax)) {

        // The step for the next iteration should be proportional to the change in λ
        lambdaStep(iso) = fabs(initialLambda(iso) - state.chemPot(iso)) * 10 + 0.001;
        // INFO("[λstep] updated %f", lambdaStep(iso));

        break;
      }
      lambdaIter(iso)++;

      //== Root finding
      // We first need two λ whose errors are of opposite signs: v2err(λmin) * v2err(λmax) < 0
      // Then we can apply a secant method

      // Update the interval
      if (v2error > 0.0)
      {
        maxLambda = state.chemPot(iso);
        maxError = v2error;
        maxErrorComputed = true;
      }
      else
      {
        minLambda = state.chemPot(iso);
        minError = v2error;
        minErrorComputed = true;
      }

      if (!minErrorComputed)
      {
        // Try to find a negative error
        state.chemPot(iso) -= lambdaStep(iso);
        // INFO("[λ<0] %f < %f < %f", minLambda, state.chemPot(iso), maxLambda);
      }
      else if (!maxErrorComputed)
      {
        // Try to find a positive error
        state.chemPot(iso) += lambdaStep(iso);
        // INFO("[λ>0] %f < %f < %f", minLambda, state.chemPot(iso), maxLambda);
      }
      else
      {
        // New value is the intersection of the secant and 0
        double slope = (maxLambda - minLambda) / (maxError - minError);
        state.chemPot(iso) = maxLambda - maxError * slope;
        // INFO("[<λ<] %f < %f < %f", minLambda, state.chemPot(iso), maxLambda);
      }

    } // lambda-iterations

    // optionally print qpStates
    // Tools::mesg("SolBro", state.qpStates(iso).info());

  } // iso


  lastState = state;

  //==============================================================================
  //==============================================================================

  if (!adjustConstraints())
  {
    DBG_RETURN(false);
  }

  //==============================================================================

  value = arma::norm(state.rho(NEUTRON) - rhoNext(NEUTRON), "inf")
        + arma::norm(state.rho(PROTON ) - rhoNext(PROTON ), "inf");

  ene = arma::accu(interaction.totalEnergy);

  state.eneQP = qpEne;
  state.vecOc = VV;

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constraints re-adjustment.
 */

bool SolverHFBBroyden::adjustConstraints(void)
{
  DBG_ENTER;

  //Useful quantities
  Basis &basis = state.basis;
  basis.calcWDN();
  arma::mat M = state.basis.HOtoOR.t();

  Multi<arma::mat> Q1, Q2;
  Multi<arma::mat> qpEne_inv;
  Multi<arma::mat> qpEne_sub;

  // Tools::info("rhoNextNeut (adjust)", rhoNext(NEUTRON));

  for (INT c = 0; c < state.constraints.size(); c++)
  {
    Q1(NEUTRON, c) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
    Q1(PROTON , c) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
    Q2(NEUTRON, c) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
    Q2(PROTON , c) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  }

  //Re-calculate multipole moments.
  multipoleOperators.calcQlm(rhoNext);

  if (fragInLoop)
  {
    GeometricalOperators adjGeometricalOperators(state);
    geometricalOperators = adjGeometricalOperators;
    geometricalOperators.calcQneck(state.rho(NEUTRON), state.rho(PROTON));
  }

  for (auto &c : state.constraints)
  {
    if (c.second.gender == Constraint::MM)
    {
      c.second.measuredVal = multipoleOperators.qlm(c.second.lm, c.second.iso) * c.second.factor;
    }
    else if (c.second.gender == Constraint::QN)
    {
      if (geometricalOperators.izNeck != -1)
      {
        c.second.measuredVal = geometricalOperators.qneck(0, c.second.iso) * c.second.factor;
      }
      else
      {
        c.second.measuredVal = c.second.val;
      }
    }

  }

  //==============================================================================
  // constraints: constraint matrix calculation

  Multi<arma::mat> DF;
  DF(NEUTRON) = arma::zeros(state.constraints.size(), state.constraints.size());
  DF(PROTON ) = arma::zeros(state.constraints.size(), state.constraints.size());

  for (INT iso: {NEUTRON, PROTON})
  {
    arma::mat eqp_mat = Tools::matFromCol(qpEne(iso), Q2(iso, 0).n_cols);

    qpEne_inv(iso) = 1 / (eqp_mat + eqp_mat.t());
    qpEne_sub(iso) = 1 / (eqp_mat - eqp_mat.t());
    qpEne_sub(iso).diag() = arma::zeros(qpEne_sub(iso).n_cols);

    INT ic = 0;

    for (auto &c : state.constraints)
    {
      if ((c.second.iso == iso) || (c.second.iso == TOTAL))
      {
        if (c.second.gender == Constraint::MM)
        {
          for (INT omega = 0 ; omega < basis.mMax ; omega++)
          {
            arma::mat U_omega = U(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
            arma::mat V_omega = V(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
            arma::mat M_omega = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
            arma::mat ql0_omega  = multipoleOperators.ql0(c.second.lm).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));

            Q1(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * ql0_omega * M_omega * U_omega
                - V_omega.t() * M_omega.t() * ql0_omega * M_omega * V_omega;
            Q2(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * ql0_omega * M_omega * V_omega
                + V_omega.t() * M_omega.t() * ql0_omega * M_omega * U_omega;
          }

          Q1(iso, ic) *= c.second.factor;
          Q2(iso, ic) *= c.second.factor;
        }
        else if (c.second.gender == Constraint::QN)
        {
          if (geometricalOperators.izNeck != -1)
          {
            for (INT omega = 0 ; omega < basis.mMax ; omega++)
            {
              arma::mat U_omega = U(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
              arma::mat V_omega = V(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
              arma::mat M_omega = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
              arma::mat qneck0_omega  = geometricalOperators.qneck0(0).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega));

              Q1(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * qneck0_omega * M_omega * U_omega
                  - V_omega.t() * M_omega.t() * qneck0_omega * M_omega * V_omega;
              Q2(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * qneck0_omega * M_omega * V_omega
                  + V_omega.t() * M_omega.t() * qneck0_omega * M_omega * U_omega;
            }

            Q1(iso, ic) *= c.second.factor;
            Q2(iso, ic) *= c.second.factor;
          }
          else
          {
            for (INT omega = 0 ; omega < basis.mMax ; omega++)
            {
              arma::mat U_omega = U(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
              arma::mat V_omega = V(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
              arma::mat M_omega = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));

              Q1(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * M_omega * U_omega
                  - V_omega.t() * M_omega.t() * M_omega * V_omega;
              Q2(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) = U_omega.t() * M_omega.t() * M_omega * V_omega
                  + V_omega.t() * M_omega.t() * M_omega * U_omega;
            }
          }
        }
      }

      ic++;
    }

    DF(iso) = arma::zeros(ic, ic);

    for (INT c = 0; c < ic; c++) 
    {
      for (INT d = 0; d < ic; d++) 
      {
        DF(iso)(c, d) = 2.0 * arma::accu(Q2(iso, c) % Q2(iso, d) % qpEne_inv(iso));
      }
    }
  } // iso

  arma::mat DF_tot = DF(NEUTRON) + DF(PROTON);

  // Tools::info("DF", DF_tot);

  //==============================================================================
  // constraints: lambda correction calculation
  arma::mat DF_inv;

  try
  {
    DF_inv = DF_tot.i();
  }
  catch (const std::runtime_error &e)
  {
    Tools::warning("singular matrix in SolverHFB::adjustConstraints()");
    DBG_RETURN(false);
  }

  arma::vec lambdaCor = arma::zeros(state.constraints.size());


  INT ic = 0;

  for (auto &c : state.constraints)
  {
    Constraint &cstc = c.second;
    lambdaCor(ic) = 0.0;

    INT id = 0;

    for (auto &d : state.constraints)
    {
      Constraint &cstd = d.second;
      // Eq. (97) p. 163, Berger's PhD Thesis
      lambdaCor(ic) += DF_inv(ic, id) * (cstd.val - cstd.measuredVal);
      id++;
    }

    cstc.lambda += lambdaCor(ic) * 1.0;

    ic++;
  }


  //==============================================================================

  // constraints: rho and kappa correction calculation
  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0 ; omega < basis.mMax ; omega++)
    {
      arma::mat M_omega = M.submat(basis.omegaIndexHO(omega), basis.omegaIndexOR(omega));
      arma::mat U_omega = U(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
      arma::mat V_omega = V(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));
      UINT dimOR = V_omega.n_rows;
      arma::mat SumQ1_omega = arma::zeros(dimOR, dimOR);
      arma::mat SumQ2_omega = arma::zeros(dimOR, dimOR);
      INT ic = 0;

      for (auto &c : state.constraints)
      {
        if ((c.second.iso == iso) || (c.second.iso == TOTAL))
        {
          SumQ1_omega += Q1(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) * lambdaCor(ic);
          SumQ2_omega += Q2(iso, ic).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) * lambdaCor(ic);
        }

        ic++;
      }

      // Eq. (I-6.22) in hfb2ct_Ia.pdf
      arma::mat S1 = SumQ1_omega % qpEne_sub(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega)) * -1.0;

      // Eq. (I-6.14) in hfb2ct_Ia.pdf
      arma::mat S2 = SumQ2_omega % qpEne_inv(iso).submat(basis.omegaIndexOR(omega), basis.omegaIndexOR(omega));


      /* Tools::info("S2", S2); */
      /* INFO("S2: min=%f max=%f", S2.min(), S2.max()); */


      // Eq. (I-6.19) in hfb2ct_Ia.pdf
      rhoNext(  iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += M_omega * ( V_omega * S2.t() * U_omega.t() + U_omega * S2     * V_omega.t() ) * M_omega.t();
      kappaNext(iso).submat(basis.omegaIndexHO(omega), basis.omegaIndexHO(omega)) += M_omega * ( U_omega * S2     * U_omega.t() - V_omega * S2.t() * V_omega.t() ) * M_omega.t();
    } // omega
  } // iso

  // Tools::info("rhoNextNeut", rhoNext(NEUTRON));

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string SolverHFBBroyden::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"state.", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"interaction", interaction.info(true)},
      {"momen.", multipoleOperators.info(true)},
      {"mixing", mixing.info(true)},
      {"status", Solver::statusStr[status]},
      {"block.", Tools::infoStr((state.blockedQP(NEUTRON) > -1) ||
                                (state.blockedQP(PROTON ) > -1))},
    }, true);
  }
  else
  {
    std::string nBlocked = "none";
    std::string pBlocked = "none";

    if (!state.qpStates.empty())
    {
      nBlocked =state.qpStates(NEUTRON).info(state.blockedQP(NEUTRON), 0, false);
      pBlocked =state.qpStates(PROTON ).info(state.blockedQP(PROTON ), 0, false);
    }


    result += Tools::treeStr(
    {
      {"SolverHFBBroyden", ""},
      {"state.", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"interaction", interaction.info(true)},
      {"momen.", multipoleOperators.info(true)},
      {"mixing", mixing.info(true)},
      {"status", Solver::statusStr[status]},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"target", Tools::infoStr(cvgTarget)},
      {"maxItL", Tools::infoStr(lambdaIterMax)},
      {"ltarg.", Tools::infoStr(lambdaMax)},
      {"bloc.n", nBlocked},
      {"bloc.p", pBlocked},
    }, false);
  }

  DBG_RETURN(result);
}
