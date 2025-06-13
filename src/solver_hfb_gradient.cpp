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

#include "solver_hfb_gradient.h"
#include "solver_hfb_broyden.h"
#include "plot.h"
#include "tools.h"

/** \file
 *  \brief Methods of the SolverHFBGradient class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > SolverHFBGradient::validKeys =
  {
    { "solver/gradient/cvgTarget"               , "Convergence target value"                                                  , "1e-6" , "D" },
    { "solver/gradient/maxIter"                 , "Maximum number of iterations"                                              , "20"   , "I" },
    { "solver/gradient/cvgTargetLambda"         , "Convergence target value for lambda-iterations"                            , "1e-5" , "D" },
    { "solver/gradient/maxIterLambda"           , "Maximum number of iterations for lambda-iterations"                        , "20"   , "I" },
    { "solver/gradient/randomSeed"              , "Seed used for the generation of initial random U and V matrices"           , "1337" , "I" },
    { "solver/gradient/cvgTargetSwitchToBroyden", "Convergence value under which we switch to the Broyden Mixing HFB solver"  , "1e-3" , "D" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

SolverHFBGradient::SolverHFBGradient(const std::string &filename) : SolverHFBGradient(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

SolverHFBGradient::SolverHFBGradient(const DataTree &dataTree) : SolverHFBGradient(dataTree, State(dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree.
 */

SolverHFBGradient::SolverHFBGradient(const DataTree &dataTree, State _state) : Solver(dataTree, _state),
  interaction(dataTree, &state),
  discrete(&state.basis, Mesh::regular(-10.0, 0.0, -15.0, 10.0, 0.0, 15.0,  101, 1, 151)),
  multipoleOperators(state)
{
  DBG_ENTER;

  // interaction.setDisabledInteractionFromDataTree(dataTree);
  // interaction.setCustomNodesFromDataTree(dataTree);

  dataTree.get(maxIter,                  "solver/gradient/maxIter",                  true);
  dataTree.get(cvgTarget,                "solver/gradient/cvgTarget",                true);
  dataTree.get(lambdaMax,                "solver/gradient/lambdaMax",                true);
  dataTree.get(lambdaIterMax,            "solver/gradient/lambdaIterMax",            true);
  dataTree.get(randomSeed,               "solver/gradient/randomSeed",               true);
  dataTree.get(cvgTargetSwitchToBroyden, "solver/gradient/cvgTargetSwitchToBroyden", true);

  Tools::setSeed(randomSeed);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the HFB loop.
 */

void SolverHFBGradient::init()
{
  DBG_ENTER;

  startTime = Tools::clock();

  converged = false;

  localNbIter = 0;

  if (maxIter < 1) maxIter = 1;

  // Chemical potentials (lagrange multipliers for the numbers of particles).
  mu = state.chemPot;

  // Calculate the basis orthogonalization.
  state.basis.calcWDN();

  // Initialize objects.
  if (!state.empty())
  {
    if (nbIter == 0)
      Tools::mesg("SolGra", "starting from previous rho and kappa");

    state.calcUVFromRhoKappa(dataTree);

    rho(NEUTRON) = state.rho(NEUTRON);
    rho(PROTON ) = state.rho(PROTON );
    kappa(NEUTRON) = state.kappa(NEUTRON);
    kappa(PROTON ) = state.kappa(PROTON );

    arma::mat Mi = state.basis.ORtoHO.t();

    // in HO representation
    U(NEUTRON) = Mi * state.U(NEUTRON);
    U(PROTON ) = Mi * state.U(PROTON );
    V(NEUTRON) = Mi * state.V(NEUTRON);
    V(PROTON ) = Mi * state.V(PROTON );

  }
  else
  {
    Tools::mesg("SolGra", "starting from random wavefunctions (" + PF_GREEN("seed = %d", randomSeed) + ")");
    Tools::warning("Starting from random wavefunctions increases the (small) probability of finding a local minimum");

    initRandomWF();
    calcRhoKappa(U, V, rho, kappa);

    state.rho(NEUTRON)   = rho(NEUTRON);
    state.kappa(NEUTRON) = kappa(NEUTRON);
    state.rho(PROTON)    = rho(PROTON);
    state.kappa(PROTON)  = kappa(PROTON);
  }

  // Calculate the multipole moment operator matrices.
  multipoleOperators.calcQl0Matrices();

  // If a constraint on geometrical operators is active, calculates fragment properties at each HFB iteration.
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

    if (geometricalOperators.izNeck != -1)
    {
      geometricalOperators.calcQNeckMatrix();
      geometricalOperators.calcQneck(state.rho(NEUTRON), state.rho(PROTON));
      // fragments.calcSepDistMatrix();
      // fragments.calcMassAsymMatrix();
    }
  }

  // Avoid empty kappa matrices.
  if (arma::norm(kappa(NEUTRON), "inf") < 1e-14)
  {
    Tools::mesg("SolGra", "Detected an empty kappa matrix [NEUT]");

    if (forceNonEmptyKappa)
    {
      Tools::mesg("SolGra", "Forcing non-empty kappa matrix [NEUT]");
      state.kappa(NEUTRON) = state.rho(NEUTRON);
    }
  }

  if (arma::norm(kappa(PROTON), "inf") < 1e-14)
  {
    Tools::mesg("SolGra", "Detected an empty kappa matrix [PROT]");

    if (forceNonEmptyKappa)
    {
      Tools::mesg("SolGra", "Forcing non-empty kappa matrix [PROT]");
      state.kappa(PROTON ) = state.rho(PROTON );
    }
  }

  // Tools::mesg("SolGra", "starting HFB iterations (Gradient method)");
  // Tools::mesg("SolGra", info(false));

  if (state.eneQP(NEUTRON).n_elem == 0)
  {
    state.eneQP(NEUTRON).ones(state.basis.HOqn.nb);
    state.eneQP(PROTON ).ones(state.basis.HOqn.nb);
  }

  G(0) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  G(1) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);

  mu = arma::vec{.0, 0.};
  nu = arma::vec{.001, .001};

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Prepare and perform the next HFB iteration.
 */

bool SolverHFBGradient::nextIter()
{
  DBG_ENTER;

  //============================================================================
  //===== Single HFB iteration --> recompute 3b field and adjust constraints
  //============================================================================

  // double startTime = Tools::clock();

  state.rho(NEUTRON)   = rho(NEUTRON);
  state.rho(PROTON )   = rho(PROTON);
  state.kappa(NEUTRON) = kappa(NEUTRON);
  state.kappa(PROTON ) = kappa(PROTON);

  if (!gradIter(G, mu, nu)) DBG_RETURN(false);

  state.rho(NEUTRON)   = rho(NEUTRON);
  state.rho(PROTON )   = rho(PROTON);
  state.kappa(NEUTRON) = kappa(NEUTRON);
  state.kappa(PROTON ) = kappa(PROTON);

  value = (arma::norm(G(0), "inf") + arma::norm(G(1), "inf")) / arma::sum(nu);

  // adjusting mu and nu
  for (INT iso: {NEUTRON, PROTON})
  {
    double    emin, emax;
    arma::vec indivEne = state.eneQP(iso);
    emin          = 2. * arma::min(arma::abs(indivEne));
    emax          = 4. * arma::max(arma::abs(indivEne));
    nu(iso)       = std::pow(2. / (std::sqrt(emin) + std::sqrt(emax)), 2);
    mu(iso) =
      std::pow((std::sqrt(emin) - std::sqrt(emax)) / (std::sqrt(emin) + std::sqrt(emax)), 2);
  }

  calcRhoKappa(U, V, rho, kappa);

  state.rho(NEUTRON)   = rho(NEUTRON);
  state.rho(PROTON )   = rho(PROTON);
  state.kappa(NEUTRON) = kappa(NEUTRON);
  state.kappa(PROTON ) = kappa(PROTON);

  multipoleOperators.calcQlm(rho);


  if (fragInLoop)
  {
    if (geometricalOperators.izNeck != -1)
    {
      // fragments.calcSepDist(rho(NEUTRON), rho(PROTON));
      // fragments.calcMassAsym(rho(NEUTRON), rho(PROTON));
      geometricalOperators.calcQneck(rho(NEUTRON), rho(PROTON));
    }
  }

  ene = interaction.totalEnergy(NEUTRON) + interaction.totalEnergy(PROTON);

  // double iterLength = Tools::clock() - startTime;

  //============================================================================

  // print iteration message

  std::string eneStr = ((ene > -9999.999 ) && (ene < 0.0)) ? PF_GREEN("%9.3f", ene) : PF_RED("%9.2e", ene);
  Tools::mesg("SolGra",
              PF("#it: %03d ", nbIter) +
              PF("cvg: %8.2e ", value) +
              PF("ene: ") + eneStr + " " +
              // PF("tot: %6.3fs (fld: %6.3fs)", iterLength, interaction.calcLength) + " " +
              PF("ln: %7.3f ", state.chemPot(NEUTRON)) +
              PF("lp: %7.3f ", state.chemPot(PROTON )) +
              PF("(#%2d)", lambdaIter)
              // PF_YELLOW(" " + interaction.getWarningStr())
             );

  // Tools::mesg("SolGra", multipoleOperators.getNiceInfo());
  // Tools::mesg("SolGra", interaction.getNiceInfo());

  //============================================================================

  if (useBokeh)
  {
    Plot::slot(0);
    Plot::curve("#iter", "Ene HFB [MeV]", nbIter, ene);

    Plot::slot(1);
    Plot::curve("#iter", "Convergence", nbIter, log10(value));

    if (plotDensities)
    {
      Plot::slot(2);
      arma::mat densn = discrete.getLocalXZ(rho(NEUTRON), true);
      Plot::map("HFB local density (neut)", densn, discrete.mesh);

      Plot::slot(3);
      arma::mat densp = discrete.getLocalXZ(rho(PROTON ), true);
      Plot::map("HFB local density (prot)", densp, discrete.mesh);
    }
  }

  //============================================================================

  // Continue or break ?
  bool mustContinue = true;

  if (localNbIter + 1 >= maxIter)
  {
    Tools::mesg("SolGra", "Maximum number of iterations reached, exiting HFB loop");
    mustContinue = false;
  }

  if ((value < cvgTarget) && (localNbIter > 1))
  {
    Tools::mesg("SolGra", "Target value reached, exiting HFB loop");
    mustContinue = false;
  }

  if (value > 1e20)
  {
    Tools::mesg("SolGra", "Convergence is too high, exiting HFB loop");
    mustContinue = false;
  }

  state.calculationLength = Tools::clock() - startTime;

  if (!mustContinue)
  {
    // Last HFB iteration
    finalize(G, 1.0);

    state.vecOc(NEUTRON) = arma::diagvec(state.V(NEUTRON) * state.V(NEUTRON).t());
    state.vecOc(PROTON ) = arma::diagvec(state.V(PROTON ) * state.V(PROTON ).t());

    bestEne = 0.0;

    if (localNbIter >= 0)
    {
      if ((value <= cvgTarget) && (fabs(ene) < 1e6))
      {
        bestEne = ene;
        converged  = true;
      }
    }
    state.nbIter = nbIter;

    DBG_RETURN(false);
  }

  //============================================================================

  nbIter++;
  localNbIter++;

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Perform one gradient descent iteration
 *
 * \param oG
 * \param m
 * \param n
 * \return Multi<arma::mat>
 */

bool SolverHFBGradient::gradIter(Multi<arma::mat> oG, arma::vec m, arma::vec n)
{
  DBG_ENTER;

  Multi<arma::mat> F;
  arma::vec           b, lda;
  arma::mat           A;
  Multi<arma::mat> nG = getGradient();
  Multi<arma::mat> Uj = U;
  Multi<arma::mat> Vj = V;

  /// Computing gradient
  if (lambdaIterMax >= 0)
  {
    calcConstraintsQP(U, V, F, b);

    /* INFO("F: " + F.info()); */
    /* for (auto &f: F) */
    /* { */
    /*   Tools::info("F:" + Tools::vecToStr(f.first),  f.second); */
    /* } */
    /**/
    /* INFO("nG: " + nG.info()); */
    /* for (auto &f: nG) */
    /* { */
    /*   Tools::info(PF("nG: %d", f.first[0]),  f.second); */
    /* } */

    A.resize(b.size(), b.size());
    lda = arma::zeros(b.size());

    lda(0) = state.chemPot(0);
    lda(1) = state.chemPot(1);

    INT i = 0;

    for (auto &c : state.constraints)
    {
      lda(i + 2) = c.second.lambda;
      i++;
    }

    for (INT i = 0; i < b.size(); i++)
    {
      b(i) = 2. * arma::accu(F(i, 0) % nG(0)) + 2. * arma::accu(F(i, 1) % nG(1));
    }


    for (INT i = 0; i < INT(b.size()); i++)
    {
      for (INT j = 0; j < INT(b.size()); j++)
      {
        A(i, j) = 2. * arma::accu(F(i, 0) % F(j, 0)) + 2. * arma::accu(F(i, 1) % F(j, 1));
      }
    }

    UVEC idx = arma::find(arma::diagvec(A) > 1e-7);
    lda(idx)       = arma::solve(A(idx, idx), b(idx));

#ifdef CLAMP_CHEMPOTS
    // clamp the chemical potential values (avoid some buggy points)
    lda(0) = MIN(MAX(lda(0), -15.0), 5.0);
    lda(1) = MIN(MAX(lda(1), -15.0), 5.0);
#endif

    for (INT i = 0; i < b.size(); i++)
    {
      for (INT iso: {NEUTRON, PROTON})
      {
        nG(iso) = nG(iso) - lda[i] * F(i, iso);
      }
    }
  }

  state.chemPot(0) = lda(0);
  state.chemPot(1) = lda(1);

  INT i = 0;

  for (auto &c : state.constraints)
  {
    c.second.lambda = lda(i + 2);
    i++;
  }



  G(0) = -n(0) * nG(0) + m(0) * oG(0);
  G(1) = -n(1) * nG(1) + m(1) * oG(1);

  /// obtaining |Phi(i,0)>
  if (!calcUV(U, V, Uj, Vj))
  {
    DBG_RETURN(false);
  }

  /// readjusting the constraints
  for (lambdaIter = 1; lambdaIter < lambdaIterMax; lambdaIter++)
  {
    calcConstraintsQP(Uj, Vj, F, b);

    if (lambdaIter > 0 && arma::norm(b, "inf") < lambdaMax) break;

    for (INT i = 0; i < b.size(); i++)
    {
      for (INT j = 0; j < b.size(); j++)
      {
        A(i, j) = 2. * arma::accu(F(i, 0) % F(j, 0)) + 2. * arma::accu(F(i, 1) % F(j, 1));
      }
    }

    lda.zeros(b.size());
    UVEC idx = arma::find(arma::diagvec(A) > 1e-7);

    if (arma::norm(b(idx), "inf") < lambdaMax) break;

    lda(idx) = arma::solve(A(idx, idx), b(idx));

#ifdef CLAMP_CHEMPOTS
    // clamp the chemical potential values (avoid some buggy points)
    lda(0) = MIN(MAX(lda(0), -15.0), 5.0);
    lda(1) = MIN(MAX(lda(1), -15.0), 5.0);
#endif

    for (INT i = 0; i < b.size(); i++)
    {
      for (INT iso: {NEUTRON, PROTON})
      {
        G(iso) = G(iso) - lda[i] * F(i, iso);
      }
    }

    if (!calcUV(U, V, Uj, Vj))
    {
      DBG_RETURN(false);
    }

    state.chemPot(0) -= lda(0) / n(0);
    state.chemPot(1) -= lda(1) / n(1);

#ifdef CLAMP_CHEMPOTS
    // clamp the chemical potential values (avoid some buggy points)
    state.chemPot(0) = MIN(MAX(state.chemPot(0), -15.0), 5.0);
    state.chemPot(1) = MIN(MAX(state.chemPot(1), -15.0), 5.0);
#endif

    INT i = 0;

    for (auto &c : state.constraints)
    {
      c.second.lambda -= .5 * lda(i + 2) / (n(0) + n(1));
      i++;
    }
  }

  // in OR representation
  U           = Uj;
  V           = Vj;

  arma::mat M = state.basis.HOtoOR.t();

  // in HO representation
  state.U(NEUTRON) = M * U(NEUTRON);
  state.U(PROTON ) = M * U(PROTON );
  state.V(NEUTRON) = M * V(NEUTRON);
  state.V(PROTON ) = M * V(PROTON );

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate next U and V matrices
 */

bool SolverHFBGradient::calcUV(Multi<arma::mat> &_U, Multi<arma::mat> &_V, Multi<arma::mat> &Uj, Multi<arma::mat> &Vj)
{
  DBG_ENTER;

  Multi<arma::mat> L;
  L(0) = arma::mat();
  L(1) = arma::mat();

  if (!getCholeskyDecomposition(G, L))
  {
    DBG_RETURN(false);
  }

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < state.basis.mMax; omega++)
    {
      UVEC &idxOR = state.basis.omegaIndexOR(omega);
      arma::mat   iL    = arma::inv(L(iso)(idxOR, idxOR));

      Uj(iso)(idxOR, idxOR) =
        (_U(iso)(idxOR, idxOR) - _V(iso)(idxOR, idxOR) * G(iso)(idxOR, idxOR)) * iL;
      Vj(iso)(idxOR, idxOR) =
        (_V(iso)(idxOR, idxOR) + _U(iso)(idxOR, idxOR) * G(iso)(idxOR, idxOR)) * iL;
    }
  }

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the HFB gradient
 */

Multi<arma::mat> SolverHFBGradient::getGradient()
{
  DBG_ENTER;

  Multi<arma::mat> h11, h20;
  // Force some field recalculations.
  interaction.clear();

  state.rho(NEUTRON) = arma::symmatu(state.rho(NEUTRON));
  state.rho(PROTON ) = arma::symmatu(state.rho(PROTON ));
  state.kappa(NEUTRON) = arma::symmatu(state.kappa(NEUTRON));
  state.kappa(PROTON ) = arma::symmatu(state.kappa(PROTON ));

  // Calculate the energy contributions (and the interaction).
  interaction.calcEnergies();

  h11(NEUTRON) = interaction.getHamiltonianContributions(NEUTRON, Field::DIRECT ) + interaction.getHamiltonianContributions(NEUTRON, Field::EXCHANGE);
  h11(PROTON ) = interaction.getHamiltonianContributions(PROTON, Field::DIRECT ) + interaction.getHamiltonianContributions(PROTON, Field::EXCHANGE);
  h20(NEUTRON) = interaction.getHamiltonianContributions(NEUTRON, Field::PAIRING);
  h20(PROTON ) = interaction.getHamiltonianContributions(PROTON, Field::PAIRING);

  arma::mat ene(state.basis.ORqn.nb, 2, arma::fill::zeros);
  Multi<arma::mat> G;

  G(0) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  G(1) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);

  // Multi<arma::mat>  tr(basis);
  // J.-P. Ebran's derivations, p. 8, Eq. (II-B.18).
  arma::mat Mtotal = state.basis.HOtoOR.t();

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < state.basis.mMax; omega++)
    {
      UVEC &idxHO = state.basis.omegaIndexHO(omega);
      UVEC &idxOR = state.basis.omegaIndexOR(omega);

      arma::mat M = Mtotal.submat(idxHO, idxOR);

      arma::mat lU     = U(iso)(idxOR, idxOR);
      arma::mat lV     = V(iso)(idxOR, idxOR);
      arma::mat lh     = M.t() * h11(iso)(idxHO, idxHO) * M;

      arma::mat lcons  = calcConstraints(iso, idxHO, M);

      arma::mat ldelta = M.t() * h20(iso)(idxHO, idxHO) * M;

      arma::mat H20, H11;
      convFromSPtoQP(lh, ldelta, H11, H20, lU, lV);

      G(iso)
      (idxOR, idxOR) = H20;

      H11 += lU.t() * lcons * lU - lV.t() * lcons * lV;

      arma::vec eigVals;
      arma::mat eigVecs;

      Tools::eig_sym(eigVals, eigVecs, arma::symmatu(H11));
      ene(idxOR, UVEC{UINT(iso)}) = eigVals;
    }
  }

  state.eneQP(NEUTRON) = ene.col(NEUTRON);
  state.eneQP(PROTON ) = ene.col(PROTON);

  DBG_RETURN(G);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Final treatment
 */

void SolverHFBGradient::finalize(Multi<arma::mat> &G, double angle)
{
  DBG_ENTER;

  Multi<arma::mat> h11, h20;

  h11(NEUTRON) = interaction.getHamiltonianContributions(NEUTRON, Field::DIRECT ) + interaction.getHamiltonianContributions(NEUTRON, Field::EXCHANGE);
  h11(PROTON ) = interaction.getHamiltonianContributions(PROTON, Field::DIRECT ) + interaction.getHamiltonianContributions(PROTON, Field::EXCHANGE);
  h20(NEUTRON) = interaction.getHamiltonianContributions(NEUTRON, Field::PAIRING);
  h20(PROTON ) = interaction.getHamiltonianContributions(PROTON, Field::PAIRING);

  // J.-P. Ebran's derivations, p. 8, Eq. (II-B.18).
  arma::mat Mtotal = state.basis.HOtoOR.t();

  arma::mat ene(state.basis.ORqn.nb, 2, arma::fill::zeros);
  Multi<arma::mat> tr;
  tr(0) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  tr(1) = arma::zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < state.basis.mMax; omega++)
    {
      UVEC &idxHO = state.basis.omegaIndexHO(omega);
      UVEC &idxOR = state.basis.omegaIndexOR(omega);

      arma::mat M = Mtotal.submat(idxHO, idxOR);

      arma::mat lU     = U(iso)(idxOR, idxOR);
      arma::mat lV     = V(iso)(idxOR, idxOR);
      arma::mat lh     = M.t() * h11(iso)(idxHO, idxHO) * M + calcConstraints(iso, idxHO, M);
      arma::mat ldelta = M.t() * h20(iso)(idxHO, idxHO) * M;

      arma::mat H20, H11;
      convFromSPtoQP(lh, ldelta, H11, H20, lU, lV);

      arma::vec eigenvalues;
      arma::mat eigenvectors;

      Tools::eig_sym(eigenvalues, eigenvectors,
                     arma::symmatu(H11) * angle +
                     (1. - angle) * arma::diagmat(ene(idxOR, UVEC{UINT(iso)})));

      tr(iso)(idxOR, idxOR) = eigenvectors;
    }

    U(iso) = U(iso) * tr(iso);
    V(iso) = V(iso) * tr(iso);
    G(iso) = tr(iso).t() * G(iso) * tr(iso);
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Make a Cholesky decomposition of the gradient
 *
 * \param G
 * \return Multi<arma::mat>
 */

bool SolverHFBGradient::getCholeskyDecomposition(Multi<arma::mat> &_G, Multi<arma::mat> &L)
{
  DBG_ENTER;

  for (auto &l : L)
  {
    l.second.zeros(state.basis.ORqn.nb, state.basis.ORqn.nb);
  }

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < state.basis.mMax; omega++)
    {
      UVEC &idxOR = state.basis.omegaIndexOR(omega);

      arma::mat X = arma::eye(idxOR.n_elem, idxOR.n_elem) + _G(iso)(idxOR, idxOR).t() * _G(iso)(idxOR, idxOR);

      arma::mat Ltemp;

      if (!arma::chol(Ltemp, X))
      {
        Tools::mesg("SolGra", "Cholesky decomposition not found in SolverHFBGradient::getCholeskyDecomposition()");

        DBG_RETURN(false);
      }

      L(iso)(idxOR, idxOR) = Ltemp;
    }
  }

  DBG_RETURN(true);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculation of the constraints.
 *
 * \param iso neutron or proton
 * \param idxHO Indices of the subblock
 *
 * \return matrix matrix containing \f$- \sum_c \lambda_c F_c \f$
 */

arma::mat SolverHFBGradient::calcConstraints(UINT iso, UVEC &idxHO, arma::mat &M)
{
  DBG_ENTER;

  arma::mat e = -state.chemPot(iso) * arma::eye(M.n_cols, M.n_cols);

  for (auto &c : state.constraints)
  {
    if (c.second.gender == Constraint::SD && geometricalOperators.izNeck == -1) continue;

    if (c.second.gender == Constraint::MA && geometricalOperators.izNeck == -1) continue;

    if (c.second.gender == Constraint::QN && geometricalOperators.izNeck == -1) continue;

    if ((c.second.iso == iso) || (c.second.iso == TOTAL))
    {
      if (c.second.itermax < 0 || localNbIter < c.second.itermax)
      {
        if (c.second.gender == Constraint::MM) e -= c.second.lambda * M.t() * multipoleOperators.ql0(c.second.lm)(idxHO, idxHO) * M;

        // if (c.second.gender == SD && fragments.izNeck != -1) e -= c.second.lambda * M.t() * fragments.sepdist0(0)(idxHO, idxHO) * M;
        // if (c.second.gender == MA && fragments.izNeck != -1) e -= c.second.lambda * M.t() * fragments.massasym0(0)(idxHO, idxHO) * M;
        if (c.second.gender == Constraint::QN) e -= c.second.lambda * M.t() * geometricalOperators.qneck0(0)(idxHO, idxHO) * M;
      }
    }
  }

  DBG_RETURN(e);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Compute rho and kappa matrices from U and V
 *
 * \param U
 * \param V
 * \param rho
 * \param kappa
 */

void SolverHFBGradient::calcRhoKappa(Multi<arma::mat> &U, Multi<arma::mat> &V, Multi<arma::mat> &rho, Multi<arma::mat> &kappa)
{
  DBG_ENTER;

  arma::mat Mtotal = state.basis.HOtoOR.t();

  for (INT iso: {NEUTRON, PROTON}) 
  {
    UINT nb = state.basis.HOqn.nb;
    rho(iso) = arma::zeros(nb, nb);
    kappa(iso) = arma::zeros(nb, nb);

    for (INT omega = 0; omega < state.basis.mMax; omega++) 
    {
      UVEC &idxOR = state.basis.omegaIndexOR(omega);
      UVEC &idxHO = state.basis.omegaIndexHO(omega);
      arma::mat   M     = Mtotal(idxHO, idxOR);

      rho(iso)(idxHO, idxHO)   = arma::symmatu(M * V(iso)(idxOR, idxOR) * V(iso)(idxOR, idxOR).t() * M.t());
      kappa(iso)(idxHO, idxHO) = arma::symmatu(M * V(iso)(idxOR, idxOR) * U(iso)(idxOR, idxOR).t() * M.t());

      // if (!Tools::checkSymmetry(rho(iso)(idxHO, idxHO), PF("rho iso:%01d omega:%02d", iso, omega)))
      //   rho(iso)(idxHO, idxHO) = arma::symmatu(rho(iso)(idxHO, idxHO));
      //
      // if (!Tools::checkSymmetry(kappa(iso)(idxHO, idxHO), PF("kappa iso:%01d omega:%02d", iso, omega)))
      //   kappa(iso)(idxHO, idxHO) = arma::symmatu(kappa(iso)(idxHO, idxHO));
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Compute the constraint operators in quasi-particle basis
 *
 * \param U
 * \param V
 * \param F
 * \param dF
 */

void SolverHFBGradient::calcConstraintsQP(Multi<arma::mat> &U, Multi<arma::mat> &V, Multi<arma::mat> &F, arma::vec &dF)
{
  DBG_ENTER;

  UINT   size = 2 + state.constraints.size();     // 2 additional places for p and n number

  Multi<arma::mat> rho, kappa;
  Multi<arma::vec> Eqp;
  Eqp = state.eneQP;

  calcRhoKappa(U, V, rho, kappa);

  arma::mat Mtotal = state.basis.HOtoOR.t();

  UINT nb = state.basis.ORqn.nb;

  for (INT i = 0; i < size; i++)
  {
    F(i, 0) = arma::zeros(nb, nb);
    F(i, 1) = arma::zeros(nb, nb);
  }

  dF = arma::zeros(size);

  arma::mat s = state.basis.ORtoHO * state.basis.ORtoHO.t();

  dF(NEUTRON) = 2. * arma::accu(s % rho(NEUTRON)) - (double)(state.sys.nNeut);
  dF(PROTON)  = 2. * arma::accu(s % rho(PROTON )) - (double)(state.sys.nProt);

  for (INT iso: {NEUTRON, PROTON})
  {
    for (INT omega = 0; omega < state.basis.mMax; omega++)
    {
      UVEC &idxOR = state.basis.omegaIndexOR(omega);

      F(iso, iso )(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * V(iso)(idxOR, idxOR) +
                                   V(iso)(idxOR, idxOR).t() * U(iso)(idxOR, idxOR);
    }
  }

  if (fragInLoop)
  {
    GeometricalOperators adjGeometricalOperators(state);
    geometricalOperators = adjGeometricalOperators;
    geometricalOperators.calcQneck(rho(NEUTRON), rho(PROTON));
  }

  INT i = 0;

  for (auto &c : state.constraints)
  {
    if (c.second.itermax > 0 && c.second.itermax < localNbIter) continue;

    if (c.second.itermax == 0) continue;

    c.second.measuredVal = 0.;

    for (INT iso: {NEUTRON, PROTON}) 
    {
      if ((c.second.iso == iso) || (c.second.iso == TOTAL))
      {
        for (INT omega = 0; omega < state.basis.mMax; omega++)
        {
          UVEC &idxOR = state.basis.omegaIndexOR(omega);
          UVEC &idxHO = state.basis.omegaIndexHO(omega);
          arma::mat   M     = Mtotal(idxHO, idxOR);

          if (c.second.gender == Constraint::MM)
          {
            F(i + 2, iso)(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * M.t() * multipoleOperators.ql0(c.second.lm)(idxHO, idxHO) *
                                          M * V(iso)(idxOR, idxOR) +
                                          V(iso)(idxOR, idxOR).t() * M.t() * multipoleOperators.ql0(c.second.lm)(idxHO, idxHO) *
                                          M * U(iso)(idxOR, idxOR);
            c.second.measuredVal += 2. * arma::trace(multipoleOperators.ql0(c.second.lm)(idxHO, idxHO) * rho(iso)(idxHO, idxHO).t());     // Constraint operator value at this stage
            dF(i + 2) = c.second.measuredVal - c.second.val;     // Difference between constrain and actual value
          }

          // if (c.second.gender == Constraint::SD)
          // {
          //   F(i + 2)(iso)(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * M.t() * fragments.sepdist0(c.second.lm)(idxHO, idxHO) *
          //                       M * V(iso)(idxOR, idxOR) +
          //                       V(iso)(idxOR, idxOR).t() * M.t() * fragments.sepdist0(c.second.lm)(idxHO, idxHO) *
          //                       M * U(iso)(idxOR, idxOR);
          //   c.second.measuredVal +=
          //     2. *
          //     arma::trace(fragments.sepdist0(c.second.lm)(idxHO, idxHO) *
          //               rho(iso)(idxHO, idxHO).t());     // Constraint operator value at this stage
          // }
          // if (c.second.gender == Constraint::MA)
          // {
          //   F(i + 2)(iso)(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * M.t() * fragments.massasym0(c.second.lm)(idxHO, idxHO) *
          //                       M * V(iso)(idxOR, idxOR) +
          //                       V(iso)(idxOR, idxOR).t() * M.t() * fragments.massasym0(c.second.lm)(idxHO, idxHO) *
          //                       M * U(iso)(idxOR, idxOR);
          //   c.second.measuredVal +=
          //     2. *
          //     arma::trace(fragments.massasym0(c.second.lm)(idxHO, idxHO) *
          //               rho(iso)(idxHO, idxHO).t());     // Constraint operator value at this stage
          // }
          if (c.second.gender == Constraint::QN)
          {
            if (geometricalOperators.izNeck != -1)
            {
              F(i + 2, iso)(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * M.t() * geometricalOperators.qneck0(0)(idxHO, idxHO)
                                            * M * V(iso)(idxOR, idxOR)
                                            + V(iso)(idxOR, idxOR).t() * M.t() * geometricalOperators.qneck0(0)(idxHO, idxHO)
                                            * M * U(iso)(idxOR, idxOR);
              c.second.measuredVal += 2. * arma::trace(geometricalOperators.qneck0(0)(idxHO, idxHO) * rho(iso)(idxHO, idxHO).t());     // Constraint operator value at this stage
              dF(i + 2) = c.second.measuredVal - c.second.val;     // Difference between constrain and actual value
            }
            else
            {
              F(i + 2, iso)(idxOR, idxOR) = U(iso)(idxOR, idxOR).t() * M.t()
                                            * M * V(iso)(idxOR, idxOR)
                                            + V(iso)(idxOR, idxOR).t() * M.t()
                                            * M * U(iso)(idxOR, idxOR);
              dF(i + 2) = 0.0;     // Until when the position of the neck is actually found, the contribution to the constraint correlation matrix is ignored.
            }
          }
        }
      }
    }

    i++;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Conversion from sp to qp
 */

void SolverHFBGradient::convFromSPtoQP(arma::mat &sp11, arma::mat &sp20, arma::mat &qp11, arma::mat &qp20,
                                       arma::mat &U, arma::mat &V)
{
  DBG_ENTER;

  qp11 = U.t() * sp11 * U - V.t() * sp11 * V - U.t() * sp20 * V - V.t() * sp20 * U;
  qp20 = U.t() * sp11 * V + V.t() * sp11 * U + U.t() * sp20 * U - V.t() * sp20 * V;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Conversion from qp to sp
 */

void SolverHFBGradient::convFromQPtoSP(arma::mat &sp11, arma::mat &sp20, arma::mat &qp11, arma::mat &qp20,
                                       arma::mat &U, arma::mat &V)
{
  DBG_ENTER;

  sp11 = U * qp11 * U.t() - V * qp11 * V.t() + U * qp20 * V.t() + V * qp20 * U.t();
  sp20 = -U * qp11 * V.t() + V * qp11 * U.t() - U * qp20 * U.t() + V * qp20 * V.t();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Build randomized density matrices
 */

void SolverHFBGradient::initRandomWF(void)
{
  DBG_ENTER;

  // The arma::randu() function is not fully deterministic.
  // Its results may depend on the compilation flags for example.
  // Thus, we implemented a deterministic equivalent function Tools::randMat().
  // Use it instead of arma::randu().
  //
  // The associated seed initializer function is Tools::setSeed().

  std::string test;

  // arma::vec v2 = .1 * arma::randu(state.basis.ORqn.nb);
  arma::vec v2 = .1 * Tools::randMat(state.basis.ORqn.nb, 1);

  arma::vec v2p = v2, v2n = v2;
  v2p.head(state.sys.nProt / 2) += 1.;
  v2n.head(state.sys.nNeut / 2) += 1.;
  v2p = arma::clamp(.5 * (double)(state.sys.nProt) / arma::sum(v2p) * v2p, 0., 1.);
  v2n = arma::clamp(.5 * (double)(state.sys.nNeut) / arma::sum(v2n) * v2n, 0., 1.);

  arma::mat un = arma::diagmat(arma::sqrt(1 - v2n));
  arma::mat up = arma::diagmat(arma::sqrt(1 - v2p));
  arma::mat vn = arma::diagmat(arma::sqrt(v2n));
  arma::mat vp = arma::diagmat(arma::sqrt(v2p));


  UINT      nb = state.basis.ORqn.nb;
  arma::mat C(nb, nb, arma::fill::zeros);
  arma::mat D(nb, nb, arma::fill::zeros);

  for (INT omega = 0; omega < state.basis.mMax; omega++)
  {
    UVEC &idx = state.basis.omegaIndexOR(omega);
    UINT        N   = idx.n_elem;
    // arma::mat   S   = .1 * arma::randu(N, N);
    arma::mat   S   = .1 * Tools::randMat(N, N);
    // arma::mat   SS  = .1 * arma::randu(N, N);
    arma::mat   SS  = .1 * Tools::randMat(N, N);
    S               = S - S.t();
    SS              = SS - SS.t();

    S  = arma::expmat(S);
    SS = arma::expmat(SS);

    D(idx, idx) = S;
    C(idx, idx) = SS;
  }

  // in OR representation
  U(PROTON)  = C * arma::diagmat(arma::sqrt(1 - v2p)) * D;
  U(NEUTRON) = C * arma::diagmat(arma::sqrt(1 - v2n)) * D;
  V(PROTON)  = C * arma::diagmat(arma::sqrt(v2p)) * D;
  V(NEUTRON) = C * arma::diagmat(arma::sqrt(v2n)) * D;

  // Tools::info("U(NEUTRON)", U(NEUTRON));
  // Tools::info("U(PROTON )", U(PROTON ));
  // Tools::info("V(NEUTRON)", V(NEUTRON));
  // Tools::info("V(PROTON )", V(PROTON ));

  arma::mat M = state.basis.HOtoOR.t();

  // in HO*OR representation
  state.U(NEUTRON) = M * U(NEUTRON);
  state.U(PROTON ) = M * U(PROTON);
  state.V(NEUTRON) = M * V(NEUTRON);
  state.V(PROTON ) = M * V(PROTON);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string SolverHFBGradient::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"bogoSt", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"interaction", interaction.info(true)},
      {"multipoleOp.", multipoleOperators.info(true)},
      //  {"geomOp.", geometricalOperators.info(true)},
    },
    true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"SolverHFBGradient", ""},
      {"bogoSt", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"interaction", interaction.info(true)},
      {"multipoleOp.", multipoleOperators.info(true)},
      // {"geomOp.", geometricalOperators.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"target", Tools::infoStr(cvgTarget)},
      {"litMax", Tools::infoStr(lambdaIterMax)},
      {"ltarg.", Tools::infoStr(lambdaMax)},
    },
    false);
  }

  DBG_RETURN(result);
}
