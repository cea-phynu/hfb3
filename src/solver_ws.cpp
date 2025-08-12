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

#include "solver_ws.h"
#include "interaction.h"
#include "tools.h"
#include "plot.h"

/** \file
 *  \brief Methods of the SolverWS class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > SolverWS::validKeys =
  {
    { "solver/ws/beta20tInit", "Temporary constraint on beta20t.", "0.1"   , "D" },
    { "solver/ws/q30tInit"   , "Temporary constraint on q30t."   , "1000.0", "D" },
    { "solver/ws/maxIter"    , "Maximum number of iterations"    , "50"    , "I" },
    { "solver/ws/cvgTarget"  , "Error target value."             , "1e-6"  , "D" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a filename.
 */

SolverWS::SolverWS(const std::string filename) : SolverWS(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree object.
 */

SolverWS::SolverWS(const DataTree &_dataTree) : SolverWS(_dataTree, State(_dataTree))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** The constructor from a DataTree object.
 */

SolverWS::SolverWS(const DataTree &_dataTree, State _state) : Solver(_dataTree, _state),
  wsInteraction(dataTree, &_state)
{
  DBG_ENTER;

  multipoleOperators = MultipoleOperators(state);
  discrete = Discrete(&state.basis, Mesh::regular(-10.0, 0.0, -15.0, 10.0, 0.0, 15.0, 51, 1, 51));

  dataTree.get(beta20tInit, "solver/ws/beta20tInit" , true);
  dataTree.get(q30tInit   , "solver/ws/q30tInit"    , true);
  dataTree.get(maxIter    , "solver/ws/maxIter"     , true);
  dataTree.get(cvgTarget  , "solver/ws/cvgTarget"   , true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Build and solve a Woods-Saxon hamiltonian.
 */

void SolverWS::calcWS()
{
  DBG_ENTER;

  wsInteraction = Interaction("WS", &state);
  calcWSWithInteraction(wsInteraction);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Build and solve a Woods-Saxon hamiltonian, using an external Interaction instance.
 */

void SolverWS::calcWSWithInteraction(Interaction &interactionWS)
{
  DBG_ENTER;

  // shortcut
  UINT HOnb = state.basis.HOqn.nb;

  state.basis.calcWDN();

  interactionWS("woodsSaxon")->setDef(currentDef);

  interactionWS.calcFields(true);

  arma::mat MTotal = state.basis.HOtoOR.t();
  arma::mat HnTotal = interactionWS.getHamiltonianContributions(NEUTRON, Field::DIRECT);
  arma::mat HpTotal = interactionWS.getHamiltonianContributions(PROTON, Field::DIRECT);

  // store energies and states for later re-ordering
  indivEnerWS = arma::ones(HOnb, 2) * 99.0;
  indivStateWS = arma::zeros(HOnb, HOnb, 2);
  UVEC toto = {NEUTRON};

  for (INT omega = 0; omega < state.basis.mMax; omega++)
  {
    arma::mat M  =  MTotal.submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexOR(omega));
    arma::mat Hn = HnTotal.submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega));
    arma::mat Hp = HpTotal.submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexHO(omega));

    Tools::checkSymmetry(Hn, "Hn in " + std::string(__PRETTY_FUNCTION__));
    Tools::checkSymmetry(Hp, "Hp in " + std::string(__PRETTY_FUNCTION__));

    arma::mat Hntilde = M.t() * Hn * M;
    arma::mat Hptilde = M.t() * Hp * M;
    arma::vec En;
    arma::vec Ep;
    arma::mat Btilde;

    Tools::checkSymmetry(Hntilde, "Hntilde in " + std::string(__PRETTY_FUNCTION__));
    Hntilde = arma::symmatu(Hntilde);
    arma::eig_sym(En, Btilde, Hntilde); // diago neutrons

    arma::mat Bn = M * Btilde;
    toto(0) = NEUTRON;
    indivEnerWS(state.basis.omegaIndexOR(omega), toto) = En;
    indivStateWS.slice(NEUTRON).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexOR(omega)) = Bn;

    Tools::checkSymmetry(Hptilde, "Hptilde in " + std::string(__PRETTY_FUNCTION__));
    Hptilde = arma::symmatu(Hptilde);
    arma::eig_sym(Ep, Btilde, Hptilde); // diago protons

    arma::mat Bp = M * Btilde;
    toto(0) = PROTON;
    indivEnerWS(state.basis.omegaIndexOR(omega), toto) = Ep;
    indivStateWS.slice(PROTON ).submat(state.basis.omegaIndexHO(omega), state.basis.omegaIndexOR(omega)) = Bp;
  }

  UVEC sortIndex;

  //indivEnerWS.print("indiv E");

  // sort individual states by increasing energy
  for (INT iso: {NEUTRON, PROTON})
  {
    sortIndex = arma::sort_index(indivEnerWS.col(iso));
    toto(0) = iso;
    indivEnerWS.col(iso) = indivEnerWS(sortIndex, toto);
    arma::mat tutu = indivStateWS.slice(iso);
    indivStateWS.slice(iso) = tutu.cols(sortIndex);
  }

  // indivEnerWS.print("indiv E (sorted)");
  // Particle numbers.
  arma::vec nPart(2);
  nPart(NEUTRON) = state.sys.nNeut;
  nPart(PROTON ) = state.sys.nProt;

  // calculate occupations
  indivOccupWS = arma::zeros(state.basis.HOqn.nb, 2);

#define WS_PAIRING
#ifdef WS_PAIRING
  // Nilsson and Ragnarsson (1995)
  double mass = nPart(NEUTRON) + nPart(PROTON);
  arma::vec pairingStrength = {20.5 / mass, 26.0 / mass};

  for (INT iso: {NEUTRON, PROTON})
  {
    // fermi level calculation
    double muMin = -500.0;
    double muMax = 500.0;

    while (muMax - muMin > 1e-4)
    {
      double mu = (muMin + muMax) / 2.0;
      double v = arma::accu(0.5 * (1.0 - (indivEnerWS.col(iso) - mu) / arma::sqrt(arma::pow(indivEnerWS.col(iso) - mu, 2) + pow(pairingStrength(iso), 2)))) - nPart(iso);

      if (v < 0.0) muMin = mu;
      else muMax = mu;
    }

    indivOccupWS.col(iso) = 0.5 * (1.0 - (indivEnerWS.col(iso) - muMin) / arma::sqrt(arma::pow(indivEnerWS.col(iso) - muMin, 2) + pow(pairingStrength(iso), 2)));
  }

#else

  for (INT iso: {NEUTRON, PROTON})
  {
    for (UINT l = 0; l < indivEnerWS.n_rows; l++)
    {
      indivOccupWS(l, iso) = l <= nPart(iso) ? 1.0 : 0.0;
    }
  }

#endif
  state.rho(NEUTRON) = (indivStateWS.slice(NEUTRON) * arma::diagmat(indivOccupWS.col(NEUTRON)) * indivStateWS.slice(NEUTRON).t()) / 2.0;
  state.rho(PROTON) = (indivStateWS.slice(PROTON ) * arma::diagmat(indivOccupWS.col(PROTON )) * indivStateWS.slice(PROTON ).t()) / 2.0;

  // TODO : find a better way
  state.kappa(NEUTRON) = state.rho(NEUTRON) * 0.1;
  state.kappa(PROTON) = state.rho(PROTON) * 0.1;

  // re-calculate multipole moments.
  multipoleOperators = MultipoleOperators(state);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the error on the constraint values.
 */

double SolverWS::getError(double scalingFactor)
{
  DBG_ENTER;

  for(auto &c : state.constraints)
  {
    c.second.measuredVal = multipoleOperators.qlm(c.second.lm, TOTAL);
  }

  double wsError = 0.0;
  for(auto &c : state.constraints)
  {
    wsError += fabs(c.second.val - c.second.measuredVal) / pow(scalingFactor, c.second.lm);
    // INFO("%s %f %f", c.first.c_str(), c.second.val, c.second.measuredVal);
  }

  DBG_RETURN(wsError);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice std::string representation of a deformation vector.
 */

std::string SolverWS::niceStr(Multi<double> &def)
{
  DBG_ENTER;

  std::string result = "(";
  bool ifirst = true;

  const std::vector<std::string> deformationNames = {"a00", "a10", "a20", "a30", "a40", "a50", "a60"};

  if (def.empty())
  {
    result += ")";
    DBG_RETURN(result);
  }

  for( auto &c : state.constraints)
  {

#ifdef WS_IGNORE_Q10_CONSTRAINT
    if (c.second.lm == 1) continue;
#endif

    if (!ifirst) result += ", ";

    result += PF("%s%s:%s%9.2e%s", Tools::color("yellow").c_str(), deformationNames[c.second.lm].c_str(), Tools::color("blue").c_str(), def(c.second.lm, 0), Tools::color().c_str());
    ifirst = false;
  }
  result += ")";

  DBG_RETURN(result);
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Change the current deformation.
 */

void SolverWS::updateCurrentDef(const arma::vec & v)
{
  DBG_ENTER;

  INT i = 0;
  for (auto &key : currentDef.getKeys())
  {
    currentDef(key) = v(i); // the keys are sorted in getKeys()
    i++;
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Initialize the Woods-Saxon optimization loop.
 */

void SolverWS::init()
{
  DBG_ENTER;

  startTime = Tools::clock();

  bool constraintOnQ20t = false;
  bool constraintOnQ30t = false;

  for (auto &c : state.constraints)
  {
#ifdef WS_IGNORE_Q10_CONSTRAINT
    if(c.second.lm == 1) continue;
#endif

    currentDef(c.second.lm, 0) = 0.0;
    // INFO("currentDef: %d %e", c.second.lm, c.second.val);

    if (c.second.lm == 2) constraintOnQ20t = true;
    if (c.second.lm == 3) constraintOnQ30t = true;
  }

  dim = currentDef.size();

  originalConstraints = state.constraints;

  if (!constraintOnQ20t) // no constraint on Q20t given by the user => we use beta20tInit as a temporary constraint
  {
    double nPart = state.sys.nProt + state.sys.nNeut;
    double q20tValue = MultipoleOperators::getQ20FromBeta(nPart, nPart, beta20tInit);
    state.constraints["q20t"] = Constraint("q20t", q20tValue);
    Tools::mesg("SolWS.", PF("Temporary constraint created <beta20t>=%.3f (converted to <q20t>=%.3f)", beta20tInit, q20tValue));

    currentDef(2, 0) = 0.0;
    dim = currentDef.size();
  }

  if (!constraintOnQ30t) // no constraint on Q30t given by the user => we use q30tInit as a temporary constraint
  {
    state.constraints["q30t"] = Constraint("q30t", q30tInit);
    Tools::mesg("SolWS.", PF("Temporary constraint created <q30t>=%.3f", q30tInit));

    currentDef(3, 0) = 0.0;
    dim = currentDef.size();
  }

  vmin  = - arma::ones(dim) * 1.00;
  vmax  = arma::ones(dim) * 1.00;
  vstep = arma::ones(dim) * 0.05;

  gradientWalk.xMin           = vmin;
  gradientWalk.xMax           = vmax;
  gradientWalk.xStep          = vstep;
  gradientWalk.initialCoords  = arma::ones(dim) * 0.0;

  gradientWalk.init();

  Tools::mesg("SolWS.", info(false));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Plot some quantities using Bokeh.
 */

void SolverWS::bokehPlot(void)
{
  DBG_ENTER;

  if (useBokeh)
  {
    Plot::slot(0);
    Plot::curve("a20", "q20", currentDef(2, 0), multipoleOperators.qlm(2));

    Plot::slot(1);
    Plot::curve("a20", "Error", currentDef(2, 0), value);

    Plot::slot(2);
    arma::mat densn = discrete.getLocalXZ(state.rho(NEUTRON), true);
    Plot::map("WS local density (neut)", densn, discrete.mesh);

    Plot::slot(3);
    arma::mat densp = discrete.getLocalXZ(state.rho(PROTON), true);
    Plot::map("WS local density (prot)", densp, discrete.mesh);

    Plot::slot(4);
    arma::mat potn = wsInteraction("woodsSaxon")->getWSPot().directPotentialNeut;
    Plot::map("WS potentian (neut)", potn, wsInteraction("woodsSaxon")->getWSPot().mesh);

    Plot::slot(5);
    arma::mat potp = wsInteraction("woodsSaxon")->getWSPot().directPotentialProt;
    Plot::map("WS potentian (prot)", potp, wsInteraction("woodsSaxon")->getWSPot().mesh);

  }

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Calculate the next iteration.
 */

bool SolverWS::nextIter()
{
  DBG_ENTER;

  arma::vec next = gradientWalk.getEval();

  if (next.empty())
  {
    Tools::mesg("SolWS.", "Minimum found, exiting WS loop");

    status = Solver::CONVERGED;
    finalize();

    DBG_RETURN(false);
  }
  else
  {
    updateCurrentDef(next);
    calcWS();
    value = getError();

    bokehPlot();

    gradientWalk.addEval(next, -value, true);
    //Tools::mesg("SolWS.", getHistTable());
  }

  bool isBest = false;

  if (value < bestError)
  {
    bestError = value;
    bestDef = currentDef;
    bestState = state;
    isBest = true;
  }

  Tools::mesg("SolWS.",
              PF("#it: %03d ", nbIter) +
              PF("def: ") + niceStr(currentDef) + " " +
              PF("err: %9.2e ", value) +
              (isBest ? PF_GREEN(" *") : "")
              );

  // optional
  // Tools::mesg(PF("WS_%03d", nbIter), multipoleOperators.getNiceInfo());

  nbIter++;
  status = Solver::ITERATING;

  if (nbIter >= maxIter)
  {
    Tools::mesg("SolWS.", PF_RED("maximum number of iterations reached (%d)", nbIter));

    status = Solver::MAXITER;

    finalize();

    DBG_RETURN(false);
  }


  if (value < cvgTarget)
  {
    Tools::mesg("SolWS.", PF_GREEN("target value reached %e < %e", value, cvgTarget));

    status = Solver::CONVERGED;

    finalize();

    state.converged = true;

    DBG_RETURN(false);
  }


  DBG_RETURN(true);
}


//==============================================================================
//==============================================================================
//==============================================================================

void SolverWS::finalize(void)
{
  DBG_ENTER;

  state = bestState;

  // Tools::mesg("SolWS.", "End of optimization");
  // Tools::mesg(PF("WS_res", nbIter), niceStr(bestDef) + PF(" err:%9.2e", bestError) + PF_GREEN(" *"));

  // Restore the original constraints.
  state.constraints = originalConstraints;

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate a Bogoliubov state solution of the WS using SolverWS.
 */

State SolverWS::calc()
{
  DBG_ENTER;

  init();
  while (nextIter());

  calculationLength = Tools::clock() - startTime;
  state.calculationLength = calculationLength;
  state.nbIter = nbIter;

  Tools::mesg("SolWS.", multipoleOperators.getNiceInfo());

  DBG_RETURN(state);
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string SolverWS::info(bool isShort) const
{
  DBG_ENTER;
  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"state ", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"maxIt.", Tools::infoStr(maxIter)},
      {"status", Solver::statusStr[status]},
    }, true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"SolverWS", ""},
      {"state ", state.info(true)},
      {"basis ", state.basis.info(true)},
      {"max.it", Tools::infoStr(maxIter)},
      {"target", Tools::infoStr(cvgTarget)},
      {"status", Solver::statusStr[status]},
      {"dim.  ", Tools::infoStr(dim)},
      {"momen.", multipoleOperators.info(true)},
      {"vmin  ", Tools::vecToStr(vmin)},
      {"vmax  ", Tools::vecToStr(vmax)},
      {"vstep ", Tools::vecToStr(vstep)},
    }, false);
  }

  DBG_RETURN(result);
}
