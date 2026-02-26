
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

#include "action.h"
#include "multipole_operators.h"
#include "solver_alternator.h"
#include "tools.h"
#include "solver_basis.h"
#include "solver_hfb_broyden.h"
#include "solver_hfb_gradient.h"
#include "solver_ws.h"
#include "fragments.h"
#include "plot.h"

/** \file
 *  \brief Methods of the Action class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Constructor from a filename.
 */

Action::Action(const std::string &filename) : Action(DataTree(filename))
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constructor from a DataTree.
 */

Action::Action(const DataTree &_dataTree) : dataTree(_dataTree), state(_dataTree)
{
  DBG_ENTER;

  dataTree = DataTree::getDefault() + dataTree;

  dataTree.get(action           , "action"                  , true);
  dataTree.get(jobName          , "action/jobName"          , true);
  dataTree.get(saveResultFiles  , "action/saveResultFiles"  , true);
  dataTree.get(nbBlockingTries  , "action/nbBlockingTries"  , true);
  dataTree.get(basisOptimization, "action/basisOptimization", true);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Execute an action corresponding to the value of the Action::action property.
 */

void Action::run(void)
{
  DBG_ENTER;

  double startTime = Tools::clock();

  //=========================
  //== CLEAR THE PLOT PAGE ==
  //=========================

  if (useBokeh)
  {
    Plot::clear();
  }

  //=====================
  //== CALCULATE STUFF ==
  //=====================

  if      (action == "hfb"      ) calcHFB();
  else if (action == "ener"     ) calcEnergies();
  else if (action == "obs"      ) calcObservables();
  else if (action == "ws"       ) calcWS();
  else if (action == "mult_exp" ) calcMultipoleExpansion();
  else
  {
    ERROR("unknown action: " + action);
  }

  double length = Tools::clock() - startTime;

  Tools::mesg("Action", "Execution length (total): " + PF_GREEN("%.3f", length) + "s");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Solve the HFB equations while optimizing the basis deformation parameters.
 */

void Action::calcHFBquiet(void)
{
  DBG_ENTER;

  //==============================================
  //== OPTIMIZE THE BASIS PARAMETERS (optional) ==
  //==============================================

  if (basisOptimization)
  {
    Tools::mesg("ActHFB", PF_GREEN("=== Optimizing basis parameters ==="));

    SolverBasis solverBasis(dataTree, state);
    solverBasis.init();
    while (solverBasis.nextIter());
    state = solverBasis.state;
  }

  //=============================
  //== SOLVE THE HFB EQUATIONS ==
  //=============================

  Tools::mesg("ActHFB", PF_GREEN("=== Fixed basis HFB calculation ==="));

  SolverAlternator solverAlternator(dataTree, state);
  solverAlternator.init();
  while (solverAlternator.nextIter());
  state = solverAlternator.state;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Solve the HFB equations while blocking some state(s).
 */

void Action::calcHFBblocking(void)
{
  DBG_ENTER;

  //=========================
  //== CLEAR THE PLOT PAGE ==
  //=========================

  if (useBokeh)
  {
    Plot::clear();
  }

  //=============================================
  //== CHECK IF AN EVEN-EVEN SOLUTION IS GIVEN ==
  //=============================================

  bool isEvenNeut = state.sys.nNeut % 2 == 0;
  bool isEvenProt = state.sys.nProt % 2 == 0;

  INT evenNeut = isEvenNeut ? state.sys.nNeut : state.sys.nNeut - 1;
  INT evenProt = isEvenProt ? state.sys.nProt : state.sys.nProt - 1;

  bool evenNeutGiven = false;
  bool evenProtGiven = false;
  bool wantedNeutGiven = false;
  bool wantedProtGiven = false;

  INT givenNeut = -1;
  INT givenProt = -1;

  Tools::mesg("ActBlo", PF_GREEN("Target system: " + System(state.sys.nProt, state.sys.nNeut).info()));

  if (!state.empty())
  {
    MultipoleOperators multipoleOperators(state);

    // INFO(multipoleOperators.getNiceInfo());

    givenNeut = INT(multipoleOperators.nPart(NEUTRON) + 0.5);
    givenProt = INT(multipoleOperators.nPart(PROTON ) + 0.5);

    Tools::mesg("ActBlo", PF_GREEN("Given  system: " + System(givenProt, givenNeut).info()));

    INT closestEvenNeut = givenNeut % 2 == 0 ? givenNeut : givenNeut - 1;
    INT closestEvenProt = givenProt % 2 == 0 ? givenProt : givenProt - 1;

    evenNeutGiven   = givenNeut == closestEvenNeut;
    evenProtGiven   = givenProt == closestEvenProt;
    wantedNeutGiven = givenNeut == state.sys.nNeut;
    wantedProtGiven = givenProt == state.sys.nProt;
  }

  State startingState = state;

  if ((!wantedNeutGiven) || (!wantedProtGiven))
  {
    // calculate the even-even solution first
    INT nNeutOrig = state.sys.nNeut;
    INT nProtOrig = state.sys.nProt;

    state.sys.nProt = evenProt;
    state.sys.nNeut = evenNeut;

    if (evenNeutGiven && evenProtGiven)
    {
      Tools::mesg("ActBlo", PF_GREEN("Neighboring nuclear system given " + System(evenProt, evenNeut).info()));

      // generate the QP states by performing one HFB iteration
      Tools::mesg("ActBlo", "(Re)calculating the QP states and U/V matrices => performing iterations of the Broyden solver");
      DataTree dt = DataTree::getDefault() + dataTree;

      SolverHFBBroyden solverHFBBroyden(dt, state);
      solverHFBBroyden.init();
      //INFO(solverHFBBroyden.info());
      while(solverHFBBroyden.nextIter());

      startingState = solverHFBBroyden.state;
    }
    else
    {
      Tools::mesg("ActBlo", PF_GREEN("Neighboring nuclear system to be calculated " + System(evenProt, evenNeut).info()));

      calcHFBquiet();
      startingState = state;
    }
    startingState.sys.nProt = nProtOrig;
    startingState.sys.nNeut = nNeutOrig;
  }

  //========================================
  //== DISPLAY INFO ON THE STARTING STATE ==
  //========================================

  Tools::mesg("ActBlo", startingState.info());

  MultipoleOperators multipoleOperators(state);
  Tools::mesg("ActBlo", multipoleOperators.getNiceInfo());



  if (state.qpStates(NEUTRON).empty())
  {
    Tools::mesg("ActBlo", "(Re)calculating the QP states and U/V matrices => performing 1 iteration of the Broyden solver");
    DataTree dt = DataTree::getDefault() + dataTree;

    std::string interactionName;
    dataTree.get(interactionName, "interaction/name");
    dt.set("interaction/name", interactionName);

    SolverHFBBroyden solver(dt, startingState);
    solver.init();
    solver.nextIter();
    startingState = solver.state;
  }

  //=========================
  //== PRINT THE QP STATES ==
  //=========================

  state = startingState;

  state.qpStates(NEUTRON).sort("energy");
  state.qpStates(PROTON ).sort("energy");

  Tools::mesg("ActBlo", state.qpStates(NEUTRON).info(20, false));
  Tools::mesg("ActBlo", state.qpStates(PROTON ).info(20, false));

  //==========================================
  //== IDENTIFY THE QP STATES TO BE BLOCKED ==
  //==========================================

  std::list<StateId> blockingCandidatesn = {{-1, -1, NEUTRON}};
  std::list<StateId> blockingCandidatesp = {{-1, -1, PROTON }};

  if (!isEvenNeut) blockingCandidatesn = startingState.qpStates(NEUTRON).getFirstEmptyStates(nbBlockingTries);
  if (!isEvenProt) blockingCandidatesp = startingState.qpStates(PROTON ).getFirstEmptyStates(nbBlockingTries);

  //===================================
  //== DEACTIVATE BASIS OPTIMIZATION ==
  //===================================

  dataTree.set("action/basisOptimization", false);

  //=====================================
  //== DO THE BLOCKED HFB CALCULATIONS ==
  //=====================================

  State bestState = startingState;
  double bestTotalEnergy = 1e99;

  blockingTries.clear();
  INT id = 0;

  for (auto &bn : blockingCandidatesn)
  {
    for (auto &bp : blockingCandidatesp)
    {
      std::string blockingName = PF("blk%02d", id);

      state = bestState;
      state.blockedQPStates = {};
      if (bn.index >= 0) state.blockedQPStates.insert(bn);
      if (bp.index >= 0) state.blockedQPStates.insert(bp);

      Tools::mesg("ActBlo", "New blocking try: " + blockingName + " " + Tools::infoStr(state.blockedQPStates));

      INT maxIterTotal = 2000;
      dataTree.get(maxIterTotal, "solver/alternator/maxIter");

      SolverHFBBroyden solverHFBBroyden(dataTree, state);
      solverHFBBroyden.init();
      while (solverHFBBroyden.nextIter());
      state = solverHFBBroyden.state;

      if (solverHFBBroyden.status == Solver::CONVERGED)
      {
        if (state.totalEnergy < bestTotalEnergy)
        {
          bestTotalEnergy = state.totalEnergy;
          bestState = state;

          Tools::mesg("ActBlo", PF_GREEN("New best state found, ene: %9.6f MeV", bestTotalEnergy));
        }
      }

      std::string blockn = " ";
      std::string blockp = " ";

      if (bn.index >= 0) blockn = PF("(%d,%d/2)", bn.index, bn.omega * 2 + 1);
      if (bp.index >= 0) blockp = PF("(%d,%d/2)", bp.index, bp.omega * 2 + 1);

      blockingTries.push_back({id,
                               blockingName,
                               blockn,
                               blockp,
                               state.totalEnergy,
                               state.nbIter,
                               (solverHFBBroyden.status == Solver::CONVERGED)});

      Tools::mesg("ActBlo", getBlockingHistTable());

      id++;
    }
  }

  state = bestState;

  if (state.converged)
  {
    Tools::mesg("ActBlo", "End of blocking loop. " + PF_GREEN("Best energy found: %9.6e MeV", bestTotalEnergy));
  }
  else
  {
    Tools::mesg("ActBlo", "End of blocking loop. " + PF_RED("No solution found"));
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a table summarizing the blocking tries.
 */

const std::string Action::getBlockingHistTable(void) const
{
  DBG_ENTER;

  std::string result;
  result += TABLE_CENTER + TABLE_BLUE
            + "Tries" + TABLE_TD
            + TABLE_YELLOW
            + "Blocked Neut." + TABLE_TD
            + "Blocked Prot." + TABLE_TD
            + "Energy"        + TABLE_TD
            + "#iter"         + TABLE_TD
            + TABLE_TR;
  result += TABLE_CENTER + TABLE_YELLOW
            + " "          + TABLE_TD
            + "(id,omega)" + TABLE_TD
            + "(id,omega)" + TABLE_TD
            + "[MeV]"      + TABLE_TD
            + "[]"         + TABLE_TD
            + TABLE_TR;

  INT idMin = 0;
  double minEne = 1e99;

  for (auto &blockingTry: blockingTries)
  {
    if (!blockingTry.converged) continue;

    if (blockingTry.ene < minEne)
    {
      idMin = blockingTry.id;
      minEne = blockingTry.ene;
    }
  }

  for (auto &blockingTry: blockingTries)
  {
    std::string minimColor = TABLE_NORM;
    if (blockingTry.id == idMin) minimColor = TABLE_BLUE;


    if (blockingTry.id == idMin)
      result += TABLE_GREEN + TABLE_RIGHT + "* " + blockingTry.name + TABLE_TD;
    else
      result += TABLE_GREEN + TABLE_RIGHT + blockingTry.name + TABLE_TD;

    result += TABLE_CENTER + minimColor + blockingTry.blockn + TABLE_TD;
    result += TABLE_CENTER + minimColor + blockingTry.blockp + TABLE_TD;

    if (blockingTry.converged)
    {
      result += PF("%6.3f", blockingTry.ene) + TABLE_TD;
    }
    else
    {
      result += TABLE_RED + "not converged" + TABLE_TD;
    }

    result += TABLE_RIGHT + PF("%d", blockingTry.nbIter) + TABLE_TD + TABLE_TR;
  }

  DBG_RETURN(Tools::printTable(result));
}


//==============================================================================
//==============================================================================
//==============================================================================

/** Solve the HFB equations while optimizing the basis deformation parameters.
 */

void Action::calcHFB(void)
{
  DBG_ENTER;

  double startTime = Tools::clock();

  //==============================================
  //== CHECK IF THE TARGET NUCLEUS IS EVEN-EVEN ==
  //==============================================

  if ((state.sys.nNeut % 2 == 0) && (state.sys.nProt % 2 == 0))
  {
    Tools::mesg("ActHFB", PF_GREEN("Even-even target nucleus"));
    calcHFBquiet();
  }
  else
  {
    Tools::mesg("ActHFB", PF_GREEN("Odd-even or odd-odd target nuclear system "
                                   + state.sys.info()));
    calcHFBblocking();
  }

  //==========================
  //== CALCULATE PROPERTIES ==
  //==========================

  if (state.converged)
  {
    Tools::mesg("ActHFB", Tools::boxed(TABLE_GREEN + "Converged solution found ! \\o/ Properties:"));

    Tools::mesg("ActHFB", state.info());

    // calculate the energy contributions
    SolverHFBBroyden solverHFBBroyden(dataTree, state);
    solverHFBBroyden.init();
    solverHFBBroyden.interaction.calcEnergies();
    Tools::mesg("ActEne", solverHFBBroyden.interaction.getNiceInfo());

    // calculate and print some observables
    calcObservables();

    if (saveResultFiles)
    {
      DataTree finalDataTree = dataTree + state.getDataTree();
      finalDataTree.clean();

      finalDataTree.save(jobName + ".msg.gz");
    }
  }
  else
  {
    Tools::mesg("ActHFB", Tools::boxed(TABLE_RED + "Converged solution not found !"));

    if (saveResultFiles)
    {
      DataTree finalDataTree = dataTree + state.getDataTree();
      finalDataTree.clean();

      finalDataTree.save(jobName + ".msg.gz");
    }

    // calculate the energy contributions
    SolverHFBBroyden solverHFBBroyden(dataTree, state);
    solverHFBBroyden.init();
    solverHFBBroyden.interaction.calcEnergies();
    // Tools::mesg("ActEne", solverHFBBroyden.interaction.getNiceInfo());
  }

  state.calculationLength = Tools::clock() - startTime;

  Tools::mesg("Action", "Execution length (calcHFB): " + PF_GREEN("%.3f", state.calculationLength) + "s");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Solve a simplified form of the Woods-Saxon equations.
 *
 *  This can be used to create an initial state for the HFB solver.
 */

void Action::calcWS(void)
{
  DBG_ENTER;

  double startTime = Tools::clock();

  // Create a Woods-Saxon solver
  SolverWS solverWS(dataTree);

  // perform the WS calculation
  solverWS.init();
  while(solverWS.nextIter());
  state = solverWS.state;

  // Merge the solution with the initial DataTree instance
  dataTree.merge(state.getDataTree());

  // save the final state
  if (saveResultFiles)
    state.getDataTree().save(jobName + ".ws.msg.gz");

  state.calculationLength = Tools::clock() - startTime;

  Tools::mesg("Action", "Execution length (calcWS): " + PF_GREEN("%.3f", state.calculationLength) + "s");

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and print HFB field energy contributions.
 */

void Action::calcEnergies(void)
{
  DBG_ENTER;

  SolverHFBBroyden solverHFBBroyden(dataTree, state);
  solverHFBBroyden.init();
  solverHFBBroyden.interaction.calcEnergies();
  Tools::mesg("ActEne", solverHFBBroyden.interaction.getNiceInfo());

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and print some observables and geometrical properties.
 */

void Action::calcObservables(void)
{
  DBG_ENTER;

  if (!state.checkSolution()) DBG_LEAVE;

  //=======================
  //== INDIVIDUAL STATES ==
  //=======================

  // INFO(state.hfStates(NEUTRON).info());
  // INFO(state.hfStates(PROTON ).info());

  //=======================
  //== MULTIPOLE MOMENTS ==
  //=======================

  // Analytical calculation of the mass multipole moments
  MultipoleOperators multipoleOperators(state);
  Tools::mesg("ActObs", multipoleOperators.getNiceInfo());

  //=========================
  //== FRAGMENT PROPERTIES ==
  //=========================

  Fragments fragments(state);
  Tools::mesg("ActObs", fragments.getNiceInfo());

  //===================
  //== NUCLEAR RADII ==
  //===================

  Tools::mesg("ActObs", fragments.getNiceInfo("radii"));

  //===========================
  //== COLLECTIVE QUANTITIES ==
  //===========================

  // calculate inertia
  std::string interactionName;
  dataTree.get(interactionName, "interaction/name");
  state.calcInertia(interactionName);

  Tools::mesg("ActObs", state.getNiceInfo("inertia"));


  //===================
  //== THE QP STATES ==
  //===================

  // state.qpStates(NEUTRON).sort("energy");
  // state.qpStates(PROTON ).sort("energy");
  //
  // Tools::mesg("State ", state.qpStates(NEUTRON).info(-1, true));
  // Tools::mesg("State ", state.qpStates(PROTON ).info(-1, true));

  //===========
  //== BASIS ==
  //===========

  // Tools::mesg("ActObs", state.basis.info());

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate and export a multipole expansion of the density.
 */

void Action::calcMultipoleExpansion(void)
{
  DBG_ENTER;

  arma::vec rVals = arma::linspace(0., 20.0, 151);
  arma::cube result = Geometry::getDensityMultipoleExpansion(state, rVals);

  std::string resultStr;

  INT ir = 0;
  for (auto r: rVals)
  {
    std::string lineStr = PF("%e ", r);

    arma::vec valsNeut = result.slice(NEUTRON).row(ir);
    arma::vec valsProt = result.slice(PROTON ).row(ir);

    INT i = 0;
    for (auto v: valsNeut)
    {
      if (i % 2 == 0) lineStr += PF("%e ", v);
      i++;
    }

    i = 0;
    for (auto v: valsProt)
    {
      if (i % 2 == 0) lineStr += PF("%e ", v);
      i++;
    }

    // INFO(lineStr);
    resultStr += lineStr + "\n";

    ir++;
  }

  std::string fileName = jobName + ".mult_exp.csv";

  Tools::mesg("Action", "multipole expansion of the density saved in file " + fileName);
  Tools::save(resultStr, fileName);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Perform the DROP method.
 */


//==============================================================================
//==============================================================================
//==============================================================================

/** Perform the LINK method.
 */


//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Action::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += PF_GREEN("Action: ");
    result += dataTree.info(true);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"Action", ""},
      {"DataTree", dataTree.info(true)},
    });
  }

  DBG_RETURN(result);
}
