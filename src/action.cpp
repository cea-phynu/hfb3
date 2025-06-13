
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

std::list<KeyStruct > Action::validKeys = {
    { "action"                  , "Action to be performed (can be: 'hfb_quiet', 'hfb', 'ener', 'obs', 'basis')", "hfb"    , "S" },
    { "action/basisOptimization", "Optimize the basis deformation parameters"                                  , "True"   , "B" },
    { "action/saveResultFiles"  , "Save result files"                                                          , "True"   , "B" },
    { "action/jobName"          , "Name of the job"                                                            , "unnamed", "S" },
    { "action/nbBlockingTrials" , "Number of blocking trials by isospin value"                                 , "4"      , "I" },
  };

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

  dataTree.get(action,            "action",                   true );
  dataTree.get(jobName,           "action/jobName",           true );
  dataTree.get(basisOptimization, "action/basisOptimization", true );
  dataTree.get(saveResultFiles,   "action/saveResultFiles",   true );
  dataTree.get(nbBlockingTrials,  "action/nbBlockingTrie",    true );

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

  if      (action == "hfb_quiet") calcHFBquiet();
  else if (action == "hfb"      ) calcHFB();
  else if (action == "ws"       ) calcWS();
  else if (action == "ws_hfb"   ) calcWSHFB();
  else if (action == "ener"     ) calcEnergies();
  else if (action == "obs"      ) calcObservables();
  else if (action == "basis"    ) calcBasis();
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

  state = SolverBasis::calcSingleHFB(dataTree, state, false, true);

  // Merge the solution with the initial DataTree instance
  dataTree = dataTree + state.getDataTree();

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

  if (!state.empty())
  {
    MultipoleOperators multipoleOperators(state);

    double nNeut = multipoleOperators.nPart(NEUTRON);
    double nProt = multipoleOperators.nPart(PROTON);

    double closestEvenNeut = (INT)(nNeut / 2.0 + 0.5) * 2.0;
    double closestEvenProt = (INT)(nProt / 2.0 + 0.5) * 2.0;

    evenNeutGiven = fabs(double(nNeut) - closestEvenNeut) < 1e-5;
    evenProtGiven = fabs(double(nProt) - closestEvenProt) < 1e-5;
    wantedNeutGiven = fabs(double(nNeut) - double(state.sys.nNeut)) < 1e-5;
    wantedProtGiven = fabs(double(nProt) - double(state.sys.nProt)) < 1e-5;
  }

  State startingState = state;

  if (wantedNeutGiven && wantedProtGiven)
  {
    Tools::mesg("ActBlo", PF_GREEN("Same nuclear system given " + state.sys.info()));
    Tools::mesg("ActBlo", PF("Generate the QP states by performing one HFB iteration"));

    SolverHFBBroyden solverHFBBroyden(dataTree, state);
    solverHFBBroyden.nextIter();
    startingState = solverHFBBroyden.state;
  }
  else
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
      Tools::mesg("ActBlo", PF("Generate the QP states by performing one HFB iteration"));
      SolverHFBBroyden solverHFBBroyden(dataTree, state);
      solverHFBBroyden.nextIter();

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

  MultipoleOperators multipoleOperators(state);
  Tools::mesg("ActBlo", multipoleOperators.getNiceInfo());

  // Tools::mesg("ActBlo", state.qpStates(NEUTRON).info());
  // Tools::mesg("ActBlo", state.qpStates(PROTON ).info());

  //==========================================
  //== IDENTIFY THE QP STATES TO BE BLOCKED ==
  //==========================================

  IVEC blockingCandidatesn = {-1};
  IVEC blockingCandidatesp = {-1};

  if (!isEvenNeut) blockingCandidatesn = startingState.qpStates(NEUTRON).getFirstEmptyStates(nbBlockingTrials);
  if (!isEvenProt) blockingCandidatesp = startingState.qpStates(PROTON ).getFirstEmptyStates(nbBlockingTrials);

  //===================================
  //== DEACTIVATE BASIS OPTIMIZATION ==
  //===================================

  dataTree.set("action/basisOptimization", false);

  //=====================================
  //== DO THE BLOCKED HFB CALCULATIONS ==
  //=====================================

  State bestState = startingState;
  double bestTotalEnergy = 1e99;

  blockingTrials.clear();
  INT id = 0;

  for (auto &bn : blockingCandidatesn)
  {
    for (auto &bp : blockingCandidatesp)
    {
      std::string blockingName = PF("blk%02d", id);

      bestState.blockedQP(NEUTRON) = bn;
      bestState.blockedQP(PROTON ) = bp;

      Tools::mesg("ActBlo", "New blocking try: " + blockingName);

      State state = SolverBasis::calcSingleHFB(dataTree, bestState, false, true);

      if (state.converged)
      {
        // === optional: print some information
        // Tools::mesg("ActBlo", state.qpStates(NEUTRON).info());
        // Tools::mesg("ActBlo", state.qpStates(PROTON ).info());

        // MultipoleOperators multipoleOperators(state);
        // Tools::mesg("ActObs", multipoleOperators.getNiceInfo());

        // for (auto &iso: {NEUTRON, PROTON})
        // {
        //   state.printCanonicalStatesInfo(state.rho(iso), iso);
        // }

        if (state.totalEnergy < bestTotalEnergy)
        {
          bestTotalEnergy = state.totalEnergy;
          bestState = state;

          Tools::mesg("ActBlo", PF_GREEN("New best state found, ene: %9.6f MeV", bestTotalEnergy));

        }
      }

      blockingTrials.push_back({id,
                               blockingName,
                               startingState.qpStates(NEUTRON).info(bn, 1, true),
                               startingState.qpStates(PROTON ).info(bp, 1, true),
                               state.totalEnergy,
                               state.nbIter,
                               state.converged});

      Tools::mesg("ActBlo", getBlockingHistTable());

      id++;
    }
  }

  Tools::mesg("ActBlo", PF("End of blocking loop. Best energy found: %9.6f MeV", bestTotalEnergy));
  state = bestState;

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
            + "Trials" + TABLE_TD
            + TABLE_YELLOW
            + "Blocked Neut." + TABLE_TD
            + "Blocked Prot." + TABLE_TD
            + "Energy"        + TABLE_TD
            + "#iter"         + TABLE_TD
            + TABLE_TR;
  result += TABLE_CENTER + TABLE_YELLOW
            + " "               + TABLE_TD
            + "(id,omega,ene.)" + TABLE_TD
            + "(id,omega,ene.)" + TABLE_TD
            + "[MeV]"           + TABLE_TD
            + "[]"              + TABLE_TD
            + TABLE_TR;

  INT idMin = 0;
  double minEne = 1e99;

  for (auto &blockingTry: blockingTrials)
  {
    if (!blockingTry.converged) continue;

    if (blockingTry.ene < minEne)
    {
      idMin = blockingTry.id;
      minEne = blockingTry.ene;
    }
  }

  for (auto &blockingTry: blockingTrials)
  {
    std::string minimColor = TABLE_NORM;
    if (blockingTry.id == idMin) minimColor = TABLE_BLUE;


    if (blockingTry.id == idMin)
      result += TABLE_GREEN + TABLE_RIGHT + "* " + blockingTry.name + TABLE_TD;
    else
      result += TABLE_GREEN + TABLE_RIGHT + blockingTry.name + TABLE_TD;

    result += TABLE_LEFT + minimColor + blockingTry.neut + TABLE_TD;
    result += TABLE_CENTER + blockingTry.prot + TABLE_TD;

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

    // calculate the energy contributions
    SolverHFBBroyden solverHFBBroyden(dataTree, state);
    solverHFBBroyden.interaction.calcEnergies();
    Tools::mesg("ActEne", solverHFBBroyden.interaction.getNiceInfo());


    // calculate and print some observables
    calcObservables();

    if (saveResultFiles)
    {
      DataTree finalDataTree = dataTree + state.getDataTree();
      finalDataTree.save(jobName + ".msg.gz");
    }
  }
  else
  {
    Tools::mesg("ActHFB", Tools::boxed(TABLE_RED + "Converged solution not found !"));

    // calculate the energy contributions
    SolverHFBBroyden solverHFBBroyden(dataTree, state);
    solverHFBBroyden.interaction.calcEnergies();
    Tools::mesg("ActEne", solverHFBBroyden.interaction.getNiceInfo());

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
  state = solverWS.calc();

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

/** Solve the simplified Woods-Saxon then the HFB equations.
 */

void Action::calcWSHFB(void)
{
  DBG_ENTER;

  double startTime = Tools::clock();

  // perform the WS calculation
  calcWS();

  // perform HFB starting from the WS state
  calcHFB();

  state.calculationLength = Tools::clock() - startTime;

  Tools::mesg("Action", "Execution length (calcWSHFB): " + PF_GREEN("%.3f", state.calculationLength) + "s");

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

  // calculate U and V matrices if needed
  state.calcUVFromRhoKappa(dataTree);

  // calculate inertia
  state.calcInertia();

  Tools::mesg("ActObs", state.getNiceInfo("inertia"));

  //===========
  //== BASIS ==
  //===========

  // Tools::mesg("ActObs", state.basis.info());

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Construct a Basis instance and print its characteristics.
 */

void Action::calcBasis(void)
{
  DBG_ENTER;

  Tools::mesg("ActBas", Basis(dataTree).info());

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
