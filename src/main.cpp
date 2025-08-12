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

#include "hfb3.h"
#include "solver_hfb_broyden.h"

/** \file
 *  \brief Main program.
 */

/*
 *    +---+   +---+  +----------+ +-------+       +-----------+
 *   /   /|  /   /|/           /|/       /+---+  /           /|
 *  +---+ | +---+ +-----------+ +-------+/   /| +-----------+ |
 *  |   | +-|   | |           | |       +---+ | |           | |
 *  |   |/  |   | |           |/|           | | |           | |
 *  |   +---+   | |   ----+---+ +   -----   | | +---+----   | |
 *  |           | |       | /   |           | |   +-|       | |
 *  |           | |       |/    |           | |  /  |       | |
 *  |   +---+   | |   +---+     |   -----   | + +---+----   | |
 *  |   | + |   | |   | +       |           |/  |           | +
 *  |   |/  |   |/|   |/        |       +---+   |           |/
 *  +---+   +---+ +---+         +-------+       +-----------+
 *
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief The main routine.
 *
 *  Let's go.
 */

INT main(INT argc, char **argv)
{
  // ===== parse the CLI options =====
  CliParser cliParser(argc, argv);

  // ===== Welcome message =====
  if (cliParser.showLogo) Tools::mesg("Logo!!", Tools::getLogoStr());
  Tools::mesg("Versio", Tools::version());
  Tools::mesg("Compil", std::string(SKILL));

#ifdef NO_DBG_STACK
  Tools::warning("The debug stack has been disabled.\nRecompile without `-DNO_DBG_STACK` to enable it.");
#endif

  DBG_ENTER;

  // ===== Initial DataTree =====
  auto dataTree = DataTree::getDefault();

  for (auto &file : cliParser.fileToLoad)
  {
    try
    {
      std::string content;

      if (file == "-")
        content = Tools::readCin();
      else
      {
        if (file[0] == '{') content = file;
        else content = Tools::readFile(file);
      }

      DataTree newDataTree = DataTree().fromContent(content);
      dataTree.merge(newDataTree);
    }
    catch (const std::string &msg)
    {
      ERROR(msg);
    }
  }

  // ===== Update the general constants/flags =====
  general = General(dataTree);

  // ===== Validate the initial dataTree =====
  dataTree.validate();

  // ===== Print the initial dataTree =====
  // INFO(dataTree.info());

  // ===== Save the initial dataTree ? =====
  if (cliParser.saveTo != "")
  {
    // dataTree.save(cliParser.saveTo);

    // State state(dataTree);
    // IOberger ioberger;
    // ioberger.saveState(state, cliParser.saveTo);

    State state(dataTree);
    state.getDataTree().save(cliParser.saveTo);

    Tools::end();
  }

  // ===== Action ! =====
  std::string actionStr;
  dataTree.get(actionStr, "action", true);

  ASSERT(actionStr != "none", "No `action` key specified.");

  if (actionStr == "")
  {
    INFO("Missing action value. Exiting.");
  }
  else if (actionStr == "test")
  {
    // insert some testing code here, get it executed with 'hfb3 "{action:test}"'
    auto solverAlternator = SolverAlternator(dataTree);

    // Tools::mesg("_main_", solverAlternator.info());

    solverAlternator.init();
    while(solverAlternator.nextIter());
    solverAlternator.finalize();
  }
  else if (actionStr == "info")
  {
    INFO(dataTree.info());
    INFO(State(dataTree).info());
    INFO(Basis(dataTree).info());
    // INFO(State(dataTree).getOmegaContributionsInfo());
  }
  else
  {
    // Create an Action instance
    Action action(dataTree);

    // Execute action
    action.run();
  }

  // Optionally save the plot page
  Plot::save("result.html");

  DBG_RETURN(0);
}

