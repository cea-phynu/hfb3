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

#include "cliparser.h"
#include "global.h"
#include "tools.h"
#include "datatree.h"

/** \file
 *  \brief Methods for the CliParser class.
 */

//==============================================================================
//==============================================================================
//==============================================================================
//
/** The constructor with CLI arguments to be parsed.
 *  CLI example: ./hfb3 toto.rhl.dat.gz '{action: hfb, HFB/Z: 100, HFB/N: 160}'
 */

CliParser::CliParser(INT argc, char **argv)
{
  if (argc < 2) help();

  bool readSaveTo = false;
  msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING}; // default value for CLI use.

  for (INT i = 1; i < argc; i++)
  {
    std::string arg = argv[i];

    if (readSaveTo)
    {
      readSaveTo = false;
      saveTo = arg;
    }
    else if ((arg[0] != '-') or (arg.size() == 1))
    {
      // '-' for std::cin
      fileToLoad.push_back(arg);
    }
    else
    {
      if      (arg == "-h" || arg == "--help")  help();
      else if (arg == "-o")                     readSaveTo = true;
      else if (arg == "-f" || arg == "--force") forceInvalidDataTree = true;
      else if (arg == "-g" || arg == "--bokeh") useBokeh = true;

      else if (arg == "--list-keys") printKeys();
      else if (arg == "--list-interactions") printInteractions();
      else if (arg == "--version")   version();

      else if (arg == "-t0")        tableStyle = 0;
      else if (arg == "-t1")        tableStyle = 1;
      else if (arg == "-t2")        tableStyle = 2;
      else if (arg == "-l")         showLogo = true;
      else if (arg == "--no-utf8")  useUtf8 = false;
      else if (arg == "--no-color") useColors = false;

      else if (arg == "-v0")        msgToOut = {MSG_ERROR, MSG_MAIN};
      else if (arg == "-v1")        msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO};
      else if (arg == "-v2")        msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING}; // default value
      else if (arg == "-v3")        msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_TIME};
      else if (arg == "-v4")        msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_TIME, MSG_DEBUG};
      else if (arg == "-v5")        msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_TIME, MSG_DEBUG, MSG_TIMELINE};
      else if (arg == "-v")         msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_TIME, MSG_DEBUG, MSG_TIMELINE};

      else ERROR("unknown option: '" + arg + "'");
    }
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print the version.
 */

void CliParser::version(void)
{
  Tools::mesg("Version", Tools::version());
  Tools::mesg("Authors", CFG_AUTHORS);
  Tools::mesg("Compil.", std::string(SKILL));
  Tools::end(0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print some help.
 */

void CliParser::help(void)
{
  std::cout << "Usage: hfb3 [OPTION|FILE|DICT]..." << std::endl;
  std::cout << "> Next-gen HFB solver." << std::endl;
  std::cout << std::endl;
  std::cout << "Possible OPTIONS:" << std::endl;
  std::cout << "  -h, --help           print this help" << std::endl;
  std::cout << "  --version            print the version of the program" << std::endl;
  std::cout << "  --list-keys          print the list of possible keys" << std::endl;
  std::cout << "  --list-interactions  print the list of possible interactions" << std::endl;
  std::cout << "  -o FILE.msg.gz       save the DataTree in FILE.msg.gz" << std::endl;
  std::cout << "  -f, --force          force the loading of an invalid DataTree" << std::endl;
  std::cout << "  -g, --bokeh          enable bokeh plots" << std::endl;
  std::cout << std::endl;
  std::cout << "  -v[0...6]    verbosity level (default 1)" << std::endl;
  std::cout << "  -t[0...2]    style of the printed tables" << std::endl;
  std::cout << "  --no-colors  no colored output" << std::endl;
  std::cout << "  --no-utf8    do not use utf8 characters in the output" << std::endl;
  std::cout << "  -l           show the HFB3 logo at start" << std::endl;
  std::cout << std::endl;
  std::cout << "Possible FILES:" << std::endl;
  std::cout << "  FILE  any file (supported formats: hfb3[.gz], msg[.gz], rhl.dat[.gz], theo[.gz])" << std::endl;
  std::cout << "  -     read content from standard input (can be .gz content)" << std::endl;
  std::cout << std::endl;
  std::cout << "DICT format:" << std::endl;
  std::cout << "  {key1:val1, key2:val2, ...}" << std::endl;
  std::cout << std::endl;
  std::cout << "Examples:" << std::endl;
  std::cout << "  * launch an HFB calculation using the file 'examples/16O_groundstate.hfb3' with a constraint on <Q20> (30fm) and 50 HFB iterations max.:" << std::endl;
  std::cout << "  bin/hfb3 examples/16O_groundstate.hfb3 '{constraints/q20t:30.0,solver/alternator/maxIter:50}'" << std::endl;
  Tools::end(0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print the list of possible interactions.
 */

void CliParser::printInteractions(void)
{
  INFO("Possible interactions: D1S, D1M, D1N, D3G3, D12S, D2, DG. Add 'X' for exact Coulomb (for example 'D1S' -> 'D1SX').");

  Tools::end(0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Print the list of possible keys.
 */

void CliParser::printKeys(void)
{
  std::string table;

  table +=
    TABLE_YELLOW + TABLE_LEFT      + "Key"           + TABLE_TD
  + TABLE_YELLOW                   + "Description"   + TABLE_TD
  + TABLE_YELLOW                   + "Default value" + TABLE_TD
  + TABLE_YELLOW                   + "Type"          + TABLE_TD + TABLE_TR;

  DataTree dt;

  std::map<std::string, std::string> fullTypeName = {
    {"B", "Boolean" },
    {"I", "Int" },
    {"D", "Double" },
    {"S", "String" },
    {"V", "Vec" },
    {"M", "Mat" },
    {"C", "Cube" },
    {"MV", "Multi Vec" },
    {"MM", "Multi Mat" },
    {"MC", "Multi Cube" },
  };

  for (auto &o : general.validKeys)
  {
    table +=
      TABLE_NORM  + o.key                + TABLE_TD
    + TABLE_GREEN + o.description        + TABLE_TD
    + TABLE_BLUE  + o.defaultValue + " " + TABLE_TD
    + TABLE_NORM  + fullTypeName[o.type] + TABLE_TD + TABLE_TR;
  }

  Tools::info(Tools::printTable(table));
  Tools::info("* " + PF_BLUE("state/") + " keys change from one HFB iteration to another");

  Tools::end(0);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string CliParser::info(bool isShort) const
{
  std::string result = "";

  result += Tools::treeStr(
  {
    {"CliParser", ""},
    {"saveTo", saveTo},
    {"toLoad", Tools::infoStr(fileToLoad)},
  }, isShort);

  return result;
}


