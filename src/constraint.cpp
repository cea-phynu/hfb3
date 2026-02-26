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

#include "constraint.h"
#include "tools.h"
#include "datatree.h"
#include "global.h"
#include "multipole_operators.h"

/** \file
 *  \brief Methods of the Constraint class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

/** Static initialization
 */

const IVEC Constraint::genderValue =
{
  MM, MM, MM, MM, MM, MM, MM, MM, SD, MA, QN,
  MM, MM, MM, MM, MM, MM, MM, MM, SD, MA, QN,
  MM, MM, MM, MM, MM, MM, MM, MM, SD, MA, QN
};

const std::vector<std::string> Constraint::types =
{
  "nNeut", "q00n", "q10n", "q20n", "q30n", "q40n", "q50n", "q60n", "sepDistn", "massAsymn", "qNeckn",
  "nProt", "q00p", "q10p", "q20p", "q30p", "q40p", "q50p", "q60p", "sepDistp", "massAsymp", "qNeckp",
  "nNucl", "q00t", "q10t", "q20t", "q30t", "q40t", "q50t", "q60t", "sepDistt", "massAsymt", "qNeckt"
};

// Using this to define the unit of measure as well. That is, the lambda value for the prefragments properties is their unit of measure, [fm^lambda].
const IVEC Constraint::lambdaValue =
{
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0,
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0,
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0
};

const IVEC Constraint::muValue =
{
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

const IVEC Constraint::isoValue =
{
  NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON, NEUTRON,
  PROTON, PROTON, PROTON, PROTON, PROTON, PROTON, PROTON, PROTON, PROTON, PROTON, PROTON,
  TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL, TOTAL
};

const arma::vec Constraint::factorValue =
{
  2 * sqrt(PI), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2 * sqrt(PI), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  2 * sqrt(PI), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

//==============================================================================
//==============================================================================
//==============================================================================

/** Dummy constructor
 */

Constraint::Constraint(void)
{
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Constructor from values.
 *
 *  \param _typeStr  Type of the constraint.
 *  \param _val      Value of the constraint.
 *  \param _lambda   Lambda of the constraint.
 */

Constraint::Constraint(const std::string &_typeStr, double _val, double _lagrangeMultiplier)
{
  DBG_ENTER;
  UINT _type = Tools::index(types, _typeStr);

  ASSERT(_type != -1, "Bad constraint type: '" + _typeStr + "'");

  type = INT(_type);
  typeStr = _typeStr;
  val = _val;
  lambda = lambdaValue[type];
  mu     = muValue[type];

  iso = isoValue[type];
  lagrangeMultiplier = _lagrangeMultiplier;
  measuredVal = val;
  gender = genderValue[type];
  factor = factorValue[type];

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Constraint::info(bool isShort) const
{
  std::string result = "";

  if (isShort)
  {
    result += Tools::treeStr(
    {
      {"type  ", types.at(type)},
      {"value ", Tools::infoStr(val)},
      {"measu.", Tools::infoStr(measuredVal)},
      {"lambda", Tools::infoStr(lambda)},
    }, isShort);
  }
  else
  {
    result += Tools::treeStr(
    {
      {"Constraint", ""},
      {"type  ", types.at(type)},
      {"value ", Tools::infoStr(val)},
      {"measu.", Tools::infoStr(measuredVal)},
      {"lambda", Tools::infoStr(lambda)},
    }, isShort);
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a DataTree instance.
 */

DataTree Constraint::getDataTree(void) const
{
  DBG_ENTER;

  DataTree result;

  result.set("constraints/" + types[type], val);
  result.set("state/constraints/lambda_" + types[type], lagrangeMultiplier);

  // remove beta2 constraints
  result.set("constraints/beta2p");
  result.set("constraints/beta2n");
  result.set("constraints/beta2t");

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Construct a list of Constraint from a DataTree instance.
 */

std::map<std::string, Constraint> Constraint::fromDataTree(const DataTree &dataTree)
{
  DBG_ENTER;

  std::map<std::string, Constraint> result;

  for (auto &type : types)
  {

    if (dataTree.contains("constraints/" + type))
    {
      double val = 0.0;

      if (dataTree.get(val, "constraints/" + type))
      {
        double lagrangeMultiplier = 0.0;
        dataTree.get(lagrangeMultiplier, "state/constraints/lambda_" + type, true);

        result[type] = Constraint(type, val, lagrangeMultiplier);
      }
    }
  }

  std::vector<std::string> equivTypes = {"beta2n", "beta2p", "beta2t"};

  for (auto &equivType : equivTypes)
  {
    if (dataTree.contains("constraints/" + equivType))
    {
      double val = 0.0;

      if (dataTree.get(val, "constraints/" + equivType))
      {
        std::string type = "";

        System system(dataTree);
        INT Z = system.nProt;
        INT N = system.nNeut;

        double A = double(Z) + double(N);

        double parts = -1.0;

        if (equivType == "beta2n") { parts = double(N); type = "q20n"; };
        if (equivType == "beta2p") { parts = double(Z); type = "q20p"; };
        if (equivType == "beta2t") { parts = double(A); type = "q20t"; };

        ASSERT(parts >= 0.0, "Wrong equivType ?");

        double equivValue = MultipoleOperators::getQ20FromBeta(parts, A, val);

        Tools::warning("Converting the " + equivType + PF(" constraint (val: %.3f)", val) + " to a " + type + PF(" constraint (val: %.3f)", equivValue));

        double lagrangeMultiplier = 0.0;
        dataTree.get(lagrangeMultiplier, "state/constraints/lambda_" + type, true);

        result[type] = Constraint(type, equivValue, lagrangeMultiplier);
      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the equivalent(?) deformation parameters for a Woods-Saxon calculation.
 */

arma::vec Constraint::toDefParameters(const std::vector<Constraint> &constraints, INT nNeut, INT nProt)
{
  DBG_ENTER;

  arma::vec result = {0.0};

  for (auto &c : constraints)
  {
    INFO(c.info());
  }

  // TODO : implement the magic
  INFO("Constraint::toDefParameters() TODO N:%f Z:%f", nNeut, nProt);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice std::string representation of an std::vector<Constraint>.
 */

const std::string Constraint::niceStr(const std::map<std::string, Constraint> &v, bool isShort)
{
  std::string result = "";

  std::vector<std::pair<std::string, std::string> > list;

  if (!isShort)
  {
    list.push_back({"Constraints", ""});
  }

  for (auto &c : v)
  {
    list.push_back({PF("%-6s", c.second.typeStr.c_str()), PF("%9.3e", c.second.val)});
  }

  if (!isShort)
  {
    for (auto &c : v)
    {
      list.push_back({PF("l%-5s", c.second.typeStr.c_str()), PF("%9.3e", c.second.lagrangeMultiplier)});
    }
  }

  result += Tools::treeStr(list, isShort);

  return result;
}
