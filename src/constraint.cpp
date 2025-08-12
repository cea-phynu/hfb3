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

std::list<KeyStruct > Constraint::validKeys =
  {
    { "constraints/q00n", "Constraint on <q00> value [neut.]", ""   , "D" },
    { "constraints/q10n", "Constraint on <q10> value [neut.]", ""   , "D" },
    { "constraints/q20n", "Constraint on <q20> value [neut.]", ""   , "D" },
    { "constraints/q30n", "Constraint on <q30> value [neut.]", ""   , "D" },
    { "constraints/q40n", "Constraint on <q40> value [neut.]", ""   , "D" },
    { "constraints/q50n", "Constraint on <q50> value [neut.]", ""   , "D" },
    { "constraints/q60n", "Constraint on <q60> value [neut.]", ""   , "D" },
    { "constraints/q00p", "Constraint on <q00> value [prot.]", ""   , "D" },
    { "constraints/q10p", "Constraint on <q10> value [prot.]", ""   , "D" },
    { "constraints/q20p", "Constraint on <q20> value [prot.]", ""   , "D" },
    { "constraints/q30p", "Constraint on <q30> value [prot.]", ""   , "D" },
    { "constraints/q40p", "Constraint on <q40> value [prot.]", ""   , "D" },
    { "constraints/q50p", "Constraint on <q50> value [prot.]", ""   , "D" },
    { "constraints/q60p", "Constraint on <q60> value [prot.]", ""   , "D" },
    { "constraints/q00t", "Constraint on <q00> value [total]", ""   , "D" },
    { "constraints/q10t", "Constraint on <q10> value [total]", "0.0", "D" },
    { "constraints/q20t", "Constraint on <q20> value [total]", ""   , "D" },
    { "constraints/q30t", "Constraint on <q30> value [total]", ""   , "D" },
    { "constraints/q40t", "Constraint on <q40> value [total]", ""   , "D" },
    { "constraints/q50t", "Constraint on <q50> value [total]", ""   , "D" },
    { "constraints/q60t", "Constraint on <q60> value [total]", ""   , "D" },
    { "constraints/beta2p", "Constraint on <q60> value [prot.]", ""   , "D" },
    { "constraints/beta2n", "Constraint on <q60> value [neut.]", ""   , "D" },
    { "constraints/beta2t", "Constraint on <q60> value [total]", ""   , "D" },

    { "constraints/sepDistn" , "Constraint on the fragment separation distance [neut.]", "", "D" },
    { "constraints/sepDistp" , "Constraint on the fragment separation distance [prot.]", "", "D" },
    { "constraints/sepDistt" , "Constraint on the fragment separation distance [total]", "", "D" },
    { "constraints/massAsymn", "Constraint on the fragment mass asymmetry [neut.]"     , "", "D" },
    { "constraints/massAsymp", "Constraint on the fragment mass asymmetry [prot.]"     , "", "D" },
    { "constraints/massAsymt", "Constraint on the fragment mass asymmetry [total]"     , "", "D" },
    { "constraints/qNeckn"   , "Constraint on <Qneck> value [prot.]"                   , "", "D" },
    { "constraints/qNeckp"   , "Constraint on <Qneck> value [neut.]"                   , "", "D" },
    { "constraints/qNeckt"   , "Constraint on <Qneck> value [total]"                   , "", "D" },
    { "constraints/nNeut"    , "Constraint on the number of neutrons"                  , "", "D" },
    { "constraints/nProt"    , "Constraint on the number of protons"                   , "", "D" },
    { "constraints/nNucl"    , "Constraint on the number of nucleons"                  , "", "D" },

    { "state/constraints/lambda_q00n", "Lagrange multiplier for the constraint on <q00> value [neut.]", "", "D" },
    { "state/constraints/lambda_q10n", "Lagrange multiplier for the constraint on <q10> value [neut.]", "", "D" },
    { "state/constraints/lambda_q20n", "Lagrange multiplier for the constraint on <q20> value [neut.]", "", "D" },
    { "state/constraints/lambda_q30n", "Lagrange multiplier for the constraint on <q30> value [neut.]", "", "D" },
    { "state/constraints/lambda_q40n", "Lagrange multiplier for the constraint on <q40> value [neut.]", "", "D" },
    { "state/constraints/lambda_q50n", "Lagrange multiplier for the constraint on <q50> value [neut.]", "", "D" },
    { "state/constraints/lambda_q60n", "Lagrange multiplier for the constraint on <q60> value [neut.]", "", "D" },
    { "state/constraints/lambda_q00p", "Lagrange multiplier for the constraint on <q00> value [prot.]", "", "D" },
    { "state/constraints/lambda_q10p", "Lagrange multiplier for the constraint on <q10> value [prot.]", "", "D" },
    { "state/constraints/lambda_q20p", "Lagrange multiplier for the constraint on <q20> value [prot.]", "", "D" },
    { "state/constraints/lambda_q30p", "Lagrange multiplier for the constraint on <q30> value [prot.]", "", "D" },
    { "state/constraints/lambda_q40p", "Lagrange multiplier for the constraint on <q40> value [prot.]", "", "D" },
    { "state/constraints/lambda_q50p", "Lagrange multiplier for the constraint on <q50> value [prot.]", "", "D" },
    { "state/constraints/lambda_q60p", "Lagrange multiplier for the constraint on <q60> value [prot.]", "", "D" },
    { "state/constraints/lambda_q00t", "Lagrange multiplier for the constraint on <q00> value [total]", "", "D" },
    { "state/constraints/lambda_q10t", "Lagrange multiplier for the constraint on <q10> value [total]", "", "D" },
    { "state/constraints/lambda_q20t", "Lagrange multiplier for the constraint on <q20> value [total]", "", "D" },
    { "state/constraints/lambda_q30t", "Lagrange multiplier for the constraint on <q30> value [total]", "", "D" },
    { "state/constraints/lambda_q40t", "Lagrange multiplier for the constraint on <q40> value [total]", "", "D" },
    { "state/constraints/lambda_q50t", "Lagrange multiplier for the constraint on <q50> value [total]", "", "D" },
    { "state/constraints/lambda_q60t", "Lagrange multiplier for the constraint on <q60> value [total]", "", "D" },

    { "state/constraints/lambda_sepDistn" , "Lagrange multiplier for the constraint on the fragment separation distance [neut.]", "", "D" },
    { "state/constraints/lambda_sepDistp" , "Lagrange multiplier for the constraint on the fragment separation distance [prot.]", "", "D" },
    { "state/constraints/lambda_sepDistt" , "Lagrange multiplier for the constraint on the fragment separation distance [total]", "", "D" },
    { "state/constraints/lambda_massAsymn", "Lagrange multiplier for the constraint on the fragment mass asymmetry [neut.]"     , "", "D" },
    { "state/constraints/lambda_massAsymp", "Lagrange multiplier for the constraint on the fragment mass asymmetry [prot.]"     , "", "D" },
    { "state/constraints/lambda_massAsymt", "Lagrange multiplier for the constraint on the fragment mass asymmetry [total]"     , "", "D" },
    { "state/constraints/lambda_qNeckn"   , "Lagrange multiplier for the constraint on <Qneck> value [prot.]"                   , "", "D" },
    { "state/constraints/lambda_qNeckp"   , "Lagrange multiplier for the constraint on <Qneck> value [neut.]"                   , "", "D" },
    { "state/constraints/lambda_qNeckt"   , "Lagrange multiplier for the constraint on <Qneck> value [total]"                   , "", "D" },
    { "state/constraints/lambda_nNeut"    , "Lagrange multiplier for the constraint on the number of neutrons"                  , "", "D" },
    { "state/constraints/lambda_nProt"    , "Lagrange multiplier for the constraint on the number of protons"                   , "", "D" },
    { "state/constraints/lambda_nNucl"    , "Lagrange multiplier for the constraint on the number of nucleons"                  , "", "D" },
  };

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

// Using this to define the unit of measure as well. That is, the lm value for the prefragments properties is their unit of measure, [fm^lm].
const IVEC Constraint::lmValue =
{
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0,
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0,
  0, 0, 1, 2, 3, 4, 5, 6, 1, 0, 0
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

const std::vector<bool> Constraint::useForInertiasValue =
{
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false, false, false, false, false, false, false, false, false,
  false, false, false,  true,  true,  true,  true,  true, false, false, false
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

Constraint::Constraint(const std::string &_typeStr, double _val, double _lambda)
{
  DBG_ENTER;
  UINT _type = Tools::index(types, _typeStr);

  ASSERT(_type != -1, "Bad constraint type: '" + _typeStr + "'");

  type = INT(_type);
  typeStr = _typeStr;
  val = _val;
  lm = lmValue[type];

  iso = isoValue[type];
  lambda = _lambda;
  measuredVal = val;
  gender = genderValue[type];
  factor = factorValue[type];
  useForInertias = useForInertiasValue[type];

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
  result.set("state/constraints/lambda_" + types[type], lambda);
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
        double lambda = 0.0;
        dataTree.get(lambda, "state/constraints/lambda_" + type, true);

        result[type] = Constraint(type, val, lambda);
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

        double lambda = 0.0;
        dataTree.get(lambda, "state/constraints/lambda_" + type, true);

        result[type] = Constraint(type, equivValue, lambda);
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
    list.push_back({PF("%-6s", c.first.c_str()), PF("%9.3e", c.second.val)});
  }

  if (!isShort)
  {
    for (auto &c : v)
    {
      list.push_back({PF("l%-5s", c.first.c_str()), PF("%9.3e", c.second.lambda)});
    }
  }

  result += Tools::treeStr(list, isShort);

  return result;
}
