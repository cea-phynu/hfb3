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

#include "interaction.h"
#include "state.h"
#include <memory>
#include "field_kinetic.h"
#include "field_cdm2.h"
#include "field_central.h"
#include "field_coulomb_slater.h"
#include "field_coulomb_exact.h"
#include "field_density.h"
#include "field_density_fr.h"
#include "field_rearrangement.h"
#include "field_rearrangement_fr.h"
#include "field_spin_orbit.h"
#include "field_ws.h"

/** \file
 *  \brief Methods of the Interaction class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

std::list<KeyStruct > Interaction::validKeys =
  {
    { "interaction/name"              , "Interaction name (use --list-interactions to list them)"              , "D1S", "S" },
    { "interaction/cdm2"              , "Use a 2-body center-of-mass correction field"                         , ""   , "B" },
    { "interaction/central"           , "Include a central field"                                              , ""   , "B" },
    { "interaction/central/p"         , "Central field: p parameter (width)"                                   , ""   , "V" },
    { "interaction/central/w"         , "Central field: w parameter"                                           , ""   , "V" },
    { "interaction/central/b"         , "Central field: b parameter"                                           , ""   , "V" },
    { "interaction/central/h"         , "Central field: h parameter"                                           , ""   , "V" },
    { "interaction/central/m"         , "Central field: m parameter"                                           , ""   , "V" },
    { "interaction/coulombExact"      , "Include a coulomb field (exact treatment)"                            , ""   , "B" },
    { "interaction/coulombSlater"     , "Include a coulomb field (slater approximation)"                       , ""   , "B" },
    { "interaction/density"           , "Include a zero-range density-dependent field"                         , ""   , "B" },
    { "interaction/density/a"         , "Zero-range density-dependent field: a parameter"                      , ""   , "V" },
    { "interaction/density/x"         , "Zero-range density-dependent field: x parameter"                      , ""   , "V" },
    { "interaction/density/t"         , "Zero-range density-dependent field: t parameter"                      , ""   , "V" },
    { "interaction/densityFR"         , "Include a finite-range density-dependent field"                       , ""   , "B" },
    { "interaction/densityFR/p"       , "Finite-range density-dependent field: p parameter (width)"            , ""   , "V" },
    { "interaction/densityFR/w"       , "Finite-range density-dependent field: w parameter"                    , ""   , "V" },
    { "interaction/densityFR/b"       , "Finite-range density-dependent field: b parameter"                    , ""   , "V" },
    { "interaction/densityFR/h"       , "Finite-range density-dependent field: h parameter"                    , ""   , "V" },
    { "interaction/densityFR/m"       , "Finite-range density-dependent field: m parameter"                    , ""   , "V" },
    { "interaction/kinetic"           , "Include a kinetic field"                                              , ""   , "B" },
    { "interaction/rearrangement"     , "Include a zero-range rearrangement field"                             , ""   , "B" },
    { "interaction/rearrangement/a"   , "Zero-range rearrangement field: a parameter"                          , ""   , "V" },
    { "interaction/rearrangement/x"   , "Zero-range rearrangement field: x parameter"                          , ""   , "V" },
    { "interaction/rearrangement/t"   , "Zero-range rearrangement field: t parameter"                          , ""   , "V" },
    { "interaction/rearrangementFR"   , "Include a finite-range rearrangement field"                           , ""   , "B" },
    { "interaction/rearrangementFR/p" , "Finite-range rearrangement field: p parameter (width)"                , ""   , "V" },
    { "interaction/rearrangementFR/w" , "Finite-range rearrangement field: w parameter"                        , ""   , "V" },
    { "interaction/rearrangementFR/b" , "Finite-range rearrangement field: b parameter"                        , ""   , "V" },
    { "interaction/rearrangementFR/h" , "Finite-range rearrangement field: h parameter"                        , ""   , "V" },
    { "interaction/rearrangementFR/m" , "Finite-range rearrangement field: m parameter"                        , ""   , "V" },
    { "interaction/spinOrbit"         , "Include a zero-range spin-orbit field"                                , ""   , "B" },
    { "interaction/spinOrbit/wso"     , "Zero-range spin-orbit field: wso parameter"                           , ""   , "V" },
    { "interaction/spinOrbitFR/wso"   , "Finite-range spin-orbit field: wso parameter"                         , ""   , "V" },
    { "interaction/woodsSaxon"        , "Include a Woods-Saxon field"                                          , ""   , "B" },
  };

//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Constructor from a pointer to a State instance and a DataTree instance.
 */

Interaction::Interaction(const std::string &_interactionName, State *_state) :
  Interaction(Interaction::getInteractionDataTree(_interactionName), _state)
{
  DBG_ENTER;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Constructor from a pointer to a State instance and a DataTree instance.
 */

Interaction::Interaction(DataTree _dataTree, State *_state) :
  dataTree(_dataTree),
  state(_state)
{
  DBG_ENTER;

  interactionName = dataTree.getS("interaction/name", false);

  DataTree interactionDataTree = getInteractionDataTree(interactionName);

  dataTree = interactionDataTree + dataTree;

  if (dataTree.getB("interaction/kinetic"        , false)) addField("kinetic"        );
  if (dataTree.getB("interaction/central"        , false)) addField("central"        );
  if (dataTree.getB("interaction/coulombSlater"  , false)) addField("coulombSlater"  );
  if (dataTree.getB("interaction/density"        , false)) addField("density"        );
  if (dataTree.getB("interaction/spinOrbit"      , false)) addField("spinOrbit"      );
  if (dataTree.getB("interaction/cdm2"           , false)) addField("cdm2"           );
  if (dataTree.getB("interaction/rearrangement"  , false)) addField("rearrangement"  );
  if (dataTree.getB("interaction/densityFR"      , false)) addField("densityFR"      );
  if (dataTree.getB("interaction/rearrangementFR", false)) addField("rearrangementFR");
  if (dataTree.getB("interaction/coulombExact"   , false)) addField("coulombExact"   );
  if (dataTree.getB("interaction/woodsSaxon"     , false)) addField("woodsSaxon"     );

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Extract field parameters from a DataTree.
 */

std::vector<Field::Parameters> Interaction::getParametersFromDataTree(DataTree &dataTree, const std::string &fieldName)
{
  DBG_ENTER;

  auto parameterKeys = dataTree.getKeys("interaction/" + fieldName + "/");

  // find the needed number of this field instances
  arma::vec p;
  for (auto &k: parameterKeys)
  {
    p = dataTree.getV(k, false);
    break; // only one parameter needed here
  }

  UINT nbInstances = p.empty() ? 1 : p.n_elem;
  // INFO("field %s: %d instances", fieldName.c_str(), nbInstances);

  std::vector<Field::Parameters> listOfParameters;

  for (UINT id = 0; id < nbInstances; id++)
  {
    arma::vec p;
    Field::Parameters parameters;

    for (auto &k: parameterKeys)
    {
      p = dataTree.getV(k, false);

      UINT vlast = k.find_last_of("/");
      std::string parameterName = k.substr(vlast + 1);

      parameters[parameterName] = p[id];
      // INFO("parameter: %s val: %f", parameterName.c_str(), p[id]);
    }

    listOfParameters.push_back(parameters);
  }

  DBG_RETURN(listOfParameters);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Add a new field.
 */

void Interaction::addField(const std::string &fieldName)
{
  DBG_ENTER;

  auto parameters = getParametersFromDataTree(dataTree, fieldName);

  // for each instance, construct the parameters struct and instantiate the Field
  for (UINT id = 0; id < parameters.size(); id++)
  {
    switch (hash2(fieldName)) {
      case "kinetic"_h        : fieldsList.push_back(std::make_shared<FieldKinetic        >(FieldKinetic        (parameters[id], state))); break;
      case "central"_h        : fieldsList.push_back(std::make_shared<FieldCentral        >(FieldCentral        (parameters[id], state))); break;
      case "coulombSlater"_h  : fieldsList.push_back(std::make_shared<FieldCoulombSlater  >(FieldCoulombSlater  (parameters[id], state))); break;
      case "density"_h        : fieldsList.push_back(std::make_shared<FieldDensity        >(FieldDensity        (parameters[id], state))); break;
      case "spinOrbit"_h      : fieldsList.push_back(std::make_shared<FieldSpinOrbit      >(FieldSpinOrbit      (parameters[id], state))); break;
      case "cdm2"_h           : fieldsList.push_back(std::make_shared<FieldCDM2           >(FieldCDM2           (parameters[id], state))); break;
      case "rearrangement"_h  : fieldsList.push_back(std::make_shared<FieldRearrangement  >(FieldRearrangement  (parameters[id], state))); break;
      case "densityFR"_h      : fieldsList.push_back(std::make_shared<FieldDensityFR      >(FieldDensityFR      (parameters[id], state))); break;
      case "rearrangementFR"_h: fieldsList.push_back(std::make_shared<FieldRearrangementFR>(FieldRearrangementFR(parameters[id], state))); break;
      case "coulombExact"_h   : fieldsList.push_back(std::make_shared<FieldCoulombExact   >(FieldCoulombExact   (parameters[id], state))); break;
      case "woodsSaxon"_h     : fieldsList.push_back(std::make_shared<FieldWS             >(FieldWS             (parameters[id], state))); break;
      default                 : INFO("unknown field type: " + fieldName); break;
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the interaction.
 */

void Interaction::calcFields(bool ignoreMissingRhoKappa)
{
  DBG_ENTER;

  if (!ignoreMissingRhoKappa)
  {
    ASSERT(state->rho.contains(NEUTRON  ), "missing rho(NEUTRON) matrix.");
    ASSERT(state->rho.contains(PROTON   ), "missing rho(PROTON) matrix.");
    ASSERT(state->kappa.contains(NEUTRON), "missing kappa(NEUTRON) matrix.");
    ASSERT(state->kappa.contains(PROTON ), "missing kappa(PROTON) matrix.");

    ASSERT(!state->rho(NEUTRON  ).empty(), "empty rho(NEUTRON) matrix");
    ASSERT(!state->rho(PROTON   ).empty(), "empty rho(PROTON) matrix");
    ASSERT(!state->kappa(NEUTRON).empty(), "empty kappa(NEUTRON) matrix");
    ASSERT(!state->kappa(PROTON ).empty(), "empty kappa(PROTON) matrix");
  }

  double startTime = Tools::clock();

  // Calculate the interaction
  for (auto &f : fieldsList)
  {
    f->calcField();
  }

  calcLength = Tools::clock() - startTime;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Clear the interaction.
 */

void Interaction::clear(void)
{
  DBG_ENTER;

  for (auto f : fieldsList)
  {
    f->clear();
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Calculate the energy contributions.
 */

void Interaction::calcEnergies(void)
{
  DBG_ENTER;

  totalEnergy = arma::zeros(2);

  calcFields();

  for (auto &f: fieldsList)
  {
    f->calcEnergy();

    if (!f->contributeToEnergy) continue;

    for (INT iso : {NEUTRON, PROTON})
    {
      for (INT type : {Field::DIRECT, Field::EXCHANGE, Field::PAIRING})
      {
        if (!f->energy.contains(iso, type)) continue;

        totalEnergy(iso) += f->energy(iso, type);
      }
    }
  }

  // Store energy contributions.
  energyContributions = getEnergyContributions();

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the energy contributions.
 */

const std::map<std::string, double> Interaction::getEnergyContributions(void) const
{
  DBG_ENTER;

  std::map<std::string, double> result;

  for (auto &f : fieldsList)
  {
    for (INT iso : {NEUTRON, PROTON})
    {
      for (INT type : {Field::DIRECT, Field::EXCHANGE, Field::PAIRING})
      {
        if (!f->energy.contains(iso, type)) continue;

        std::string name =
          f->shortName
        + "_"
        + Tools::strIsospin(iso)
        + "_"
        + Field::typeStr.at(type);

        if (result.count(name) == 0)
          result[name] = f->energy(iso, type);
        else
          result[name] += f->energy(iso, type);

      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/**
 * Accessor field of type "name"
 */

std::shared_ptr<Field> Interaction::operator()(std::string name, INT id)
{
  DBG_ENTER;

  // INFO("new try name='%s' id=%d", name.c_str(), id);

  INT idc = 0;
  for (auto f: fieldsList)
  {
    if (f->name == name)
    {
      if (id == idc)
      {
        DBG_RETURN(f);
      }

      idc++;
    }
  }

  std::string names;
  bool ifirst = true;

  for (auto n: fieldsList)
  {
    if (!ifirst) names += ',';
    names += n->name;
    ifirst = false;
  }

  ERROR(PF("No field with the name '%s' and index %d. Existing fields: [%s]", name.c_str(), id, names.c_str()));

  // must not happen
  DBG_RETURN(fieldsList[0]);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return the total hamiltonian.
 */

const arma::mat Interaction::getHamiltonianContributions(INT iso, INT type)
{
  DBG_ENTER;

  arma::mat result = arma::zeros(state->basis.HOqn.nb, state->basis.HOqn.nb);

  for (auto &f : fieldsList)
  {
    if (!f->field.contains(iso, type)) continue;

    ASSERT(!f->field(iso, type).empty(), "Empty field: " + f->name);

    result += f->field(iso, type);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return a nice table displaying the energy contributions.
 */

const std::string Interaction::getNiceInfo(void) const
{
  DBG_ENTER;

  std::list<std::list<std::string> > table;

  arma::vec total = arma::zeros(3);

  for (auto &field : fieldsList)
  {
    for (auto type: {Field::DIRECT, Field::EXCHANGE, Field::PAIRING}) {
      std::string fname = field->name;

      std::list<std::string> line;

      double pEnergy = 0.;
      double nEnergy = 0.;
      double tEnergy = 0.;

      // Skip over contributions that dont exist
      if (!field->energy.contains(PROTON, type)) continue;

      nEnergy = field->energy(NEUTRON, type);
      pEnergy = field->energy(PROTON , type);
      tEnergy = nEnergy + pEnergy;

      // Add to the totals
      if (field->contributeToEnergy) {
        total(NEUTRON) += nEnergy;
        total(PROTON ) += pEnergy;
        total(TOTAL) += tEnergy;
      }
      else
      {
        // in case of a field non-contributing to the total energy
        fname += " (non-contrib.)";
      }

      // Create the line
      line.push_back(fname + " " + Field::typeStr[type]);
      line.push_back(PF("%9.3f", nEnergy));
      line.push_back(PF("%9.3f", pEnergy));
      line.push_back(PF("%9.3f", tEnergy));

      table.push_back(line);

    }
  }

  // Total energies
  std::list<std::string> values;
  values.push_back(TABLE_RED + "TOTAL (HFB Energy)");
  values.push_back(TABLE_RED + PF("%7.3f", total(NEUTRON)));
  values.push_back(TABLE_RED + PF("%7.3f", total(PROTON )));
  values.push_back(TABLE_RED + PF("%7.3f", total(TOTAL  )));
  table.push_back(values);

  std::string result = Tools::valueTable(
    "Energies",
    {"Neutron", "Proton", "Neut+Prot"},
    {"[MeV]"  , "[MeV]" , "[MeV]"    },
    table
  );

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

/** Return some info.
 */

const std::string Interaction::info(bool isShort) const
{
  DBG_ENTER;

  std::string result = "";

  if (isShort)
  {
    result += interactionName + " (";
    bool ifirst = true;

    for (auto &f : fieldsList)
    {
      if (ifirst) ifirst = false;
      else result += ",";

      result += PF_YELLOW(f->shortName);
    }

    result += ")";
  }
  else
  {
    std::vector<std::pair<std::string, std::string> > list = {{"Interaction", ""}};

    for (auto &f : fieldsList)
    {
      list.push_back({f->shortName, f->info(true)});
    }

    result += Tools::treeStr(list, isShort);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const DataTree Interaction::getInteractionDataTree(const std::string& interactionName)
{
  DBG_ENTER;

  DataTree noInteractionDataTree = DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
cdm2                False
central             False
coulombExact        False
coulombSlater       False
density             False
densityFR           False
kinetic             False
rearrangement       False
rearrangementFR     False
spinOrbit           False
spinOrbitFR         False
tensorFR            False
woodsSaxon          False
)toto");

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "") DBG_RETURN(noInteractionDataTree);

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "WS") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                WS
kinetic             True
woodsSaxon          True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1S") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1S
cdm2                True
central             True
coulombSlater       True
density             True
kinetic             True
rearrangement       True
spinOrbit           True
central/p           (    0.700,    1.200)
central/w           (-1720.300,  103.639)
central/b           ( 1300.000, -163.483)
central/h           (-1813.530,  162.812)
central/m           ( 1397.600, -223.934)
density/a           (0.333333)
density/x           (1.000)
density/t           (1390.600)
rearrangement/a     (0.333333)
rearrangement/x     (1.000)
rearrangement/t     (1390.600)
spinOrbit/wso       (130.000)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1M") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1M
cdm2                True
central             True
coulombSlater       True
density             True
kinetic             True
rearrangement       True
spinOrbit           True
central/p           (     0.5000000,    1.0000000)
central/w           (-12769.9005956,  488.7353877)
central/b           ( 14021.2751562, -751.9060996)
central/h           (-15099.3874268,  671.2598083)
central/m           ( 11944.6321518, -692.5926358)
density/a           (0.333333)
density/x           (1.000)
density/t           (1562.8351962)
rearrangement/a     (0.333333)
rearrangement/x     (1.000)
rearrangement/t     (1562.8351962)
spinOrbit/wso       (115.4633181)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1N") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1N
cdm2                True
central             True
coulombSlater       True
density             True
kinetic             True
rearrangement       True
spinOrbit           True
central/p           (    0.800,    1.200)
central/w           (-2047.610,  293.020)
central/b           ( 1700.000, -300.780)
central/h           (-2414.930,  414.590)
central/m           ( 1519.350, -316.840)
density/a           (    0.333333)
density/x           (    1.000)
density/t           ( 1609.460)
rearrangement/a     (    0.333333)
rearrangement/x     (    1.000)
rearrangement/t     ( 1609.460)
spinOrbit/wso       (  115.000)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D3G3") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D3G3
cdm2                True
central             True
coulombSlater       True
density             True
kinetic             True
rearrangement       True
spinOrbit           True
central/p           (     0.47018,     0.74857,   1.96704)
central/w           ( -7543.80481,   590.47032,   4.62570)
central/b           ( 13485.56724, -1751.39902,  -9.25140)
central/h           (-14708.98783,  1582.84160,   9.25140)
central/m           (  6669.46303,  -909.26152, -18.50280)
density/a           (0.33333333333333333)
density/x           (1.000)
density/t           (1400.000)
rearrangement/a     (0.33333333333333333)
rearrangement/x     (1.000)
rearrangement/t     (1400.000)
spinOrbit/wso       (115.1361)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D12S") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D12S
cdm2                True
central             True
coulombSlater       True
densityFR           True
kinetic             True
rearrangementFR     True
spinOrbit           True
central/p           (    0.700,    1.200)
central/w           (-1720.300,  103.639)
central/b           ( 1300.000, -163.483)
central/h           (-1813.530,  162.812)
central/m           ( 1397.600, -223.934)
densityFR/p         (0.00001)
densityFR/w         (1390.600)
densityFR/b         (1390.600)
densityFR/h         (0.000000)
densityFR/m         (0.000000)
densityFR/a         (0.33333333333333333)
rearrangementFR/p   (0.00001)
rearrangementFR/w   (1390.600)
rearrangementFR/b   (1390.600)
rearrangementFR/h   (0.000000)
rearrangementFR/m   (0.000000)
rearrangementFR/a   (0.33333333333333333)
spinOrbit/wso       (130.000)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D2") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D2
cdm2                True
central             True
coulombSlater       True
densityFR           True
kinetic             True
rearrangementFR     True
spinOrbit           True
central/p           (    0.800,    1.300)
central/w           (-1176.440,   93.741)
central/b           (  800.000, -162.161)
central/h           ( -927.366,  122.414)
central/m           ( 1115.573, -223.859)
densityFR/p         (   0.600)
densityFR/w         (1800.000)
densityFR/b         ( 600.000)
densityFR/h         ( 400.000)
densityFR/m         (-600.000)
densityFR/a         (0.33333333333333333)
rearrangementFR/p   (   0.600)
rearrangementFR/w   (1800.000)
rearrangementFR/b   ( 600.000)
rearrangementFR/h   ( 400.000)
rearrangementFR/m   (-600.000)
rearrangementFR/a   (0.33333333333333333)
spinOrbit/wso       (130.000)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "DG") DBG_RETURN(noInteractionDataTree + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                DG
cdm2                True
central             True
coulombSlater       True
densityFR           True
kinetic             True
rearrangementFR     True
spinOrbitFR         True
tensorFR            True
central/p           (    0.800,     1.24)
central/w           (-1190.016,  109.179)
central/b           (  800.000, -191.226)
central/h           ( -877.422,  133.441)
central/m           ( 1198.923, -277.509)
densityFR/p         (    0.600)
densityFR/w         ( 1836.200)
densityFR/b         (  581.600)
densityFR/h         (  377.600)
densityFR/m         ( -633.220)
densityFR/a         (    0.33333333333333333)
rearrangementFR/p   (    0.600)
rearrangementFR/w   ( 1836.200)
rearrangementFR/b   (  581.600)
rearrangementFR/h   (  377.600)
rearrangementFR/m   ( -633.220)
rearrangementFR/a   (    0.33333333333333333)
spinOrbitFR/p       (    0.200)
spinOrbitFR/w       (  146.483)
spinOrbitFR/h       (   29.634)
tensorFR/p          (    1.100)
tensorFR/w          ( -392.544)
tensorFR/h          ( -196.481)
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1SX") DBG_RETURN(getInteractionDataTree("D1S") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1SX
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1NX") DBG_RETURN(getInteractionDataTree("D1N") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1NX
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D1MX") DBG_RETURN(getInteractionDataTree("D1M") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D1MX
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D3G3X") DBG_RETURN(getInteractionDataTree("D3G3") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D3G3X
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D12SX") DBG_RETURN(getInteractionDataTree("D12S") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D12SX
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "D2X") DBG_RETURN(getInteractionDataTree("D2") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                D2X
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  if (interactionName == "DGX") DBG_RETURN(getInteractionDataTree("DG") + DataTree::fromContent(
R"toto(#!HFB3!#

[interaction]
name                DGX
coulombSlater       False
coulombExact        True
)toto"));

  //============================================================================
  //============================================================================
  //============================================================================

  ERROR("Invalid interaction type");

  DBG_RETURN(DataTree());
};

