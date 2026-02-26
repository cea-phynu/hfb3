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

#ifndef INTERACTION_H
#define INTERACTION_H

/** \file
 *  \brief Headers for the Interaction class.
 */

#include <memory>
#include "global.h"
#include "generic.h"
#include "field.h"

class State;

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define INTERACTION_VALID_KEYS \
{ "interaction/name"              , "Interaction name (use --list-interactions to list them)", "D1S", "S" }, \
{ "interaction/cdm2"              , "Use a 2-body center-of-mass correction field"           , "", "B" }, \
{ "interaction/central"           , "Include a central field"                                , "", "B" }, \
{ "interaction/central/p"         , "Central field: p parameter"                             , "", "V" }, \
{ "interaction/central/w"         , "Central field: w parameter"                             , "", "V" }, \
{ "interaction/central/b"         , "Central field: b parameter"                             , "", "V" }, \
{ "interaction/central/h"         , "Central field: h parameter"                             , "", "V" }, \
{ "interaction/central/m"         , "Central field: m parameter"                             , "", "V" }, \
{ "interaction/coulombExact"      , "Include a coulomb field (exact treatment)"              , "", "B" }, \
{ "interaction/coulombSlater"     , "Include a coulomb field (slater approximation)"         , "", "B" }, \
{ "interaction/density"           , "Include a zero-range density-dependent field"           , "", "B" }, \
{ "interaction/density/a"         , "Zero-range density-dependent field: a parameter"        , "", "V" }, \
{ "interaction/density/x"         , "Zero-range density-dependent field: x parameter"        , "", "V" }, \
{ "interaction/density/t"         , "Zero-range density-dependent field: t parameter"        , "", "V" }, \
{ "interaction/densityFR"         , "Include a finite-range density-dependent field"         , "", "B" }, \
{ "interaction/densityFR/a"       , "Finite-range density-dependent field: a parameter"      , "", "V" }, \
{ "interaction/densityFR/p"       , "Finite-range density-dependent field: p parameter"      , "", "V" }, \
{ "interaction/densityFR/w"       , "Finite-range density-dependent field: w parameter"      , "", "V" }, \
{ "interaction/densityFR/b"       , "Finite-range density-dependent field: b parameter"      , "", "V" }, \
{ "interaction/densityFR/h"       , "Finite-range density-dependent field: h parameter"      , "", "V" }, \
{ "interaction/densityFR/m"       , "Finite-range density-dependent field: m parameter"      , "", "V" }, \
{ "interaction/kinetic"           , "Include a kinetic field"                                , "", "B" }, \
{ "interaction/rearrangement"     , "Include a zero-range rearrangement field"               , "", "B" }, \
{ "interaction/rearrangement/a"   , "Zero-range rearrangement field: a parameter"            , "", "V" }, \
{ "interaction/rearrangement/x"   , "Zero-range rearrangement field: x parameter"            , "", "V" }, \
{ "interaction/rearrangement/t"   , "Zero-range rearrangement field: t parameter"            , "", "V" }, \
{ "interaction/rearrangementFR"   , "Include a finite-range rearrangement field"             , "", "B" }, \
{ "interaction/rearrangementFR/a" , "Finite-range rearrangement field: a parameter"          , "", "V" }, \
{ "interaction/rearrangementFR/p" , "Finite-range rearrangement field: p parameter"          , "", "V" }, \
{ "interaction/rearrangementFR/w" , "Finite-range rearrangement field: w parameter"          , "", "V" }, \
{ "interaction/rearrangementFR/b" , "Finite-range rearrangement field: b parameter"          , "", "V" }, \
{ "interaction/rearrangementFR/h" , "Finite-range rearrangement field: h parameter"          , "", "V" }, \
{ "interaction/rearrangementFR/m" , "Finite-range rearrangement field: m parameter"          , "", "V" }, \
{ "interaction/spinOrbit"         , "Include a zero-range spin-orbit field"                  , "", "B" }, \
{ "interaction/spinOrbit/wso"     , "Zero-range spin-orbit field: wso parameter"             , "", "V" }, \
{ "interaction/spinOrbitFR"       , "Include a finite-range spin-orbit field"                , "", "B" }, \
{ "interaction/spinOrbitFR/p"     , "Finite-range spin-orbit field: p parameter"             , "", "V" }, \
{ "interaction/spinOrbitFR/w"     , "Finite-range spin-orbit field: w parameter"             , "", "V" }, \
{ "interaction/spinOrbitFR/h"     , "Finite-range spin-orbit field: h parameter"             , "", "V" }, \
{ "interaction/tensorFR"          , "Include a finite-range tensor field"                    , "", "B" }, \
{ "interaction/tensorFR/p"        , "Finite-range tensor field: p parameter"                 , "", "V" }, \
{ "interaction/tensorFR/w"        , "Finite-range tensor field: w parameter"                 , "", "V" }, \
{ "interaction/tensorFR/h"        , "Finite-range tensor field: h parameter"                 , "", "V" }, \
{ "interaction/woodsSaxon"        , "Include a Woods-Saxon field"                            , "", "B" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Store the Field instances.
 *
 * This class stores the instances of the Field* classes.
 */

class Interaction : public Generic
{
public :

  Interaction(const std::string &_interactionName = "D1S", State *_state = NULL); // #TEST#
  Interaction(DataTree dataTree, State *_state = NULL);                // #TEST#
  void clear(void);                                                    // #TEST#
  void calcFields(bool ignoreMissingRhoKappa = false);                 // #TEST#
  void calcEnergies(void);                                             // #TEST#
  const std::map<std::string, double> getEnergyContributions(void) const; // #TEST#
  const arma::mat getHamiltonianContributions(INT iso, INT type);      // #TEST#
  static const DataTree getInteractionDataTree(const std::string& interactionName);

  std::shared_ptr<Field> operator()(std::string name, INT id = 0);

  void addField(const std::string &fieldName);                         // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  const std::string getNiceInfo(void) const;                           // #TEST#
  std::vector<Field::Parameters> getParametersFromDataTree(DataTree &dataTree, const std::string &fieldName);

  //============================================================================
  //============================================================================
  //============================================================================

  /// The interaction name.
  std::string interactionName = "";

  /// A DataTree instance
  DataTree dataTree;

  /// A pointer to a State instance.
  State *state;

  /// Total energies.
  arma::vec totalEnergy;

  /// List of Field instances.
  std::vector<std::shared_ptr<Field>> fieldsList;

  /// Table of energy contributions.
  std::map<std::string, double> energyContributions;

  /// Total calculation length.
  double calcLength = 0.0;
};

#endif // INTERACTION_H
