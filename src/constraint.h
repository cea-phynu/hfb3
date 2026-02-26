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

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

/** \file
 *  \brief Headers for the Constraint class.
 */

#include "global.h"
#include "datatree.h"

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define CONSTRAINT_VALID_KEYS \
{ "constraints/q00n"                  , "Constraint on <q00> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q10n"                  , "Constraint on <q10> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q20n"                  , "Constraint on <q20> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q30n"                  , "Constraint on <q30> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q40n"                  , "Constraint on <q40> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q50n"                  , "Constraint on <q50> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q60n"                  , "Constraint on <q60> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/q00p"                  , "Constraint on <q00> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q10p"                  , "Constraint on <q10> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q20p"                  , "Constraint on <q20> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q30p"                  , "Constraint on <q30> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q40p"                  , "Constraint on <q40> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q50p"                  , "Constraint on <q50> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q60p"                  , "Constraint on <q60> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/q00t"                  , "Constraint on <q00> value [total]"                                                 , ""   , "D" }, \
{ "constraints/q10t"                  , "Constraint on <q10> value [total]"                                                 , "0.0", "D" }, \
{ "constraints/q20t"                  , "Constraint on <q20> value [total]"                                                 , ""   , "D" }, \
{ "constraints/q30t"                  , "Constraint on <q30> value [total]"                                                 , ""   , "D" }, \
{ "constraints/q40t"                  , "Constraint on <q40> value [total]"                                                 , ""   , "D" }, \
{ "constraints/q50t"                  , "Constraint on <q50> value [total]"                                                 , ""   , "D" }, \
{ "constraints/q60t"                  , "Constraint on <q60> value [total]"                                                 , ""   , "D" }, \
{ "constraints/beta2p"                , "Constraint on <q60> value [prot.]"                                                 , ""   , "D" }, \
{ "constraints/beta2n"                , "Constraint on <q60> value [neut.]"                                                 , ""   , "D" }, \
{ "constraints/beta2t"                , "Constraint on <q60> value [total]"                                                 , ""   , "D" }, \
{ "constraints/sepDistn"              , "Constraint on the fragment separation distance [neut.]"                            , ""   , "D" }, \
{ "constraints/sepDistp"              , "Constraint on the fragment separation distance [prot.]"                            , ""   , "D" }, \
{ "constraints/sepDistt"              , "Constraint on the fragment separation distance [total]"                            , ""   , "D" }, \
{ "constraints/massAsymn"             , "Constraint on the fragment mass asymmetry [neut.]"                                 , ""   , "D" }, \
{ "constraints/massAsymp"             , "Constraint on the fragment mass asymmetry [prot.]"                                 , ""   , "D" }, \
{ "constraints/massAsymt"             , "Constraint on the fragment mass asymmetry [total]"                                 , ""   , "D" }, \
{ "constraints/qNeckn"                , "Constraint on <Qneck> value [prot.]"                                               , ""   , "D" }, \
{ "constraints/qNeckp"                , "Constraint on <Qneck> value [neut.]"                                               , ""   , "D" }, \
{ "constraints/qNeckt"                , "Constraint on <Qneck> value [total]"                                               , ""   , "D" }, \
{ "constraints/nNeut"                 , "Constraint on the number of neutrons"                                              , ""   , "D" }, \
{ "constraints/nProt"                 , "Constraint on the number of protons"                                               , ""   , "D" }, \
{ "constraints/nNucl"                 , "Constraint on the number of nucleons"                                              , ""   , "D" }, \
{ "state/constraints/lambda_q00n"     , "Lagrange multiplier for the constraint on <q00> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q10n"     , "Lagrange multiplier for the constraint on <q10> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q20n"     , "Lagrange multiplier for the constraint on <q20> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q30n"     , "Lagrange multiplier for the constraint on <q30> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q40n"     , "Lagrange multiplier for the constraint on <q40> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q50n"     , "Lagrange multiplier for the constraint on <q50> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q60n"     , "Lagrange multiplier for the constraint on <q60> value [neut.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q00p"     , "Lagrange multiplier for the constraint on <q00> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q10p"     , "Lagrange multiplier for the constraint on <q10> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q20p"     , "Lagrange multiplier for the constraint on <q20> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q30p"     , "Lagrange multiplier for the constraint on <q30> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q40p"     , "Lagrange multiplier for the constraint on <q40> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q50p"     , "Lagrange multiplier for the constraint on <q50> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q60p"     , "Lagrange multiplier for the constraint on <q60> value [prot.]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q00t"     , "Lagrange multiplier for the constraint on <q00> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q10t"     , "Lagrange multiplier for the constraint on <q10> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q20t"     , "Lagrange multiplier for the constraint on <q20> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q30t"     , "Lagrange multiplier for the constraint on <q30> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q40t"     , "Lagrange multiplier for the constraint on <q40> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q50t"     , "Lagrange multiplier for the constraint on <q50> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_q60t"     , "Lagrange multiplier for the constraint on <q60> value [total]"                     , ""   , "D" }, \
{ "state/constraints/lambda_sepDistn" , "Lagrange multiplier for the constraint on the fragment separation distance [neut.]", ""   , "D" }, \
{ "state/constraints/lambda_sepDistp" , "Lagrange multiplier for the constraint on the fragment separation distance [prot.]", ""   , "D" }, \
{ "state/constraints/lambda_sepDistt" , "Lagrange multiplier for the constraint on the fragment separation distance [total]", ""   , "D" }, \
{ "state/constraints/lambda_massAsymn", "Lagrange multiplier for the constraint on the fragment mass asymmetry [neut.]"     , ""   , "D" }, \
{ "state/constraints/lambda_massAsymp", "Lagrange multiplier for the constraint on the fragment mass asymmetry [prot.]"     , ""   , "D" }, \
{ "state/constraints/lambda_massAsymt", "Lagrange multiplier for the constraint on the fragment mass asymmetry [total]"     , ""   , "D" }, \
{ "state/constraints/lambda_qNeckn"   , "Lagrange multiplier for the constraint on <Qneck> value [prot.]"                   , ""   , "D" }, \
{ "state/constraints/lambda_qNeckp"   , "Lagrange multiplier for the constraint on <Qneck> value [neut.]"                   , ""   , "D" }, \
{ "state/constraints/lambda_qNeckt"   , "Lagrange multiplier for the constraint on <Qneck> value [total]"                   , ""   , "D" }, \
{ "state/constraints/lambda_nNeut"    , "Lagrange multiplier for the constraint on the number of neutrons"                  , ""   , "D" }, \
{ "state/constraints/lambda_nProt"    , "Lagrange multiplier for the constraint on the number of protons"                   , ""   , "D" }, \
{ "state/constraints/lambda_nNucl"    , "Lagrange multiplier for the constraint on the number of nucleons"                  , ""   , "D" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Constraint to be used in a constrained HFB solver.
 *
 * This class represents a single constraint.
 */

class Constraint
{
public:

  /// Gender types.
  enum { MM, SD, MA, QN };

  Constraint(void);                                                    // #TEST#
  Constraint(const std::string &_typeStr,                              // #TEST#
             double _val = 0.0,
             double _lagrangeMultiplier = 0.0);

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  DataTree getDataTree(void) const;                                    // #TEST#
  static std::map<std::string, Constraint> fromDataTree(const DataTree &dataTree); // #TEST#
  static arma::vec toDefParameters(const std::vector<Constraint> &,    // #TEST#
                                   INT,
                                   INT);
  static const std::string niceStr(const std::map<std::string, Constraint> &v,     // #TEST#
                                   bool isShort = false);

  //============================================================================
  //============================================================================
  //============================================================================

  /// Possible types of constraints.
  static const std::vector<std::string> types;

  /// Type of constraint: multipole moment, geometrical, ...
  static const IVEC genderValue;

  /// Associated multipolar moment lambda value for each constraint type
  static const IVEC lambdaValue;

  /// Associated multipolar moment mu value for each constraint type
  static const IVEC muValue;

  /// Associated isospin value for each constraint type
  static const IVEC isoValue;

  /// Associated scaling factor for each constraint type
  static const arma::vec factorValue;

  /// Type
  INT type = -1;

  /// Value
  double val = 1e99;

  /// Lagrange multiplier value
  double lagrangeMultiplier = 1e99;

  /// Measured value
  double measuredVal = 1e99;

  /// Scaling factor
  double factor = 1.0;

  /// Type
  std::string typeStr = "";

  /// Multipole moment lambda value
  INT lambda = -1;

  /// Multipole moment mu value
  INT mu = -1;

  /// Isospin
  INT iso = -1;

  /// Isospin
  INT gender = -1;

  /// Maximum number of constraint iterations
  INT itermax = -1;
};

#endif // CONSTRAINT_H
