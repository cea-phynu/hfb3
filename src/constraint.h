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

/** \brief Constraint to be used in a constrained HFB solver.
 *
 * This class represents a single constraint.
 */

class Constraint
{
public :

  /// Gender types.
  enum { MM, SD, MA, QN };

  Constraint(void);                                                    // #TEST#
  Constraint(const std::string &_typeStr,                              // #TEST#
             double _val,
             double _lambda = 0.0);

  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#

  DataTree getDataTree(void) const;                                    // #TEST#
  static std::map<std::string, Constraint> fromDataTree(               // #TEST#
    const DataTree &dataTree);
  static arma::vec toDefParameters(const std::vector<Constraint> &,    // #TEST#
                                   INT,
                                   INT);
  static const std::string niceStr(const std::map<std::string,         // #TEST#
                                   Constraint> &v,
                                   bool isShort = false);

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

  /// Possible types of constraints.
  static const std::vector<std::string> types;


  /// Type of constraint: multipole moment, geometrical, ...
  static const IVEC genderValue;

  /// Associated multipolar moment index for each constraint type
  static const IVEC lmValue;

  /// Associated isospin value for each constraint type
  static const IVEC isoValue;

  /// Associated scaling factor for each constraint type
  static const arma::vec factorValue;

  /// Should each constraint be used for inertia calculation ?
  static const std::vector<bool> useForInertiasValue;

  /// Type
  INT type = -1;

  /// Value.
  double val = 1e99;

  /// Lambda value.
  double lambda = 1e99;

  /// Measured value.
  double measuredVal = 1e99;

  /// Scaling factor.
  double factor = 1.0;

  /// Type.
  std::string typeStr = "";

  /// Multipole moment index
  INT lm = -1;

  /// Isospin
  INT iso = -1;

  /// Isospin
  INT gender = -1;

  /// Maximum number of constraint iterations
  INT itermax = -1;

  /// Should this constraint be used for inertia calculation ?
  bool useForInertias = false;
};

#endif // CONSTRAINT_H
