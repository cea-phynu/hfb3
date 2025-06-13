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

#ifndef GEOMETRICAL_OPERATORS_H
#define GEOMETRICAL_OPERATORS_H

/** \file
 *  \brief Headers for the GeometricalOperators class.
 */

#include "generic.h"
#include "multi.h"
#include "state.h"
#include "multipole_operators.h"

/** \brief This class analytically calculates the geometrical operator (neck, separation distance, mass asymmetry) mean value.
 */

class GeometricalOperators : public Generic
{
public :

  // Constructors
  GeometricalOperators();
  GeometricalOperators(State);

  /// The initial State instance.
  State state;

  // Iz there a neck?
  double izNeck;

  // The position of the neck.
  double zNeck;

  /// A MultipoleOperators object.
  MultipoleOperators multipoleOperators;

  //============================================================================
  // void calcSepDistMatrix(void);   // TODO: Calculation of HO representation of the geometrical operator <d>
  // void calcMassAsymMatrix(void);  // TODO: Calculation of HO representation of the geometrical operator <xi>
  void calcQNeckMatrix(void);        // Calculation of HO representation of the geometrical operator <qNeck>

  // void calcSepDist(const arma::mat &rhon, const arma::mat &rhop);   // TODO: Calculation of the expectation value of the geometrical operator <d>
  // void calcMassAsym(const arma::mat &rhon, const arma::mat &rhop);  // TODO: Calculation of the expectation value of the geometrical operator <xi>
  void calcQneck(const arma::mat &rhon, const arma::mat &rhop);        // Calculation of the expectation value of the geometrical operator <qNeck>

  //============================================================================
  //============================================================================
  //============================================================================

  // /// The \f$\langle m,n,n_z,d,s|\hat{d}|m',n',n_z',d',s'\rangle\f$ matrices.
  // Multi<arma::mat, 1> sepdist0;
  //
  // /// The \f$\langle m,n,n_z,d,s|\hat{\xi}|m',n',n_z',d',s'\rangle\f$ matrices.
  // Multi<arma::mat, 1> massasym0;

  /// The \f$\langle m,n,n_z,d,s|\hat{Q}_N|m',n',n_z',d',s'\rangle\f$ matrices.
  Multi<arma::mat> qneck0;

  /// The \f$\psi|d|\psi\rangle\f$ value.
  // Multi<arma::mat, 2> matD;

  // /// The \f$\psi|\xi|\psi\rangle\f$ value.
  // Multi<arma::mat, 2> matXI;

  /// The \f$\psi|Q_N|\psi\rangle\f$ value.
  Multi<arma::mat> matQN;

  // /// The \f$\psi|\hat{d}|\psi\rangle\f$ value.
  // arma::mat sepdist;
  //
  // /// The \f$\psi|\hat{\xi}|\psi\rangle\f$ value.
  // arma::mat massasym;

  /// The \f$\psi|\hat{Q}_{N}|\psi\rangle\f$ value.
  arma::mat qneck;

  // /// The \f$\langle m,n|\hat{d}|m',n'\rangle\f$ matrices.
  // arma::mat matDk;

  // /// The \f$\langle m,n|\hat{\xi}|m',n'\rangle\f$ matrices.
  // arma::mat matXIk;

  /// The \f$\langle m,n|qNeck|m',n'\rangle\f$ matrices.
  arma::mat matQNk;

  /// The \f$\langle m,n|r^{0}|m',n'\rangle\f$ matrices.
  arma::mat matR0;

  /// The \f$\langle s|s'\rangle\f$ arma::mat.
  arma::mat matS0;

  //============================================================================
  //============================================================================
  //============================================================================

private:

  // void calcDk(void);
  // void calcXIk(void);
  void calcQNk(void);

};

#endif // GEOMETRICAL_OPERATORS_H
