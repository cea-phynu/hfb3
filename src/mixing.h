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

#ifndef MIXING_H
#define MIXING_H

/** \file
 *  \brief Headers for the Mixing class.
 *
 *  Created on July 25th 2019
 *  Mikael Frosini
 */

#include "global.h"
#include "generic.h"
#include "state.h"

// REF : https://arxiv.org/pdf/0805.4446.pdf
// https://www.sciencedirect.com/science/article/pii/S0021999196900595

// TODO Make a proper benchmark
// TODO Fine-tune w0 and size
// TODO Flatten matrices block-wise to avoid useless zeros

//==============================================================================
//==============================================================================
//==============================================================================

/// Implementation of the Broyden Mixing method.
class Mixing : public Generic
{
public:

  Mixing();                                                            // #TEST#
  Mixing(const DataTree &);                                            // #TEST#
  Mixing(INT _mixingSize, double _linearFactor, double _mixingFactor); // #TEST#

  bool newVec(const arma::vec &);                                      // #TEST#
  const arma::vec & getVec();                                          // #TEST#
  const std::string info(bool isShort = USE_SHORT_INFO) const;         // #TEST#
  bool isLinearMode(void);

  //============================================================================
  //== Externally set parameters ===============================================
  //============================================================================

  /// Number of configurations.
  INT                mixingSize;

  /// Linear mixing factor.
  double             linearFactor;

  /// General mixing factor.
  double             mixingFactor;

  //============================================================================
  //============================================================================
  //============================================================================

  /// List of keys used by this class.
  static std::list<KeyStruct > validKeys;

private:
  INT                m = 0;     // m <= M gives back linear mixing
  INT                n = -1; // density matrices have size n * n
  INT                N = -1; // n*n
  INT                n_const = 0;
  double w0 = 0.1;     // Parameter preventing singular systems. maybe to be put into
  arma::vec x;     //
  arma::vec F;     //
  arma::mat dx;     // \Delta x(l) = x(l+1) - x(l) <-> density matrices   : (4N, M)
  arma::mat dF;     // F(l)  = y(l) - x(l)                                : (4N, M)
  arma::mat dFdF;     // <F(l)-F(l-i)|F(l)-F(l-j)>                        : (M, M)
  arma::vec dFF;     // <F(l)-F(l-i)|F(l)>                                : (M)
  arma::vec gamma;     // State to the equation
  bool firstBroyden = true;
  bool firstLinear  = true;
};

#endif // MIXING_H
