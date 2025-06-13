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

#ifndef PLOT_H
#define PLOT_H

/** \file
 *  \brief Headers for the Plot class.
 */

#include "global.h"
#include "generic.h"

class Mesh;

/// A class containing static plotting functions (not a namespace because swig ignores namespaces).
class Plot : public Generic
{
public:
  /// Set the active slot.
  static void slot(INT id);                                           // #TEST#

  /// Clear the plotting instance.
  static void clear(const std::string &name = "");                    // #TEST#

  /// Save the plot.
  static void save(const std::string filename);                       // #TEST#

  /// Plot xy values.
  static void curve(const std::string &xlabel,
                    const std::string &ylabel,
                    double x,
                    double y);

  /// Plot temporal values.
  static void curve(const std::string &xlabel,
                    const std::string &ylabel,
                    double y);

  /// Plot a curve.
  static void curve(const std::string &xlabel,
                    const std::string &ylabel,
                    const arma::vec &x,
                    const arma::vec &y);

  /// Plot a generic field.
  static void field(const std::string &title,
                    const arma::mat &mat,
                    const Mesh &mesh,
                    const std::string &xlabel = "z [fm]",
                    const std::string &ylabel = "r [fm]");


  /// Plot a local density in the (x, z) plane.
  static void map(const std::string &title,
                  const arma::mat &mat,
                  const Mesh &mesh,
                  const std::string &xlabel = "z [fm]",
                  const std::string &ylabel = "r [fm]");

  /// Plot a matrix.
  static void mat(const std::string &title,
                  const arma::mat &mat);

  /// Send a command to the ploting FIFO.
  static void sendFifo(std::string mesg = "");                        // #TEST#

  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;
};

#endif // PLOT_H
