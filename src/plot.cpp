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

#include <ctime>
#include <fcntl.h>

#include "tools.h"
#include "plot.h"
#include "mesh.h"

/** \file
 *  \brief Methods of the Plot class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::save(const std::string filename)
{
  Plot::sendFifo("@save \"" + filename + "\"");
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::clear(const std::string &name)
{
  if (name.empty())
    Plot::sendFifo("@clear");
  else
    Plot::sendFifo("@clear \"" + name + "\"");
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::slot(INT id)
{
  Plot::sendFifo(PF("@slot %d", id));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::sendFifo(std::string mesg)
{
  if (!useBokeh) return;

  INT fifo = open("/tmp/bokeh.fifo", O_WRONLY); //| O_NONBLOCK);

  if (fifo == -1) return;

  mesg += "\n";
retry:
  UINT nbBytes = write(fifo, mesg.c_str(), mesg.size());

  if (nbBytes != UINT(mesg.size()))
  {
    Tools::warning("incomplete write in Plot::sendFifo(): " + std::to_string(nbBytes));
    goto retry;
  }

  close(fifo);
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::curve(const std::string &xlabel, const std::string &ylabel, double x, double y)
{
  // sendFifo(PF("curve \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" %21.14e %21.14e", idStr.c_str(), title.c_str(), xlabel.c_str(), ylabel.c_str(), cname.c_str(), x, y));
  sendFifo(PF("@xy_data \"%s\" \"%s\" %21.14e %21.14e", xlabel.c_str(), ylabel.c_str(), x, y));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::curve(const std::string &xlabel, const std::string &ylabel, double y)
{
  double x = Tools::clock();

  // sendFifo(PF("curve \"%s\" \"%s\" \"%s\" \"%s\" \"%s\" %21.14e %21.14e", idStr.c_str(), title.c_str(), xlabel.c_str(), ylabel.c_str(), cname.c_str(), x, y));
  sendFifo(PF("@xy_data \"%s\" \"%s\" %21.14e %21.14e", xlabel.c_str(), ylabel.c_str(), x, y));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::curve(const std::string &xlabel, const std::string &ylabel, const arma::vec &x, const arma::vec &y)
{
  ASSERT(x.n_rows == y.n_rows, "incompatible dimensions: %d vs. %d.", x.n_rows, y.n_rows);

  for (UINT i = 0; i < x.n_rows; i++)
  {
    curve(xlabel, ylabel, x(i), y(i));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::mat(const std::string &title, const arma::mat &mat)
{
  for (UINT i = 0; i < mat.n_rows; i++)
  {
    std::string cmd = PF("@mat \"%s\" %d %d %d ", title.c_str(), mat.n_rows, mat.n_cols, i);

    for (UINT j = 0; j < mat.n_cols; j++)
    {
      cmd += PF("%21.14e ", mat(i, j));
    }

    sendFifo(cmd);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::map(const std::string &title, const arma::mat &mat, const Mesh &mesh, const std::string &xlabel, const std::string &ylabel)
{
  ASSERT(mesh.ax.p.n_rows == mat.n_rows, "incompatible dimensions: %d (x.rows) vs. %d (mat.rows).", mesh.ax.p.n_rows, mat.n_rows);
  ASSERT(mesh.az.p.n_rows == mat.n_cols, "incompatible dimensions: %d (z.rows) vs. %d (mat.cols).", mesh.az.p.n_rows, mat.n_cols);

  for (UINT i = 0; i < mat.n_rows; i++)
  {
    std::string cmd = PF("@map \"%s\" \"%s\" \"%s\" %f %f %f %f %d %d %d ", title.c_str(), xlabel.c_str(), ylabel.c_str(), mesh.az.p[0], mesh.ax.p[0], mesh.az.p[mesh.az.nb - 1], mesh.ax.p[mesh.ax.nb - 1], mat.n_cols, mat.n_rows, i);

    for (UINT j = 0; j < mat.n_cols; j++)
    {
      cmd += PF("%21.14e ", mat(i, j));
    }

    sendFifo(cmd);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Plot::field(const std::string &title, const arma::mat &mat, const Mesh &mesh, const std::string &xlabel, const std::string &ylabel)
{
  for (UINT i = 0; i < mat.n_cols; i++)
  {
    std::string cmd = PF("@map \"%s\" \"%s\" \"%s\" %f %f %f %f %d %d %d ",
                         title.c_str(),
                         xlabel.c_str(),
                         ylabel.c_str(),
                         mesh.ax.p[0],
                         mesh.ay.p[0],
                         mesh.ax.p[mesh.ax.nb - 1],
                         mesh.ay.p[mesh.ay.nb - 1],
                         mat.n_rows, mat.n_cols, i);

    for (UINT j = 0; j < mat.n_rows; j++)
    {
      cmd += PF("%21.14e ", mat(j, i));
    }

    sendFifo(cmd);
  }
}
