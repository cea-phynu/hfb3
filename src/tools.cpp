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

#include <chrono>
#include <ctime>
#include <ostream>
#include <sstream>
#include <fstream>
#include <fcntl.h>
#include <random>
#include <iomanip>

#include "global.h"
#include "tools.h"
#include "mesh.h"
#include "md5.h"
#include "zlib.h"
#include "gzstream.h"
#include "plot.h"
#include "multi.h"

/** \file
 *  \brief Methods of the Tools class.
 */

//==============================================================================
//==============================================================================
//==============================================================================

// Global quantities.

static INT timer_level = 0;
static arma::wall_clock wclock;
static std::string current_mesg;
static arma::vec timerStartTime;
static double progressStartTime;
static double now = 0.0;
static std::string progressCurrentMesg;
static INT progressCurrent = 0;
static std::mt19937 rngGenerator(1337);
static std::uniform_real_distribution<double> rngDistribution(0.0, 1.0);
static std::vector<StackItem> stack;

//==============================================================================
//==============================================================================
//==============================================================================


std::vector<std::string> utf8Chars = {std::string("\xE2\x94\x82\x00"),
                                      std::string("\xE2\x94\x80\x00"),
                                      std::string("\xE2\x94\x94\x00"),
                                      std::string("\xE2\x94\x8C\x00"),
                                      std::string("\xE2\x94\x90\x00"),
                                      std::string("\xE2\x94\x98\x00"),
                                      std::string("\xE2\x94\x9C\x00"),
                                      std::string("\xE2\x94\xB4\x00"),
                                      std::string("\xE2\x94\xA4\x00"),
                                      std::string("\xE2\x94\xAC\x00"),
                                      std::string("\xE2\x94\xBC\x00")
                                     };

std::vector<std::string> asciiChars = {std::string("|"),
                                       std::string("-"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+"),
                                       std::string("+")
                                      };

//==============================================================================
//==============================================================================
//==============================================================================

UINT hash2(std::string s)
{
  return hash2(s.c_str(), s.length());
}

//==============================================================================
//== Math
//==============================================================================

double Tools::accu(arma::vec in, const char* fileName, INT lineNumber)
{
  // 1. get the usual result
  double res = arma::accu(in);

  // shuffle the input
  arma::vec inshfl = arma::shuffle(in);

  // 2. accu the shuffle vector
  double resshfl = arma::accu(inshfl);

  // 3. precisely sum the shuffled vector
  __float128 bigres = 0;
  for (__float128 x: in) bigres += x;
  double preciseRes = (double) bigres;

  double err = fabs(res - preciseRes);
  if (err > 1e-10)
  {
    warning(PF(
      "[%s:%d] accu is not precise: accu %e, accu_shfl %e, precise %e (err: %e) [%e %e]x%d",
      fileName, lineNumber, res, resshfl, preciseRes, err, arma::min(arma::abs(in)), arma::max(in), in.n_elem
    ));
  }

  return preciseRes;
}

//==============================================================================
//==============================================================================
//==============================================================================

double Tools::getPfaffian(const arma::mat m)
{
  DBG_ENTER;

  arma::mat mat = m;
  UINT ndim = mat.n_rows;
  double t1 = 1.0;

  for (UINT j = 0; j < ndim / 2; j++)
  {
    t1 *= mat(0, 1);

    if (j > ndim / 2 - 1)
    {
      break;
    }

    UINT ndimr = mat.n_rows;

    for (UINT i = 2; i < ndimr; i++)
    {
      if (mat(0, 1) != 0.0)
      {
        arma::rowvec tv1 = mat.row(1) * mat(i, 0) / mat(1, 0);
        mat.row(i) -= tv1;
        arma::colvec tv2 = mat.col(1) * mat(0, i) / mat(0, 1);
        mat.col(i) -= tv2;
      }
      else
      {
        ERROR("need to pivot.");
      }
    }

    mat = mat(2, 2, arma::size(ndimr - 2, ndimr - 2));
  }

  DBG_RETURN(t1);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::humanSize(UINT bytes)
{
  UINT gb = 1024 * 1024 * 1024;
  UINT mb = 1024 * 1024;
  UINT kb = 1024;
  if( bytes >= gb) return std::to_string( (float)bytes/gb ) +  " GB ";
  if( bytes >= mb) return std::to_string( (float)bytes/mb ) +  " MB ";
  if( bytes >= kb) return std::to_string( (float)bytes/kb ) +  " KB ";
  return std::to_string(bytes) + " B ";
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::strIsospin(INT iso)
{
  std::string result = "unknown isospin !?";

  if (iso == NEUTRON) result = "neutron";

  if (iso ==  PROTON) result = "proton";

  if (iso ==   TOTAL) result = "total";

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::matToCrosses(const arma::mat &m, double epsilon)
{
  DBG_ENTER;

  std::stringstream ss;

  for (UINT i = 0; i < m.n_rows; i++)
  {
    for (UINT j = 0; j < m.n_cols; j++)
    {
      if (fabs(m(i, j)) > epsilon)
        ss << "x";
      else
        ss << ".";
    }

    ss << "\n";
  }

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::matToTxt(const arma::mat &m)
{
  DBG_ENTER;

  std::stringstream ss;

  for (UINT i = 0; i < m.n_rows; i++)
  {
    for (UINT j = 0; j < m.n_cols; j++)
    {
      ss << m(i, j) << " ";
    }

    ss << "\n";
  }

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::matToCsv(const arma::mat &m)
{
  DBG_ENTER;

  std::string result;

  for (UINT i = 0; i < m.n_rows; i++)
  {
    result += PF("%d,", i);

    for (UINT j = 0; j < m.n_cols; j++)
    {
      result += PF("%e", m(i, j));

      if (j != m.n_cols - 1) result += ",";
    }

    result += "\n";
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::cubeToDf3(const arma::cube &m)
{
  DBG_ENTER;

  std::stringstream ss(std::stringstream::out | std::stringstream::binary);
  INT nx = INT(m.n_rows);
  INT ny = INT(m.n_cols);
  INT nz = INT(m.n_slices);
  ss.put((unsigned char)(nx >> 8));
  ss.put((unsigned char)(nx & 0xff));
  ss.put((unsigned char)(ny >> 8));
  ss.put((unsigned char)(ny & 0xff));
  ss.put((unsigned char)(nz >> 8));
  ss.put((unsigned char)(nz & 0xff));
  double theMin = 0.0;
  double theMax = m.max();

  for (UINT k = 0; k < m.n_slices; k++)
  {
    for (UINT j = 0; j < m.n_cols; j++)
    {
      for (UINT i = 0; i < m.n_rows; i++)
      {
        UINT v = UINT(255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin));
        ss.put((unsigned char)v);
      }
    }
  }

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::matToVTK(const std::string &name, const arma::mat &m, const Mesh &mesh)
{
  DBG_ENTER;

  arma::cube tempCube(m.n_rows, m.n_cols, 1);
  tempCube.slice(0) = m;

  DBG_RETURN(cubeToVTK(name, tempCube, mesh));
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::cubeToVTK(const std::string &name, const arma::cube &m, const Mesh &mesh)
{
  DBG_ENTER;

  std::stringstream ss;
  ss << "<VTKFile type='ImageData' version='1.0' byte_order='LittleEndian' header_type='uint64'>" << "\n";
  ss << " <ImageData WholeExtent='";
  ss <<        0 << " " << mesh.ax.nb - 1;
  ss << " " << 0 << " " << mesh.ay.nb - 1;
  ss << " " << 0 << " " << mesh.az.nb - 1;
  ss << "' Origin='";
  ss <<        mesh.ax.p(0);
  ss << " " << mesh.ay.p(0);
  ss << " " << mesh.az.p(0);
  ss << "' Spacing='";
  double dx = 0.0;
  double dy = 0.0;
  double dz = 0.0;

  if (mesh.ax.nb != 1) dx = mesh.ax.p(1) - mesh.ax.p(0);

  if (mesh.ay.nb != 1) dy = mesh.ay.p(1) - mesh.ay.p(0);

  if (mesh.az.nb != 1) dz = mesh.az.p(1) - mesh.az.p(0);

  ss <<        dx;
  ss << " " << dy;
  ss << " " << dz;
  ss << "'>" << "\n";
  ss << "  <Piece Extent='";
  ss <<        0 << " " << mesh.ax.nb - 1;
  ss << " " << 0 << " " << mesh.ay.nb - 1;
  ss << " " << 0 << " " << mesh.az.nb - 1;
  ss << "'>" << "\n";
  ss << "   <PointData Scalars='density'>" << "\n";
  ss << "    <DataArray type='Float64' Name='" << name << "' format='appended' RangeMin='" << m.min() << "' RangeMax='" << m.max() << "' offset='0'>" << "\n";
  ss << "    </DataArray>" << "\n";
  ss << "   </PointData>" << "\n";
  ss << "   <CellData>" << "\n";
  ss << "   </CellData>" << "\n";
  ss << "  </Piece>" << "\n";
  ss << " </ImageData>" << "\n";
  UINT nb = m.n_cols * m.n_rows * m.n_slices;
#ifdef VTK_EXPORT_USE_BASE64
  ss << " <AppendedData encoding='base64'>" << "\n";
  ss << "_";
  INT lentoto = nb * sizeof(double) + 8;
  char *toto = new char[lentoto];
  memcpy(toto, (char *)( &nb ) + 4, 4);
  memcpy(toto + 4, (const char *)(&nb), 4);
  memcpy(toto + 8, (char *)(m.memptr()), nb * sizeof(double));
  ss << Base64::encode(toto, lentoto);
  delete[] toto;
#else
  ss << " <AppendedData encoding='raw'>" << "\n";
  ss << "_";
  ss.write(((char *)( &nb )) + 4, 4);
  ss.write(((char *)( &nb )), 4);
  ss.write((char *)( m.memptr() ), nb * sizeof(double));
#endif
  ss << "</AppendedData>" << "\n";
  ss << "</VTKFile>" << "\n";

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::cubeToRaw(const arma::cube &m)
{
  DBG_ENTER;

  std::stringstream ss(std::stringstream::out | std::stringstream::binary);
  double theMin = 0.0;
  double theMax = m.max();

  for (UINT k = 0; k < m.n_slices; k++)
  {
    for (UINT j = 0; j < m.n_cols; j++)
    {
      for (UINT i = 0; i < m.n_rows; i++)
      {
        UINT v = UINT(255 * (fabs(m(i, j, k)) - theMin) / (theMax - theMin));
        ss.put((unsigned char)v);
      }
    }
  }

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::trueColor(INT r, INT g, INT b)
{
  DBG_ENTER;
  DBG_RETURN(Tools::_printf("\033[38;2;%d;%d;%dm", r, g, b));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::progress(const std::string &mesg, INT p, INT total)
{
  DBG_ENTER;

  if (progressCurrentMesg != mesg)
  {
    progressCurrentMesg = mesg;
  }

  if (p < total)
  {
    UINT barLength = 40;
    UINT labelLength = 28;
    INT nb1 = INT((double)p * (double)barLength / (double)(total - 1));
    std::string mesgTemp;

    if (stringLen(mesg) > labelLength)
    {
      mesgTemp = mesg.substr(0, labelLength);
    }
    else
    {
      mesgTemp = mesg + std::string(labelLength - stringLen(mesg), ' ');
    }

    double length = clock() - progressStartTime;
    std::string colorTemp = color();

    if (length > 0.1) colorTemp = color("green");

    if (length > 1.0) colorTemp = color("yellow");

    if (length > 2.0) colorTemp = color("red");

    if (useColors)
    {
      // clear current line
      std::cout << "\33[2K\r";
      // draw progress bar
      std::cout << colorTemp << mesgTemp << color() << " [" << color("blue") << std::string(nb1, '=') << color() << std::string(barLength - nb1, ' ') << "] " << colorTemp << std::fixed << std::setprecision(2) << std::setw(5) << length << "s" << color() << std::flush;
    }

    if (p > total - 2)
    {
      // if no colors are used, only prINT the completed progress bar
      if (!useColors) std::cout << mesgTemp << " [" << std::string(nb1, '=') << std::string(barLength - nb1, ' ') << "] " << std::fixed << std::setprecision(2) << std::setw(5) << length << "s" << "\n";
      else std::cout << "\n";

      progressCurrent = -1;
      progressCurrentMesg = "";
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::progress(const std::string &mesg, INT total)
{
  DBG_ENTER;

  if (progressCurrentMesg != mesg)
  {
    progressCurrent = 0;
    progressCurrentMesg = mesg;
    progressStartTime = clock();
  }

  progress(progressCurrentMesg, progressCurrent, total);
  progressCurrent++;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::mat Tools::matFromCol(const arma::vec &v, UINT nb)
{
  DBG_ENTER;

  arma::mat result;

  if (nb == 0) nb = v.n_elem;

  result = arma::zeros(v.n_elem, nb);

  for (UINT i = 0; i < nb; i++)
    result.col(i) = v;

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::mat Tools::matFromRow(const arma::rowvec &v, UINT nb)
{
  DBG_ENTER;

  arma::mat result;

  if (nb == 0) nb = v.n_elem;

  result = arma::zeros(nb, v.n_elem);

  for (UINT i = 0; i < nb; i++)
    result.row(i) = v;

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::cube Tools::matTimesVec(const arma::mat &m, const arma::vec &v)
{
  DBG_ENTER;

  arma::cube result(m.n_rows, m.n_cols, v.n_rows);

  for (UINT i = 0; i < v.n_rows; i++)
    result.slice(i) = m * v(i);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::cube Tools::cubeTimesVec(const arma::cube &c, const arma::vec &v)
{
  DBG_ENTER;

  ASSERT(v.n_rows == c.n_slices, PF("wrong dimensions %d != %d", v.n_rows, c.n_slices));

  arma::cube result = arma::zeros(c.n_rows, c.n_cols, c.n_slices);

  for (UINT i = 0; i < v.n_rows; i++)
    result.slice(i) = c.slice(i) * v(i);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::cube Tools::cubeTimesMat(const arma::cube &c, const arma::mat &m)
{
  DBG_ENTER;

  ASSERT(c.n_rows == m.n_rows, PF("wrong dimensions %d != %d", c.n_rows, m.n_rows));
  ASSERT(c.n_cols == m.n_cols, PF("wrong dimensions %d != %d", c.n_cols, m.n_cols));

  arma::cube result = arma::zeros(c.n_rows, c.n_cols, c.n_slices);

  for (UINT i = 0; i < c.n_slices; i++)
    result.slice(i) = c.slice(i) % m;

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::mat Tools::cubeToMat(const arma::cube &c, INT iaxis, INT n)
{
  DBG_ENTER;

  arma::mat result;

  switch (iaxis)
  {
    case 0:
      result = arma::zeros(c.n_cols, c.n_slices);

      for (UINT i = 0; i < c.n_cols; i++)
        for (UINT j = 0; j < c.n_slices; j++)
          result(i, j) = c(n, i, j);

      break;

    case 1:
      result = arma::zeros(c.n_rows, c.n_slices);

      for (UINT i = 0; i < c.n_rows; i++)
        for (UINT j = 0; j < c.n_slices; j++)
          result(i, j) = c(i, n, j);

      break;

    case 2:
      result = c.slice(n);
      break;

    default:
      break;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg, const UMAT &m, bool dump)
{
  DBG_ENTER;

  info(mesg, arma::conv_to<IMAT >::from(m), dump);

  DBG_LEAVE;
}


//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg, const IMAT &m, bool dump)
{
  DBG_ENTER;

  if (msgToOut.count(MSG_INFO) == 0) DBG_LEAVE;

  info("[" + color("red") + mesg + color() + "] " + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols));

  if (dump)
  {
    std::stringstream ss;
    ss << TABLE_TD;
    ss.width(5);
    ss.setf(std::ios::scientific);

    for (UINT i = 0; i < m.n_cols; i++)
    {
      ss << TABLE_CENTER << TABLE_BLUE << i << TABLE_TD;
    }

    ss << TABLE_TR;

    for (UINT i = 0; i < m.n_rows; i++)
    {
      ss << TABLE_CENTER << TABLE_BLUE << i << TABLE_TD << TABLE_RIGHT << TABLE_NORM;

      for (UINT j = 0; j < m.n_cols; j++)
      {
        ss << m(i, j) << TABLE_TD;
      }

      ss << TABLE_TR;
    }

    info(printTable(ss.str()));
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg, const arma::mat &m, bool dump)
{
  DBG_ENTER;

  if (msgToOut.count(MSG_INFO) == 0) DBG_LEAVE;

  info(matToStr(mesg, m, dump));

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::matToStr(const std::string &mesg, const arma::mat &m, bool dump, const std::vector<std::string> lineLabels, const std::vector<std::string> colLabels)
{
  DBG_ENTER;

  std::string result;

  if (dump)
  {
    std::string tableStr;
    tableStr += TABLE_BLUE + TABLE_LEFT + mesg + TABLE_TD;

    for (UINT i = 0; i < m.n_cols; i++)
    {
      std::string label = colLabels.empty() ? PF("%d", i) : colLabels.at(i);
      tableStr += TABLE_CENTER + TABLE_YELLOW + label + TABLE_TD;
    }

    tableStr += TABLE_TR;

    for (UINT i = 0; i < m.n_rows; i++)
    {
      std::string label = lineLabels.empty() ? PF("%d", i) : lineLabels.at(i);
      tableStr += TABLE_LEFT + TABLE_GREEN + label + TABLE_TD + TABLE_RIGHT + TABLE_NORM;

      for (UINT j = 0; j < m.n_cols; j++)
      {
        tableStr += PF("%14.6e", m(i, j)) + TABLE_TD;
      }

      tableStr += TABLE_TR;
    }

    result += Tools::printTable(tableStr);
  }
  else
  {
    result += "[" + color("red") + mesg + color() + "] ";
    result += std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols) + " ";

    if (m.empty()) result += PF_RED("empty matrix");
    else
    {
      result += "MD5: " + color("green") + armaMD5(m) + " " + color();
      result += PF("range: [%e, %e] -- norm [%e]", m.min(), m.max(), arma::norm(m, "inf"));
    }

    result += "\n";
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg, const arma::cube &m, bool dump)
{
  DBG_ENTER;

  if (msgToOut.count(MSG_INFO) == 0) DBG_LEAVE;

  info("[" + color("red") + mesg + color() + "] " + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols) + "x" + std::to_string(m.n_slices) + " MD5: " + color("green") + armaMD5(m) + color() + PF(" range: [%e, %e]", m.min(), m.max()));

  if (dump)
  {
    for (INT i = 0; i < m.n_slices; i++)
    {
      Tools::info(PF("slice %d", i), m.slice(i), true);
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg, const ICUBE &m, bool dump)
{
  DBG_ENTER;

  if (msgToOut.count(MSG_INFO) == 0) DBG_LEAVE;

  info("[" + color("red") + mesg + color() + "] " + std::to_string(m.n_rows) + "x" + std::to_string(m.n_cols) + "x" + std::to_string(m.n_slices) + color() + PF(" range: [%e, %e]", m.min(), m.max()));

  if (dump)
  {
    for (INT i = 0; i < m.n_slices; i++)
    {
      Tools::info(PF("slice %d", i), m.slice(i), true);
    }
  }

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::matPrint(const arma::mat &m)
{
  DBG_ENTER;

  std::stringstream ss;

  for (UINT i = 0; i < m.n_rows; i++)
  {
    for (UINT j = 0; j < m.n_cols; j++)
    {
      if (fabs(m(i, j)) > 1e-8)
        ss << 'X';
      else
        ss << '.';
    }

    ss << "\n";
  }

  DBG_RETURN(ss.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::vec Tools::gaussian(const arma::vec &g, const arma::vec &x)
{
  DBG_ENTER;

  arma::vec result = arma::zeros(arma::size(x));
  result = arma::exp(-arma::square(x - g(1)) / g(2)) * g(0);

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

INT Tools::irand(INT imax)
{
  DBG_ENTER;

  DBG_RETURN((INT)(drand() * (double)imax));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::setSeed(INT seed)
{
  DBG_ENTER;

  rngGenerator.seed(seed);

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::growIVec(IVEC &m, INT v)
{
  DBG_ENTER;

  m.resize(m.size() + 1);
  m(m.size() - 1) = v;

  DBG_LEAVE;
}

//==============================================================================
//==============================================================================
//==============================================================================

const arma::mat Tools::randMat(UINT nbr, UINT nbc)
{
  DBG_ENTER;

  arma::mat result(nbr, nbc);

  for (UINT ir = 0; ir < nbr; ir++)
    for (UINT ic = 0; ic < nbc; ic++)
      result(ir, ic) = drand();

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::compressString(const std::string &str)
{
  DBG_ENTER;

  z_stream zs;
  memset(&zs, 0, sizeof(zs));

  if (deflateInit2(
        &zs,
        Z_BEST_COMPRESSION,
        Z_DEFLATED,
        MOD_GZIP_ZLIB_WINDOWSIZE + 16,
        MOD_GZIP_ZLIB_CFACTOR,
        Z_DEFAULT_STRATEGY
      ) != Z_OK)
  {
    ERROR("deflateInit2 failed while compressing.");
  }

  zs.next_in = (Bytef *)str.data();
  zs.avail_in = (unsigned int)(str.size());
  INT ret;
  char outbuffer[32768];
  std::string outstring;

  do
  {
    zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
    zs.avail_out = sizeof(outbuffer);
    ret = deflate(&zs, Z_FINISH);

    if (outstring.size() < zs.total_out)
    {
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  }
  while (ret == Z_OK);

  deflateEnd(&zs);

  if (ret != Z_STREAM_END) ERROR("zlib compression error. (%d)", ret);

  DBG_RETURN(outstring);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::readFile(const std::string &filename)
{
  DBG_ENTER;

  std::stringstream ss;
  std::ifstream stream(filename.c_str());

  if (stream.is_open())
  {
    while (stream.peek() != EOF)
    {
      ss << (char) stream.get();
    }

    stream.close();
  }
  else
  {
    ERROR("can not open file `" + filename + "`.");
  }

  std::string result = ss.str();

  // Decompress if needed
  try
  {
    std::string deflated = Tools::decompressString(result);
    result = deflated;
  }
  catch (...)
  {
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::readCin(void)
{
  DBG_ENTER;

  std::stringstream ss;
  std::string line;

  while (std::getline(std::cin, line))
  {
    ss << line << "\n";
  }

  std::string result = ss.str();

  // Decompress if needed
  try
  {
    std::string deflated = Tools::decompressString(result);
    result = deflated;
  }
  catch (...)
  {
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::armaMD5(const arma::vec &vec)
{
  DBG_ENTER;

  MD5 md5;
  md5.update((unsigned char *)(vec.memptr()), (unsigned int)(vec.n_elem * sizeof(double)));
  md5.finalize();

  DBG_RETURN(md5.hexdigest());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::armaMD5(const arma::mat &mat)
{
  DBG_ENTER;

  MD5 md5;
  md5.update((unsigned char *)(mat.memptr()), (unsigned int)(mat.n_elem * sizeof(double)));
  md5.finalize();

  DBG_RETURN(md5.hexdigest());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::armaMD5(const arma::cube &cube)
{
  DBG_ENTER;

  MD5 md5;
  md5.update((unsigned char *)(cube.memptr()), (unsigned int)(cube.n_elem * sizeof(double)));
  md5.finalize();

  DBG_RETURN(md5.hexdigest());
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::end(INT code)
{
  while (!stack.empty()) Tools::timerEnd();

  exit(code);
}

//==============================================================================
//==============================================================================
//==============================================================================

UINT Tools::stringLen(const std::string &s)
{
  INT len = 0;
  auto ptr = s.begin();

  while (*ptr) len += (*ptr++ & 0xC0) != 0x80;

  return len;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string operator*(const std::string &s, const INT &n)
{
  std::string result;

  if (n < 1)
  {
    return result;
  };

  for (UINT i = 0; i < n; i++)
    result += s;

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::boxed(std::string input)
{
  return Tools::printTable(input + TABLE_TD + TABLE_TR);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::printTable(std::string input, INT style)
{
  DBG_ENTER;

  if (style == -1)
  {
    style = tableStyle;
  }

  bool lines  = true;
  bool border = true;
  bool paddingLeft  = true;
  bool paddingRight = true;

  if (style == 1)
  {
    lines  = false;
    border = true;
    paddingLeft  = true;
    paddingRight = true;
  }

  if (style == 2)
  {
    lines = false;
    border = false;
    paddingLeft  = false;
    paddingRight = true;
  }

  std::stringstream result;
  std::vector<std::vector<std::string> > table;
  std::vector<std::vector<std::string> > align;
  std::vector<std::vector<std::string> > color;
  std::size_t icurrent = 0;
  bool last = false;
  UINT i = 0;
  UINT j = 0;
  UINT maxcols = 0;
  std::vector<std::string> inputs;

  while (!last)
  {
    if (input.substr(icurrent, 2) == "##")
    {
      std::string tag = input.substr(icurrent, 4);
      icurrent += 4;
      inputs.push_back(tag);
    }
    else
    {
      std::size_t ipos = input.find("##", icurrent);

      if (ipos == std::string::npos)
      {
        last = true;
        ipos = input.size();
      }

      if (ipos != icurrent)
      {
        std::string text = input.substr(icurrent, ipos - icurrent);
        icurrent += ipos - icurrent;
        inputs.push_back(text);
      }
    }
  }

  std::string currentAlign = TABLE_RIGHT; // Default horizontal alignment.
  std::string currentColor = ""; // Default color.

  for (UINT k = 0; k < inputs.size(); k++)
  {
    if (i + 1 > table.size())
    {
      table.push_back(std::vector<std::string>());
      align.push_back(std::vector<std::string>());
      color.push_back(std::vector<std::string>());
    }

    if (j + 1 > table[i].size())
    {
      table[i].push_back("");
      align[i].push_back("");
      color[i].push_back("");
      maxcols = MAX(maxcols, table[i].size());
    }

    if (inputs[k] == TABLE_TD)
    {
      j++;
    }
    else if (inputs[k] == TABLE_TR)
    {
      i++;
      j = 0;
    }
    else if (inputs[k] == TABLE_CENTER)
    {
      currentAlign = TABLE_CENTER;
    }
    else if (inputs[k] == TABLE_LEFT)
    {
      currentAlign = TABLE_LEFT;
    }
    else if (inputs[k] == TABLE_RIGHT)
    {
      currentAlign = TABLE_RIGHT;
    }
    else if (inputs[k] == TABLE_NORM)
    {
      currentColor = "";
    }
    else if (inputs[k] == TABLE_RED)
    {
      currentColor = "red";
    }
    else if (inputs[k] == TABLE_GREEN)
    {
      currentColor = "green";
    }
    else if (inputs[k] == TABLE_BLUE)
    {
      currentColor = "blue";
    }
    else if (inputs[k] == TABLE_YELLOW)
    {
      currentColor = "yellow";
    }
    else
    {
      table[i][j] = inputs[k];
      align[i][j] = currentAlign;
      color[i][j] = currentColor;
    }
  }

  table.push_back(std::vector<std::string>());
  align.push_back(std::vector<std::string>());
  color.push_back(std::vector<std::string>());
  /* maxcols++; */

  // All rows must have the same number of columns.
  for (UINT i = 0; i < table.size(); i++)
  {
    while (table[i].size() < maxcols)
    {
      table[i].push_back("");
      align[i].push_back("");
      color[i].push_back("");
    }
  }

  UVEC emptyLine = arma::ones<UVEC>(table.size() * 2 + 1);

  // Calculate the width of the columns.
  IVEC width = arma::zeros<IVEC>(maxcols);
  UMAT fill = arma::zeros<UMAT>(table.size() + 2, maxcols + 2);

  for (UINT j = 0; j < maxcols; j++)
  {
    for (UINT i = 0; i < table.size(); i++)
    {
      width(j) = MAX(width(j), INT(stringLen(table[i][j])));

      if (table[i][j] != "")
      {
        fill(i + 1, j + 1) = 1;
        emptyLine(i * 2 + 1) = 0;
      }
    }
  }

  // Apply the alignment tags.
  for (UINT i = 0; i < table.size(); i++)
  {
    for (UINT j = 0; j < maxcols; j++)
    {
      INT left = 0;
      INT right = 0;
      INT w = width(j);
      INT s = INT(stringLen(table[i][j]));

      if (align[i][j] == TABLE_RIGHT)
      {
        left = w - s;
        right = w - s - left;
      }

      if (align[i][j] == TABLE_LEFT)
      {
        right = w - s;
        left = w - s - right;
      }

      if (align[i][j] == TABLE_CENTER)
      {
        right = (w - s) / 2;
        left = w - s - right;
      }

      table[i][j] = (paddingLeft ? " " : "")
                  + std::string(left, ' ')
                  + table[i][j]
                  + std::string(right, ' ')
                  + (paddingRight ? " " : "");
    }
  }

  for (UINT j = 0; j < maxcols; j++)
  {
    if (paddingLeft)  width(j) += 1;

    if (paddingRight) width(j) += 1;
  }

  // Create and fill the final table.
  Multi<std::string> finalTable;


  for (INT i = 0; i < table.size(); i++) 
  {
    for (INT j = 0; j < maxcols; j++) 
    {
      if (table[i][j] != "")
      {
        finalTable(i * 2 + 1, j * 2 + 1) = Tools::color(color[i][j])
                                         + std::string(" ") * (width(j) - INT(stringLen(table[i][j])))
                                         + table[i][j] + Tools::color();
      }
      else
      {
        finalTable(i * 2 + 1, j * 2 + 1) = "";
      }
    }
  }

  if (lines)
  {
    for (INT i = 0; i < table.size(); i++) 
    {
      for (INT j = 0; j < maxcols; j++) 
      {
        if (fill(i + 1, j + 1) == 1)
        {
          if (fill(i, j + 1) == 1) finalTable(i * 2 + 0, j * 2 + 1) = getChar(CHAR_RL) * width(j);

          if (fill(i + 2, j + 1) == 1) finalTable(i * 2 + 2, j * 2 + 1) = getChar(CHAR_RL) * width(j);

          if (fill(i + 1, j    ) == 1) finalTable(i * 2 + 1, j * 2 + 0) = getChar(CHAR_TB);

          if (fill(i + 1, j + 2) == 1) finalTable(i * 2 + 1, j * 2 + 2) = getChar(CHAR_TB);

          finalTable(i * 2 + 0, j * 2 + 0) = "x";
          finalTable(i * 2 + 2, j * 2 + 0) = "x";
          finalTable(i * 2 + 0, j * 2 + 2) = "x";
          finalTable(i * 2 + 2, j * 2 + 2) = "x";
        }
      }
    }
  }

  if (border)
  {
    for (INT i = 0; i < table.size(); i++) 
    {
      for (INT j = 0; j < maxcols; j++) 
      {
        if (fill(i + 1, j + 1) == 1)
        {
          if (fill(i, j + 1) == 0) finalTable(i * 2 + 0, j * 2 + 1) = getChar(CHAR_RL) * width(j);

          if (fill(i + 2, j + 1) == 0) finalTable(i * 2 + 2, j * 2 + 1) = getChar(CHAR_RL) * width(j);

          if (fill(i + 1, j    ) == 0) finalTable(i * 2 + 1, j * 2 + 0) = getChar(CHAR_TB);

          if (fill(i + 1, j + 2) == 0) finalTable(i * 2 + 1, j * 2 + 2) = getChar(CHAR_TB);

          finalTable(i * 2 + 0, j * 2 + 0) = "x";
          finalTable(i * 2 + 2, j * 2 + 0) = "x";
          finalTable(i * 2 + 0, j * 2 + 2) = "x";
          finalTable(i * 2 + 2, j * 2 + 2) = "x";
        }
      }
    }
  }

  for (INT i = 0; i < table.size(); i++) 
  {
    for (INT j = 0; j < maxcols; j++) 
    {
      if (finalTable(i * 2 + 0, j * 2 + 0) != "x") continue;

      bool it = false;
      bool ir = false;
      bool ib = false;
      bool il = false;

      if ((i > 0               ) && (finalTable(i * 2 - 1, j * 2) != "")) it = true;

      if ((j < maxcols - 1     ) && (finalTable(i * 2, j * 2 + 1) != "")) ir = true;

      if ((i < table.size() - 1) && (finalTable(i * 2 + 1, j * 2) != "")) ib = true;

      if ((j > 0               ) && (finalTable(i * 2, j * 2 - 1) != "")) il = true;

      std::string toto = "";

      if ( it &&   ir && !ib && !il)
      {
        toto = getChar(CHAR_TR);
        emptyLine(i * 2 + 0) = 0;
      }

      if ( it &&  !ir &&  ib && !il)
      {
        toto = getChar(CHAR_TB);
      }

      if ( it &&  !ir && !ib &&  il)
      {
        toto = getChar(CHAR_TL);
        emptyLine(i * 2 + 0) = 0;
      }

      if (!it &&   ir &&  ib && !il)
      {
        toto = getChar(CHAR_RB);
        emptyLine(i * 2 + 0) = 0;
      }

      if (!it &&   ir && !ib &&  il)
      {
        toto = getChar(CHAR_RL);
        emptyLine(i * 2 + 0) = 0;
      }

      if (!it &&  !ir &&  ib &&  il)
      {
        toto = getChar(CHAR_BL);
        emptyLine(i * 2 + 0) = 0;
      }

      if ( it &&   ir &&  ib && !il)
      {
        toto = getChar(CHAR_TRB);
        emptyLine(i * 2 + 0) = 0;
      }

      if ( it &&   ir && !ib &&  il)
      {
        toto = getChar(CHAR_TRL);
        emptyLine(i * 2 + 0) = 0;
      }

      if ( it &&  !ir &&  ib &&  il)
      {
        toto = getChar(CHAR_TBL);
        emptyLine(i * 2 + 0) = 0;
      }

      if (!it &&   ir &&  ib &&  il)
      {
        toto = getChar(CHAR_RBL);
        emptyLine(i * 2 + 0) = 0;
      }

      if ( it &&   ir &&  ib &&  il)
      {
        toto = getChar(CHAR_TRBL);
        emptyLine(i * 2 + 0) = 0;
      }

      finalTable(i * 2 + 0, j * 2 + 0) = toto;
    }
  }

  // Re-calculate the width of the columns.
  UVEC finalWidth = arma::zeros<UVEC>(maxcols * 2 + 1);

  for (UINT j = 0; j < maxcols; j++)
  {
    finalWidth(j * 2 + 1) = width(j);
  }

  for (INT j = 0; j < maxcols * 2 + 1; j += 2) 
  {
    finalWidth(j) = 0;

    for (INT i = 0; i < table.size() * 2 + 1; i++) 
    {
      finalWidth(j) = MAX(finalWidth(j), stringLen(finalTable(i, j)));
    }
  }

  // Fill with spaces if needed.
  for (INT i = 0; i < table.size() * 2 + 1; i++) 
  {
    for (INT j = 0; j < maxcols * 2 + 1; j++) 
    {
      while (stringLen(finalTable(i, j)) < finalWidth(j))
        finalTable(i, j) += " ";
    }
  }

  for (INT i = 0; i < table.size() * 2 - 1; i++) 
  {
    std::string toto = "";

    for (INT j = 0; j < maxcols * 2 - 1; j++)
    {
      toto += finalTable(i, j);
    }

    toto = trim_e(toto);

    if (emptyLine(i) == 0)
    {
      result << toto;

      if (i != table.size() * 2 - 2) result << "\n";
    }
  }

  DBG_RETURN(result.str());
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::valueTable(const std::string &title,
                              const std::list<std::string> &columnLabel,
                              const std::list<std::string> &unit,
                              const std::list<std::list<std::string> > &values)
{
  DBG_ENTER;

  std::string result;

  result += TABLE_LEFT   + TABLE_BLUE + title + TABLE_TD;
  result += TABLE_CENTER + TABLE_YELLOW;

  for (const auto &c : columnLabel) result += c + TABLE_TD;

  bool useUnits = false;

  for (const auto &c : unit) if (!c.empty())
    {
      useUnits = true;
      break;
    }

  if (useUnits)
  {
    result += TABLE_TR;
    result += TABLE_LEFT   + TABLE_YELLOW + " " + TABLE_TD;
    result += TABLE_CENTER + TABLE_YELLOW;

    for (const auto &c : unit) result += c + TABLE_TD;
  }

  for (const auto &p : values)
  {
    result += TABLE_NORM + TABLE_TR;

    bool first = true;

    for (const auto &v : p)
    {
      if (first)
      {
        first = false;
        result += TABLE_LEFT + TABLE_GREEN + v + TABLE_RIGHT + TABLE_NORM + TABLE_TD;
      }
      else
      {
        result += v + TABLE_TD;
      }

    }
  }

  DBG_RETURN(printTable(result));
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::printOut(const std::string &mesg)
{
  std::cout.setf(std::ios::unitbuf);
  std::string preMesg = color();

  if (msgToOut.count(MSG_PID) != 0)
  {
    preMesg += Tools::_printf("[%06ld] ", getpid());
  }

  if (msgToOut.count(MSG_TIME) != 0)
  {
    preMesg += Tools::_printf("%08.3f ", clock());
  }

#define USE_LOGGER
#ifdef  USE_LOGGER
  logger.log(preMesg + mesg);
#else
  std::cout << preMesg << mesg << color() << std::endl;
#endif
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::debug(const std::string &mesg)
{
  if (msgToOut.count(MSG_DEBUG) == 0) return;

  std::vector<std::string> lines = stringSplit(mesg, '\n');

  for (UINT l = 0; l < lines.size(); l++)
  {
    if (!lines[l].empty()) printOut(color("white") + "#DEBUG# " + lines[l]);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::warning(const std::string &mesg)
{
  if (msgToOut.count(MSG_WARNING) == 0) return;

  std::vector<std::string> lines = stringSplit(mesg, '\n');

  for (UINT l = 0; l < lines.size(); l++)
  {
    if (!lines[l].empty()) printOut(color("orange") + "WARNING " + lines[l] + color());
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::info(const std::string &mesg)
{
  if (msgToOut.count(MSG_INFO) == 0) return;

  std::vector<std::string> lines = stringSplit(mesg, '\n');

  for (UINT l = 0; l < lines.size(); l++)
  {
    if (!lines[l].empty()) printOut(lines[l]);
    //    Plot::log("log", lines[l]);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::mesg(const std::string &cat, const std::string &mesg)
{
  if (msgToOut.count(MSG_MAIN) == 0) return;

  std::vector<std::string> lines = stringSplit(mesg, '\n');

  for (UINT l = 0; l < lines.size(); l++)
  {
    if (!lines[l].empty()) printOut(Tools::color("magenta") + cat + ": " + Tools::color()  + Tools::trim_e(lines[l]));

    //    Plot::log("log", lines[l]);
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::getStack(void)
{
  std::string result;

  INT indentLevel = 0;

  for (auto &s : stack)
  {
    std::string name = PF_YELLOW(s.name);
    std::string indent;

    for (INT i = 0; i < indentLevel - 1; i++)
    {
      indent += "  ";
    }

    if (indentLevel == 0) result += name + "\n";
    else                  result += indent + getChar(CHAR_TR) + " " + name + " " + s.location + "\n";

    indentLevel++;
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::_error(const std::string &mesg)
{
  std::string finalMesg = mesg + "\n";
  finalMesg += getStack();

  if (exitOnError)
  {
    std::vector<std::string> lines = stringSplit(finalMesg, '\n');

    for (UINT l = 0; l < lines.size(); l++)
    {
      if (!lines[l].empty()) printOut(color("red") + "ERROR!: " + color() + lines[l]);
    }

    exit(-1);
  }

  throw std::string(color("red") + "ERROR!? " + color() + finalMesg);
}

//==============================================================================
//==============================================================================
//==============================================================================

INT Tools::ipow(INT base, INT exp)
{
  DBG_ENTER;

  INT result = 1;

  while (exp)
  {
    if (exp & 1)
      result *= base;

    exp >>= 1;
    base *= base;
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::vector<double> Tools::vecToStd(const arma::vec &v)
{
  DBG_ENTER;

  DBG_RETURN(arma::conv_to<std::vector<double> >::from(v));
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::cube Tools::vecToCube(arma::vec X, UINT nRows, UINT nCols)
{
  DBG_ENTER;

  arma::cube c = arma::cube(nRows, nCols, X.n_elem);

  for (UINT i = 0; i < nRows; i++)
  {
    for (UINT j = 0; j < nCols; j++)
    {
      for (UINT k = 0; k < X.n_elem; k++)
        c(i, j, k) = X(k);
    }
  }

  DBG_RETURN(c);
}

//==============================================================================
//==============================================================================
//==============================================================================

arma::cube Tools::matToCube(arma::mat m, UINT nSlices)
{
  DBG_ENTER;

  arma::cube c = arma::cube(m.n_rows, m.n_cols, nSlices);
  UINT i;

  for (i = 0; i < nSlices; i++)
  {
    c.slice(i) = m;
  }

  DBG_RETURN(c);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::date(void)
{
  DBG_ENTER;

  auto nowTime = std::chrono::system_clock::now();
  std::time_t now_time = std::chrono::system_clock::to_time_t(nowTime);

  DBG_RETURN(std::ctime(&now_time));
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::save(const std::string &result, const std::string &filename, bool compress)
{
  DBG_ENTER;

  if (compress)
  {
    ogzstream fp(filename.c_str());
    fp << result;
  }
  else
  {
    std::filebuf fb;
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << result;
    fb.close();
  }

  DBG_RETURN(filename);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::vecToStr(const arma::vec &vec)
{
  DBG_ENTER;

  std::string result = "[";

  for (UINT i = 0; i < vec.n_rows; i++)
  {
    if (i != 0) result += ", ";

    result += Tools::_printf("%7.4e", vec(i));
  }

  result += "]";

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::ivecToStr(const IVEC &vec)
{
  DBG_ENTER;

  std::string result = "[";

  for (UINT i = 0; i < vec.n_rows; i++)
  {
    if (i != 0) result += ", ";

    result += Tools::_printf("%d", vec(i));
  }

  result += "]";

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

bool Tools::checkSymmetry(const arma::mat &mat, const std::string &mesg)
{
  DBG_ENTER;

  bool result = mat.is_symmetric(1e-12);

  if (!result)
  {
    double vError = arma::abs(mat - mat.t()).max();

    debug(PF("Matrix is not 100%% symmetric (%e): ", vError) + mesg);
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::version(void)
{
  DBG_ENTER;

  std::string result = "";

  result += Tools::color("green") + "HFB3 (" + std::string(CFG_GIT_VERSION) + ")" + Tools::color();

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::treeStr(const std::vector<std::pair<std::string, std::string> > values, bool isShort)
{
  DBG_ENTER;

  std::string result;

  UINT nb = values.size();

  if (isShort)
  {
    result += "(";

    for (UINT i = 0; i < nb; i++)
    {
      result += PF_YELLOW(trim(values.at(i).first)) + ":" + PF_BLUE(trim(values.at(i).second));

      if (i != nb - 1) result += ",";
    }

    result += ")";
  }
  else
  {
    if (nb == 0)
    {
      DBG_RETURN("");
    }

    if (nb == 1)
    {
      result += PF_YELLOW(values.at(0).first) + ": " + values.at(0).second;
      DBG_RETURN(result);
    }
    else
    {
      result += PF_GREEN(values.at(0).first) + " " + values.at(0).second;

      for (INT i = 1; i < nb; i++)
      {
        std::string barStr = getChar(CHAR_TB);

        if (i == nb - 1) barStr = getChar(CHAR_TR);

        std::vector<std::string> lines = stringSplit(values.at(i).second, '\n');

        if (lines.size() == 0)
        {
          result += "\n" + PF_GREEN(barStr) + " " + PF_YELLOW(values.at(i).first) + ": \"\"";
        }

        if (lines.size() == 1)
        {
          result += "\n" + PF_GREEN(barStr) + " " + PF_YELLOW(values.at(i).first) + ": " + values.at(i).second;
        }


        if (lines.size() > 1)
        {
          for (UINT l = 0; l < lines.size(); l++)
          {
            std::string barStr2 = getChar(CHAR_TB);

            if ((i == nb - 1) && (l == lines.size() - 1)) barStr2 = getChar(CHAR_TR);

            result += "\n" + PF_GREEN(barStr2) + " " + PF_YELLOW(lines[l]);
          }
        }
      }
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const States &, bool)
{
  DBG_ENTER;

  std::string result = "";

  // TODO

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const bool &value)
{
  DBG_ENTER;

  std::string result = "";

  if (value) result = PF_GREEN("True");
  else       result = PF_RED("False");

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const std::vector<std::string> &list)
{
  DBG_ENTER;

  std::string result = "[";

  for (UINT i = 0; i < list.size(); i++)
  {
    result += "\"" + PF_BLUE(list.at(i)) + "\"";

    if (i != list.size() - 1) result += ", ";
  }

  result += "]";

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

bool Tools::eig_sym(arma::vec &eigVals, arma::mat &eigVecs, const arma::mat &mat)
{
  DBG_ENTER;

  bool result = arma::eig_sym(eigVals, eigVecs, arma::symmatu(mat));

  for (UINT i = 0; i < eigVals.n_elem; i++)
  {
    if (arma::index_min(eigVecs.col(i)) > arma::index_max(eigVecs.col(i)))
    {
      eigVecs.col(i) *= -1;
    }
  }

  DBG_RETURN(result);
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const INT &v, bool)
{
  DBG_ENTER;

  DBG_RETURN(PF("%d", v));
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const UINT &v, bool)
{
  DBG_ENTER;

  DBG_RETURN(PF("%d", v));
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const double &v, bool)
{
  DBG_ENTER;

  DBG_RETURN(PF("%.6e", v));
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const std::string &v, bool)
{
  DBG_ENTER;

  DBG_RETURN(PF(v));
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const arma::vec   &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%d", v.n_elem));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%d (%d bytes)", v.n_elem, v.n_elem * sizeof(double)) + ", " + PF("%9.3e", arma::norm(v, "inf")) + ", " + PF_GREEN(armaMD5(v)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const arma::mat   &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%d", v.n_rows, v.n_cols));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_rows * v.n_cols * sizeof(double)) + ", " + PF("%9.3e", arma::norm(v, "inf")) + ", " + PF_GREEN(armaMD5(v)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const arma::cube  &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d", v.n_rows, v.n_cols, v.n_slices));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_slices, v.n_rows * v.n_cols * v.n_slices * sizeof(double)) + ", " + PF_GREEN(armaMD5(v)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const IVEC  &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%d", v.n_elem));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%d (%d bytes)", v.n_elem, v.n_elem * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const UVEC  &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%d", v.n_elem));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%d (%d bytes)", v.n_elem, v.n_elem * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const IMAT  &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%d", v.n_rows, v.n_cols));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_rows * v.n_cols * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const UMAT  &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%d", v.n_rows, v.n_cols));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_rows * v.n_cols * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const ICUBE &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d", v.n_rows, v.n_cols, v.n_slices));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_slices, v.n_rows * v.n_cols * v.n_slices * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::infoStr(const UCUBE &v, bool isShort)
{
  DBG_ENTER;

  if (isShort)
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d", v.n_rows, v.n_cols, v.n_slices));
  }
  else
  {
    DBG_RETURN(PF_BLUE("%dx%dx%d (%d bytes)", v.n_rows, v.n_cols, v.n_slices, v.n_rows * v.n_cols * v.n_slices * sizeof(INT)));
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

INT Tools::findColInMat(const arma::vec &v, const arma::mat &m)
{
  DBG_ENTER;

  bool found = false;

  INT j = 0;

  for (j = 0; j < m.n_cols; j++)
  {
    if (arma::norm(v - m.col(j), "inf") < 1e-12)
    {
      found = true;
      break;
    }
  }

  DBG_RETURN(found ? j : -1);
}

//==============================================================================
//==============================================================================
//==============================================================================

// No DBG_ENTER / DBG_LEAVE / DBG_RETURN in the following methods

//==============================================================================
//==============================================================================
//==============================================================================

double Tools::clock(void)
{
  now += wclock.toc();
  wclock.tic();
  return now;
}

//==============================================================================
//==============================================================================
//==============================================================================

void Tools::timer(const std::string &label, const std::string &location, const std::string &_color)
{
  // std::cout << "timer    : " << label << std::endl;

  if (msgToOut.count(MSG_TIMELINE) == 0)
  {
    stack.push_back({label, location, ""});
    return;
  }

  std::string colorStr = (_color == "") ? randomColor(timer_level) : _color;

  timer_level += 1;
  stack.push_back({label, location, colorStr});

  timerStartTime.resize(timer_level);
  timerStartTime(timer_level - 1) = clock();
  std::string indent;

  for (INT i = 0; i < timer_level - 1; i++)
  {
    indent += color(stack[i].color);
    indent += getChar(CHAR_TB);
    indent += color();
    //    indent += " ";
  }

  Tools::printOut(indent + color(colorStr) + getChar(CHAR_RB) + color() + " " + label);
}

//==============================================================================
//==============================================================================
//==============================================================================

bool Tools::timerEnd(const std::string &label, const std::string &msg)
{
  // std::cout << "timer_end: " << label << std::endl;

  if (msgToOut.count(MSG_TIMELINE) == 0)
  {
    if (!stack.empty()) stack.pop_back();
    return false;
  }

  std::string &lastName = stack.back().name;

  if (label == "") ASSERT(false, "empty label in Tools::timerEnd()");

  if (label != lastName)
  {
    std::cout << "'" << label << "' != '" << lastName << "'" << std::endl;

    ASSERT(false, "missing matching timerEnd() in " + lastName);
    exit(-1);
    Tools::timerEnd(lastName);
  }

  if (timer_level == 0) return false;

  std::string indent;

  for (INT i = 0; i < timer_level - 1; i++)
  {
    indent += color(stack[i].color);
    indent += getChar(CHAR_TB);
    indent += color();
    //    indent += " ";
  }

  double length = clock() - timerStartTime(timer_level - 1);
  std::string colorTemp = color();

  if (length > 0.1) colorTemp = color("green");

  if (length > 1.0) colorTemp = color("yellow");

  if (length > 2.0) colorTemp = color("red");

  std::string mesg = "";

  if (msg != "") mesg = " " + msg;

  Tools::printOut(colorTemp + indent + color(stack[timer_level - 1].color) + getChar(CHAR_TR) + " " + colorTemp + stack.back().name + color() + " [" + colorTemp + Tools::_printf("%.4f", length) + " s" + color() + "]" + mesg);

  stack.pop_back();
  timer_level -= 1;
  return true;
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::randomColor(INT fixedColor)
{
  INT i = (fixedColor == -1) ? irand(7) : (fixedColor % 7);

  if (i == 0) return "red"    ;

  if (i == 1) return "green"  ;

  if (i == 2) return "yellow" ;

  if (i == 3) return "blue"   ;

  if (i == 4) return "magenta";

  if (i == 5) return "cyan"   ;

  return "";
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::color(const std::string &c)
{
  if (useColors)
  {

#ifdef ANSI_TRUECOLORS

    if (darkMode)
    {
      if (c == "red"    ) DBG_RETURN(trueColor(255, 200, 200));

      if (c == "green"  ) DBG_RETURN(trueColor(200, 255, 200));

      if (c == "yellow" ) DBG_RETURN(trueColor(200, 255, 255));

      if (c == "blue"   ) DBG_RETURN(trueColor(200, 200, 255));

      if (c == "magenta") DBG_RETURN(trueColor(255, 200, 255));

      if (c == "cyan"   ) DBG_RETURN(trueColor(200, 255, 255));

      if (c == "orange" ) DBG_RETURN(trueColor(255, 200,   0));
    }
    else
    {
      if (c == "red"    ) DBG_RETURN(trueColor(200,   0,   0));

      if (c == "green"  ) DBG_RETURN(trueColor(  0, 200,   0));

      if (c == "yellow" ) DBG_RETURN(trueColor(  0, 200, 200));

      if (c == "blue"   ) DBG_RETURN(trueColor(  0,   0, 200));

      if (c == "magenta") DBG_RETURN(trueColor(200,   0, 200));

      if (c == "cyan"   ) DBG_RETURN(trueColor(  0, 200, 200));

      if (c == "orange" ) DBG_RETURN(trueColor(200, 100,   0));
    }

#endif

#ifdef ANSI_256_COLORS

    if (c == "red"    ) return "\033[38;5;9m";   // "256 colors" red

    if (c == "green"  ) return "\033[38;5;10m";  // "256 colors" green

    if (c == "yellow" ) return "\033[38;5;11m";  // "256 colors" yellow

    if (c == "blue"   ) return "\033[38;5;12m";  // "256 colors" blue

    if (c == "magenta") return "\033[38;5;13m";  // "256 colors" magenta

    if (c == "cyan"   ) return "\033[38;5;14m";  // "256 colors" cyan

    if (c == "orange" ) return "\033[38;5;208m"; // "256 colors" orange

#endif

    if (c == "red"    ) return "\033[1;31m"; // default red

    if (c == "green"  ) return "\033[1;32m"; // default green

    if (c == "yellow" ) return "\033[1;33m"; // default yellow

    if (c == "blue"   ) return "\033[1;34m"; // default blue

    if (c == "magenta") return "\033[1;35m"; // default magenta

    if (c == "cyan"   ) return "\033[1;36m"; // default cyan

    if (c == "orange" ) return "\033[1;37m"; // default orange

    return "\033[0m";                          // default fg
  }

  return "";
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::getChar(INT id)
{
  ASSERT((id >= 0) && (id < CHAR_NB), "wrong char id.");

  if (useUtf8)
  {
    return utf8Chars[id];
  }
  else
  {
    return asciiChars[id];
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::trim(const std::string str)
{
  return trim_b(trim_e(str));
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::trim_b(const std::string str)
{
  if (str.size() == 0)
  {
    return str;
  }

  INT cfirst = INT(str.find_first_not_of(" "));
  std::string result = str.substr(cfirst);

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::vector<std::string> Tools::stringSplit(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;

  while (std::getline(ss, item, delim))
    elems.push_back(item);

  return elems;
}

//==============================================================================
//==============================================================================
//==============================================================================

double Tools::drand(double dmax)
{
  return rngDistribution(rngGenerator) * dmax;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::decompressString(const std::string &str)
{
  std::string outstring;
  z_stream zs;
  memset(&zs, 0, sizeof(zs));

  if (inflateInit2(&zs, MOD_GZIP_ZLIB_WINDOWSIZE + 16) != Z_OK)
  {
    ERROR("inflateInit2 failed while decompressing.");
  }

  zs.next_in = (Bytef *)str.data();
  zs.avail_in = (unsigned int)(str.size());
  INT ret;
  char outbuffer[32768];

  do
  {
    zs.next_out = reinterpret_cast<Bytef *>(outbuffer);
    zs.avail_out = sizeof(outbuffer);
    ret = inflate(&zs, 0);

    if (outstring.size() < zs.total_out)
    {
      outstring.append(outbuffer, zs.total_out - outstring.size());
    }
  }
  while (ret == Z_OK);

  inflateEnd(&zs);

  if (ret != Z_STREAM_END) throw(PF("zlib decompression error. (%d)", ret));

  return outstring;
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::trim_e(const std::string str)
{
  if (str.size() == 0)
  {
    return str;
  }

  INT clast = INT(str.find_last_not_of(" "));
  std::string result = str.substr(0, clast + 1);

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::_printf(const char *format, ...)
{
  char buffer[512]; // TODO: get rid of this limit ?
  va_list args;
  va_start (args, format);
  vsnprintf(buffer, sizeof(buffer) - 1, format, args);
  va_end (args);

  return buffer;
}

//==============================================================================
//==============================================================================
//==============================================================================

const std::string Tools::_printf(const std::string &string)
{
  return string;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::size_t Tools::ivec2hash(const IVEC &vec)
{
  std::size_t result = 0;

  for (UINT i = 0; i < vec.n_elem; i++)
  {
    result ^= std::hash<INT> {}(vec(i)) << i; // "<< i" to avoid hash collisions for symmetric vectors
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::size_t Tools::vec2hash(const arma::vec &vec)
{
  std::size_t result = 0;

  for (UINT i = 0; i < vec.n_elem; i++)
  {
    result ^= std::hash<double> {}((double)(vec(i))) << i; // "<< i" to avoid hash collisions for symmetric vectors
  }

  return result;
}

//==============================================================================
//==============================================================================
//==============================================================================

std::string Tools::getLogoStr(void)
{
  std::string result;
  result += "\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m  \ue0ba██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba████████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba█████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m   \x1b[38;2;42;127;255m \x1b[38;2;229;128;255m\ue0ba████████\x1b[48;2;229;128;255m\x1b[38;2;136;0;170m\ue0ba\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0bc\x1b[0;49m\x1b[38;2;128;179;255m██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m█\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0bc\x1b[0;49m\x1b[38;2;128;179;255m██\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m███\x1b[48;2;212;42;255m\x1b[38;2;229;128;255m▁▁▁\x1b[0;49m\x1b[38;2;212;42;255m███\x1b[38;2;136;0;170m█\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m█   \x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;128;179;255m \x1b[38;2;229;128;255m   \x1b[38;2;212;42;255m██████\x1b[38;2;136;0;170m█\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m\ue0bc   \x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m \x1b[38;2;229;128;255m\ue0ba██\x1b[48;2;212;42;255m▁▁▁\x1b[0;49m\x1b[38;2;212;42;255m███\x1b[38;2;136;0;170m█\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█      \x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m█\x1b[0m\n";
  result += "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc      \x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;128;179;255m   \x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m\ue0bc\x1b[0m\n";

  return result;
}
