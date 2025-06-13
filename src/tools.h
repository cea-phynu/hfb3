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

#ifndef TOOLS_H
#define TOOLS_H

/** \file
 *  \brief Headers for the Tools class.
 */

#include <list>

#include "global.h"
#include "generic.h"
#include "mesh.h"
#include "states.h"

#ifdef CHECK_ACCU
#define ACCU(V) Tools::accu(arma::vectorise(V), __FILE__, __LINE__)
#else
#define ACCU(V) arma::accu(V)
#endif

// template <class T> class Multi;

/// Enable all output messages
#define ALL_OUT msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_DEBUG}

/// Enable all output messages plus call tree
#define ALL_OUT_DBG msgToOut = {MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING, MSG_DEBUG, MSG_TIMELINE}

/// GZIP constants
/// @{
#define MOD_GZIP_ZLIB_WINDOWSIZE 15
#define MOD_GZIP_ZLIB_CFACTOR 9
#define MOD_GZIP_ZLIB_BSIZE 8096
/// @}

/// A shortcut to check equality between armadillo objects.
//#define ARMA_EQUAL(X,Y) ((arma::size(X) == arma::size(Y))&&(arma::accu(X != Y) == 0))

/// Shorcuts for colored Tools::_printf()
/// @{
#define PF(...)         Tools::_printf(__VA_ARGS__)
#define PF_YELLOW(...)  Tools::color("yellow" ) + Tools::_printf(__VA_ARGS__) + Tools::color()
#define PF_GREEN(...)   Tools::color("green"  ) + Tools::_printf(__VA_ARGS__) + Tools::color()
#define PF_MAGENTA(...) Tools::color("magenta") + Tools::_printf(__VA_ARGS__) + Tools::color()
#define PF_BLUE(...)    Tools::color("blue"   ) + Tools::_printf(__VA_ARGS__) + Tools::color()
#define PF_RED(...)     Tools::color("red"    ) + Tools::_printf(__VA_ARGS__) + Tools::color()
#define PF_ORANGE(...)  Tools::color("orange" ) + Tools::_printf(__VA_ARGS__) + Tools::color()
/// @}

/// A shorcut for Tools::info(Tools::_printf()).
#define INFO(...) Tools::info(PF(__VA_ARGS__))

/// A debugging shorcut
#define DEBUG(...) Tools::info(PF_RED("=== DEBUG %s:%d ===\n", __FILE__, __LINE__) + Tools::_printf(__VA_ARGS__))

/// Print colored separation lines.
/// @{
#define PFLR INFO(Tools::color("red"  ) + "================================================================================")
#define PFLG INFO(Tools::color("green") + "================================================================================")
#define PFLB INFO(Tools::color("blue" ) + "================================================================================")
#define PFLW INFO(Tools::color("white") + "================================================================================")
/// @}

#ifdef NO_DBG_STACK

#define ASSERT(A, ...) {}
#define DBG_ENTER {}
#define DBG_LEAVE {return;}
#define DBG_RETURN(A) return A

#else

#define ASSERT(A, ...) {if (!(A)) Tools::_error(PF_RED(__VA_ARGS__) + PF(" [%s:%d]", __FILE__, __LINE__));}
#define DBG_ENTER Tools::timer(PF("%s", __PRETTY_FUNCTION__), PF("[%s:%d]", __FILE__, __LINE__))
#define DBG_LEAVE {Tools::timerEnd(__PRETTY_FUNCTION__); return;}
#define DBG_RETURN(A) {Tools::timerEnd(__PRETTY_FUNCTION__); return A;}

#endif

/// Error macro define
#define ERROR(...) {Tools::_error(PF_RED(__VA_ARGS__) + PF(" [%s:%d]", __FILE__, __LINE__));}

/// Hash `char*` data.
constexpr UINT hash2(const char* data, size_t const size) {
    uint32_t hash = 5381;

    for(const char *c = data; c < data + size; ++c)
        hash = ((hash << 5) + hash) + (unsigned char) *c;

    return hash;
}

/// Hack to allow `switch(std::string)`.
constexpr UINT operator "" _h(const char* data, size_t const size) {
  return hash2(data, size);
}

/// Hash `std::string` data.
UINT hash2(std::string s);

//========================
//========================
//========================

/// An entry in the stack.
struct StackItem
{
  std::string name;
  std::string location;
  std::string color;
};

//========================
//========================
//========================

/// A class containing various tools (not a namespace because swig ignores namespaces).
class Tools : public Generic
{
public:

  /// Enum for table characters.
  enum {CHAR_TB,
        CHAR_RL,
        CHAR_TR,
        CHAR_RB,
        CHAR_BL,
        CHAR_TL,
        CHAR_TRB,
        CHAR_TRL,
        CHAR_TBL,
        CHAR_RBL,
        CHAR_TRBL,
        CHAR_NB
       };

  //========================
  //=== numerical checks ===
  //========================

  static double accu(arma::vec in, const char*, INT lineNumber);

  //========================
  //=== string functions ===
  //========================

  /// Remove whitespaces from the start and the end of a string.
  static const std::string trim(const std::string);                   // #TEST#

  /// Remove whitespaces from the start of a string.
  static const std::string trim_b(const std::string);                 // #TEST#

  /// Remove whitespaces from the end of a string.
  static const std::string trim_e(const std::string);                 // #TEST#

  /// Split a std::string into a std::vector<std::string> using a separator.
  static std::vector<std::string> stringSplit(const std::string &,    // #TEST#
      char);

  /// Construct a std::string from a printf formated char*.
  static const std::string _printf(const char *format, ...);          // #TEST#

  /// Construct a std::string from a std::string.
  static const std::string _printf(const std::string &mesg);          // #TEST#

  /// Find the index of an occurence of an item in an std container.
  template<typename T, typename A>
  static UINT index(std::vector<T, A> const &cont, T item)            // #TEST#
  {
    auto it = std::find(cont.begin(), cont.end(), item);

    if (it == cont.end()) return -1;

    return std::distance(cont.begin(), it);
  }

  static std::string humanSize(UINT bytes);

  //=========================
  //=== armadillo objects ===
  //=========================

  /// Insert value at the end
  static void growIVec(IVEC &m, INT v);                               // #TEST#

  /// Return a uniform (0 < v < 1) random matrix.
  static const arma::mat randMat(UINT nbr, UINT nbc);                 // #TEST#

  /// Convert an arma::mat object to .txt format.
  static std::string matToTxt(const arma::mat &m);                    // #TEST#

  /// Convert an arma::mat object to .csv format.
  static std::string matToCsv(const arma::mat &m);                    // #TEST#

  /// Print an arma::mat object as dot/crosses.
  static std::string matToCrosses(const arma::mat &m,                 // #TEST#
                                  double epsilon = 1e-12);

  /// Convert an arma::cube object to .df3 format.
  static std::string cubeToDf3(const arma::cube &m);                  // #TEST#

  /// Convert an arma::mat object to .vtk format.
  static std::string matToVTK(const std::string &name,                // #TEST#
                              const arma::mat &m,
                              const Mesh &mesh);

  /// Convert an arma::cube object to .vtk format.
  static std::string cubeToVTK(const std::string &name,               // #TEST#
                               const arma::cube &m,
                               const Mesh &mesh);

  /// Convert an arma::cube object to .raw format (usable by POVray).
  static std::string cubeToRaw(const arma::cube &m);                  // #TEST#

  /// Construct an arma::mat from an arma::vec.
  static arma::mat matFromCol(const arma::vec &v, UINT nb = 0);       // #TEST#

  /// Construct an arma::mat from an arma::rowvec.
  static arma::mat matFromRow(const arma::rowvec &v, UINT nb = 0);    // #TEST#

  /// Calculate the element-wise product of an arma::mat and an arma::vec.
  static arma::cube matTimesVec(const arma::mat &m,                   // #TEST#
                                const arma::vec &v);

  /// Calculate the "product" of an arma::cube and an arma::vec.
  static arma::cube cubeTimesVec(const arma::cube &c,                 // #TEST#
                                 const arma::vec &v);

  /// Calculate the "product" of an arma::cube and an arma::mat.
  static arma::cube cubeTimesMat(const arma::cube &c,                 // #TEST#
                                 const arma::mat &m);

  /// Extract an arma::mat from an arma::cube.
  static arma::mat cubeToMat(const arma::cube &c, INT iaxis, INT n);  // #TEST#

  /// Print an arma::mat object.
  static void info(const std::string &mesg,                           // #TEST#
                   const arma::mat &m,
                   bool dump = false);

  static const std::string matToStr(const std::string &mesg,          // #TEST#
                                    const arma::mat &m,
                                    bool dump = false,
                                    const std::vector<std::string> lineLabels = {},
                                    const std::vector<std::string> colLabels = {}
                                   );

  /// Print an IMAT object.
  static void info(const std::string &mesg,                           // #TEST#
                   const IMAT &m,
                   bool dump = false);

  /// Print an UMAT object.
  static void info(const std::string &mesg,                           // #TEST#
                   const UMAT &m,
                   bool dump = false);

  /// Print an arma::cube object.
  static void info(const std::string &mesg,                           // #TEST#
                   const arma::cube &c,
                   bool dump = false);

  /// Print an ICUBE object.
  static void info(const std::string &mesg,                           // #TEST#
                   const ICUBE &c,
                   bool dump = false);

  /// Print an 'X' and '.' representation of an arma::mat object.
  static std::string matPrint(const arma::mat &m);                    // #TEST#

  /// Evaluate a gaussian.
  static arma::vec gaussian(const arma::vec &g, const arma::vec &x);  // #TEST#

  /// Calculate MD5 sum.
  static std::string armaMD5(const arma::vec   &mat);                 // #TEST#
  static std::string armaMD5(const arma::mat   &mat);                 // #TEST#
  static std::string armaMD5(const arma::cube  &mat);                 // #TEST#

  /// Calculate hash.
  static std::size_t vec2hash(const arma::vec &vec);                  // #TEST#
  static std::size_t ivec2hash(const IVEC &vec);                // #TEST#

  /// Calculate the pfaffian of a skew-symmetric matrix m
  static double getPfaffian(const arma::mat m);                       // #TEST#

  /// Construct an std::vector<double> from an arma::vec.
  static const std::vector<double> vecToStd(const arma::vec &v);      // #TEST#

  /// Construct an arma::cube from an arma::vec.
  static arma::cube vecToCube(arma::vec X,                            // #TEST#
                              UINT nRows,
                              UINT nCols);

  /// Construct an arma::cube from an arma::mat.
  static arma::cube matToCube(arma::mat m,                            // #TEST#
                              UINT nSlices);

  /// Return a string representation of a vector (example: "[1.2, 3.4]")
  static std::string vecToStr(const arma::vec &vec);                  // #TEST#

  /// Return a string representation of a vector of integers
  static std::string ivecToStr(const IVEC &vec);                 // #TEST#

  /// Check that a matrix is symmetric
  static bool checkSymmetry(const arma::mat &mat,                     // #TEST#
                            const std::string &mesg);

  /// More robust implementation of arma::eig_sym().
  static bool eig_sym(arma::vec &val, arma::mat &vec,                 // #TEST#
                      const arma::mat &mat);

  /// Find a column in a matrix
  static INT findColInMat(const arma::vec &v, const arma::mat &m);    // #TEST#

  //=================================
  //=== Timer / call tree builder ===
  //=================================

  /// Get date and time.
  static std::string date(void);                                      // #TEST#

  /// Wall clock.
  static double clock(void);                                          // #TEST#

  /// Set a timer.
  static void timer(const std::string &label = "",                    // #TEST#
                    const std::string &location = "",
                    const std::string &color = "");

  /// End of a function.
  static bool timerEnd(const std::string &label = "",                 // #TEST#
                       const std::string &msg = "");

  /// Progress bar.
  static void progress(const std::string &, INT, INT);                // #TEST#

  /// Progress bar.
  static void progress(const std::string &, INT);                     // #TEST#

  //============
  //=== zlib ===
  //============

  /// Compress a std::string.
  static std::string compressString(const std::string &str);          // #TEST#

  /// Decompress a std::string.
  static std::string decompressString(const std::string &);           // #TEST#

  /// Read a file.
  static std::string readFile(const std::string &filename);           // #TEST#

  /// Read std::cin.
  static std::string readCin(void);                                   // #TEST#

  /// Save content to file.
  static std::string save(const std::string &result,                  // #TEST#
                          const std::string &filename,
                          bool compress = false);

  //============
  //=== Misc ===
  //============

  /// Get the name of the isospin.
  static std::string strIsospin(INT iso);

  /// Get a random integer < imax.
  static INT irand(INT imax);

  /// Get a random double < dmax.
  static double drand(double dmax = 1.0);

  /// Specify a seed for the RNG.
  static void setSeed(INT seed);

  /// Exit the program.
  static void end(INT c = 0);

  /// Print a table.
  static std::string printTable(std::string, INT style = -1);

  /// Print a boxed message.
  static std::string boxed(std::string input);

  /// Print tabulated values.
  static std::string valueTable(const std::string &title,
                                const std::list<std::string> &columnLabel,
                                const std::list<std::string> &unit,
                                const std::list<std::list<std::string> > &values);

  /// Return special characters for CLI drawing.
  static const std::string getChar(INT id);

  /// Get the length (number of characters) of a UTF-8 string.
  static UINT stringLen(const std::string &s);

  /// Get random color name.
  static const std::string randomColor(INT fixedColor = -1);

  /// Get ANSI color code.
  static std::string color(const std::string &c = "");

  /// Get ANSI truecolor sequence.
  static std::string trueColor(INT r, INT g, INT b);

  /// Print on standard output.
  static void printOut(const std::string &mesg);

  /// Print a debug message.
  static void debug(const std::string &mesg);

  /// Print an important message.
  static void mesg(const std::string &cat, const std::string &mesg);

  /// Print an info message.
  static void info(const std::string &mesg);

  /// Print a warning.
  static void warning(const std::string &mesg);

  /// Print an error.
  static void _error(const std::string &mesg);

  /// Calculate the power of integers.
  static INT ipow(INT base, INT exp);

  /// Return a version string.
  static const std::string version(void);

  /// Return a tree representation of strings.
  static const std::string treeStr(const std::vector<std::pair<std::string, std::string> >, bool isShort = false);

  /// Return the call stack.
  static const std::string getStack(void);

  /// Get the HFB3 logo string.
  static std::string getLogoStr(void);

  /// Return a std::string representation of various object types.
  static const std::string infoStr(const std::vector<std::string> &list);
  static const std::string infoStr(const INT &, bool isShort = false);
  static const std::string infoStr(const UINT &, bool isShort = false);
  static const std::string infoStr(const double &, bool isShort = false);
  static const std::string infoStr(const std::string &, bool isShort = false);
  static const std::string infoStr(const arma::vec &, bool isShort = false);
  static const std::string infoStr(const arma::mat &, bool isShort = false);
  static const std::string infoStr(const arma::cube &, bool isShort = false);
  static const std::string infoStr(const IVEC &, bool isShort = false);
  static const std::string infoStr(const UVEC &, bool isShort = false);
  static const std::string infoStr(const IMAT &, bool isShort = false);
  static const std::string infoStr(const UMAT &, bool isShort = false);
  static const std::string infoStr(const ICUBE &, bool isShort = false);
  static const std::string infoStr(const UCUBE &, bool isShort = false);
  static const std::string infoStr(const States &, bool isShort = false);
  static const std::string infoStr(const bool &);
};

/// Multiplication operator between `std::string` and integer.
std::string operator*(const std::string &s, const INT &n);

#endif // TOOLS_H
