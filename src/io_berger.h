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

#ifndef IO_BERGER_H
#define IO_BERGER_H

/** \file
 *  \brief Headers for the IOberger class.
 */

#include "global.h"
#include "state.h"
#include "zlib.h"

/// Double nubers are stored in haxadecimal format.
#define ALWAYS_HEXA

//==============================================================================
//==============================================================================
//==============================================================================

/// Key, description,  optional default value and type for each input or output.
#define IO_BERGER_VALID_KEYS \
{ "berger/cvg", "convergence reached", "", "D" }, \
{ "berger/iterx", "number of x-iteractions performed", "", "I" }, \
{ "berger/itert", "number of t-iteractions performed", "", "I" }, \
{ "berger/qsciss", "mean-value of Qneck operator", "", "D" }, \
{ "berger/parent", "parent (starting solution)", "", "S" }, \
{ "berger/enerNeut_Kinetic"         ,"Energy contribution: Neutron Kinetic"         , "" , "D" }, \
{ "berger/enerNeut_Coulomb_Direct"  ,"Energy contribution: Neutron Coulomb_Direct"  , "" , "D" }, \
{ "berger/enerNeut_Central_Direct"  ,"Energy contribution: Neutron Central_Direct"  , "" , "D" }, \
{ "berger/enerNeut_Central_Exchange","Energy contribution: Neutron Central_Exchange", "" , "D" }, \
{ "berger/enerNeut_Spin-orbit"      ,"Energy contribution: Neutron Spin-orbit"      , "" , "D" }, \
{ "berger/enerNeut_Density"         ,"Energy contribution: Neutron Density"         , "" , "D" }, \
{ "berger/enerNeut_Density_D2"      ,"Energy contribution: Neutron Density D2"      , "" , "D" }, \
{ "berger/enerNeut_2-Body_COM_Cor." ,"Energy contribution: Neutron 2-Body_COM_Cor." , "" , "D" }, \
{ "berger/enerNeut_Coulomb_Exchange","Energy contribution: Neutron Coulomb_Exchange", "" , "D" }, \
{ "berger/enerNeut_Pairing_Central" ,"Energy contribution: Neutron Pairing_Central" , "" , "D" }, \
{ "berger/enerNeut_Rearrangement"   ,"Energy contribution: Neutron Rearrangement"   , "" , "D" }, \
{ "berger/enerNeut_Rearrangement_D2","Energy contribution: Neutron Rearrangement D2", "" , "D" }, \
{ "berger/enerProt_Kinetic"         ,"Energy contribution: Proton Kinetic"         , "" , "D" }, \
{ "berger/enerProt_Coulomb_Direct"  ,"Energy contribution: Proton Coulomb_Direct"  , "" , "D" }, \
{ "berger/enerProt_Central_Direct"  ,"Energy contribution: Proton Central_Direct"  , "" , "D" }, \
{ "berger/enerProt_Central_Exchange","Energy contribution: Proton Central_Exchange", "" , "D" }, \
{ "berger/enerProt_Spin-orbit"      ,"Energy contribution: Proton Spin-orbit"      , "" , "D" }, \
{ "berger/enerProt_Density"         ,"Energy contribution: Proton Density"         , "" , "D" }, \
{ "berger/enerProt_Density_D2"      ,"Energy contribution: Proton Density D2"      , "" , "D" }, \
{ "berger/enerProt_2-Body_COM_Cor." ,"Energy contribution: Proton 2-Body_COM_Cor." , "" , "D" }, \
{ "berger/enerProt_Coulomb_Exchange","Energy contribution: Proton Coulomb_Exchange", "" , "D" }, \
{ "berger/enerProt_Pairing_Central" ,"Energy contribution: Proton Pairing_Central" , "" , "D" }, \
{ "berger/enerProt_Rearrangement"   ,"Energy contribution: Proton Rearrangement"   , "" , "D" }, \
{ "berger/enerProt_Rearrangement_D2","Energy contribution: Proton Rearrangement D2", "" , "D" }, \
{ "berger/enerCorrections", "Energy corrections", "", "V" }, \
{ "berger/inertia0"  , "inertia tensor 0", "", "V" }, \
{ "berger/inertia23" , "inertia tensor 2*3", "", "V" }, \
{ "berger/inertia24" , "inertia tensor 2*4", "", "V" }, \
{ "berger/inertia34" , "inertia tensor 3*4", "", "V" }, \
{ "berger/inertia234", "inertia tensor 2*3*4", "", "V" }

//==============================================================================
//==============================================================================
//==============================================================================

/** \brief Interface to berger2ct result files.
 *
 * This class provides methods to open rhl.dat files (results from the berger2ct code).
 */

class IOberger
{
public:

  IOberger(void);                                                      // #TEST#
  ~IOberger();                                                         // #TEST#

  DataTree fromContent(const std::string &content);                    // #TEST#
  //TODO  const std::string info(bool isShort = USE_SHORT_INFO) const;
  static bool checkFileType(const std::string &content);               // #TEST#
  void saveState(const State &,                    // #TEST#
                           const std::string &filename);


  //============================================================================
  //============================================================================
  //============================================================================

private:

  void readBasis(void);
  void readMisc(void);
  void readRhoKappa(void);
  void hexaToDouble(void *doubleVal, INT size, const char *hexaVal);
  INT noeIndex();
  INT noeOpen();
  bool noeRead(const std::string label, const std::string type, INT n, ...);
  INT noeSearch(const std::string label) const;
  void binWrite(std::ofstream &fp, const std::string &mesg = "");
  template<class T> void bufferAppend(T val);
  void bufferAppend(const std::string &);


  //============================================================================
  //============================================================================
  //============================================================================

  /// Buffer to be written.
  char *buffer = NULL;

  /// Size of the buffer (bytes).
  UINT bufferSize = 0;

  /// The content of a rhl.dat file.
  std::string content;

  /// The \f$\alpha\f$ parameter.
  double alpha;

  /// The \f$\beta\f$ parameter.
  double beta;

  /// The \f$d_0\f$ parameter.
  double d_0;

  /// The two-center flag.
  bool v2ct;

  /// The number of oscillator shells.
  INT nOscil;

  /// The maximum value for \f$m\f$.
  INT mxMax = -1;

  /// The maximum imposed value for \f$n_z\f$.
  INT n_zMaxImposed;

  /// The "deformation" of the basis.
  double gq;

  /// The number of labels in the current rhl.dat file.
  INT nbLabels;

  /// The positions of the labels of the current rhl.dat file.
  std::map< std::string, INT > labelPos;

  /// The jdx value.
  INT jdx;

  /// Basis object built from the rhl.dat file.
  Basis basis;

  /// The \f$\rho\f$ matrices.
  Multi<arma::mat> rho;

  /// The \f$\kappa\f$ matrices.
  Multi<arma::mat> kappa;

  /// The \f$D_n\f$ matrix.
  arma::mat matDn;

  /// The \f$D_p\f$ matrix.
  arma::mat matDp;

  /// The HF energies for neutrons.
  arma::vec energyn;

  /// The HF energies for protons.
  arma::vec energyp;

  /// The QP energies for neutrons.
  arma::vec eneQPn;

  /// The QP energies for protons.
  arma::vec eneQPp;

  /// The \f$V_n\f$ vector.
  arma::vec vecVn;

  /// The \f$V_p\f$ vector.
  arma::vec vecVp;

  /// The \f$\Omega_n\f$ vector.
  IVEC vecOmegan;

  /// The \f$\Omega_p\f$ vector.
  IVEC vecOmegap;

  /// The number of protons.
  double zNumber;

  /// The number of neutrons.
  double nNumber;

  /// Chemical potential associated with N/Z constraints
  arma::vec chemPot;

  /// The total binding energy.
  double eneTot;

  /// The convergence.
  double convergence;

  /// Number of iteration to get HFB state
  INT niterx;

  /// Number of iterations to optimize basis within Berger stuff
  INT nitert;

  /// The \f$Q_{\textrm{sciss}}\f$ value.
  double qsciss;

  /// The name of the job.
  std::string jobname;

  /// The number of constraints.
  INT nbcon;

  /// The type of the constraints.
  std::vector<std::string> constraint_type;

  /// The value of the constraints.
  arma::vec constraint_val;

  /// The ponderation of the constraints.
  arma::vec constraint_pond;

  /// The starting iteration of the constraints.
  IVEC constraint_iter1c;

  /// The lagrange multiplier of the constraints.
  arma::vec constraint_lambda;

  /// The energy corrections.
  arma::vec corrections;

  /// The inertia without coupling.
  arma::vec inertia0;

  /// The inertia with Q20 x Q30 coupling.
  arma::vec inertia23;

  /// The inertia with Q20 x Q40 coupling.
  arma::vec inertia24;

  /// The inertia with Q30 x Q40 coupling.
  arma::vec inertia34;

  /// The inertia with Q20 x Q30 x Q40 coupling.
  arma::vec inertia234;

  /// The parent file.
  std::string parent;

  /// The energy contributions [neutrons].
  std::map<std::string, double> energiesNeut;

  /// The energy contributions [protons].
  std::map<std::string, double> energiesProt;
};

#endif // IO_BERGER_H
