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

//==============================================================================
//==============================================================================
//==============================================================================

%module(directors="1") hfb3
%include "typemaps.i"
%include "stl.i"
%include "std_map.i"
%include "std_unordered_map.i"
%include "std_vector.i"
%include "std_set.i"
%include "std_list.i"
%include "std_string.i"
%include "std_shared_ptr.i"
%include "exception.i"

//==============================================================================
//==============================================================================
//==============================================================================

%exception
{
  try
  {
    $action
  }
  catch (const std::runtime_error& e) {
    SWIG_exception(SWIG_RuntimeError, e.what());
  }
  catch (const std::logic_error& e) {
    SWIG_exception(SWIG_IndexError, e.what());
  }
  catch (const std::string s) {
    SWIG_exception(SWIG_RuntimeError, s.c_str());
  }
}

//==============================================================================
//==============================================================================
//==============================================================================

%{
#define SWIG_FILE_WITH_INIT

#include "hfb3.h"
#include <memory>

using namespace std;
%}

//==============================================================================
//==============================================================================
//==============================================================================

%include "armanpy.i"
%fragment("armanpy_vec_typemaps");
%fragment("armanpy_mat_typemaps");
%fragment("armanpy_cube_typemaps");

%copyctor State;
%copyctor Basis;
%copyctor Interaction;
%copyctor Field;
%copyctor Constraint;
%copyctor States;

//==============================================================================
//==============================================================================
//==============================================================================

%shared_ptr(Generic)
%shared_ptr(Action)
%shared_ptr(Axis)
%shared_ptr(Base64)
%shared_ptr(Basis)
%shared_ptr(CliParser)
%shared_ptr(DataTree)
%shared_ptr(Discrete)
%shared_ptr(Field)
%shared_ptr(FieldCDM2)
%shared_ptr(FieldCentral)
%shared_ptr(FieldCoulombSlater)
%shared_ptr(FieldCoulombExact)
%shared_ptr(FieldDensity)
%shared_ptr(FieldDensityFR)
%shared_ptr(FieldKinetic)
%shared_ptr(FieldRearrangement)
%shared_ptr(FieldRearrangementFR)
%shared_ptr(FieldSpinOrbit)
%shared_ptr(FieldWS)
%shared_ptr(Fragments)
%shared_ptr(GeometricalOperators)
%shared_ptr(Geometry)
%shared_ptr(GradientWalk)
%shared_ptr(States)
%shared_ptr(Interaction)
%shared_ptr(IOamedee)
%shared_ptr(IOberger)
%shared_ptr(IOhfb3)
%shared_ptr(IOjson)
%shared_ptr(IOmsgp)
%shared_ptr(Mesh)
%shared_ptr(Mixing)
%shared_ptr(MultipoleOperators)
/* %shared_ptr(Multi) */
%shared_ptr(Plot)
%shared_ptr(Qnumbers)
%shared_ptr(State)
%shared_ptr(Solver)
%shared_ptr(SolverBasis)
%shared_ptr(SolverHFBBroyden)
%shared_ptr(SolverHFBGradient)
%shared_ptr(SolverWS)
%shared_ptr(System)
%shared_ptr(WSPot)
%shared_ptr(Tools)


//==============================================================================
//==============================================================================
//==============================================================================

namespace std
{
  %template(vec_qnumbers) vector<Qnumbers>;
  %template(map_constraints) map<string, Constraint>;
  %template(FieldMap) map<int, shared_ptr<Field>>;
  %template(FieldVec) vector<shared_ptr<Field>>;
  %template(map_string_int) map<string, int>;
  %template(map_string_double) map<string, double>;
  %template(map_string_string) map<string, string>;
  %template(vec_string) vector<string>;
  %template(vec_double) vector<double>;
  %template(set_int) set<int>;
}

%copyctor map_string_double;

//==============================================================================
//==============================================================================
//==============================================================================

// Ignored methods and properties. Feel free to un-ignore if needed.

%ignore DataTree::get(arma::mat &, const std::string&) const;
%ignore DataTree::get(arma::vec &, const std::string&) const;
%ignore DataTree::get(arma::vec &, const std::string&) const;
%ignore DataTree::get(arma::cube &, const std::string&) const;
%ignore DataTree::get(IVEC &, const std::string&) const;
%ignore DataTree::get(IMAT &, const std::string&) const;
%ignore DataTree::get(ICUBE &, const std::string&) const;
%ignore DataTree::get(UVEC &, const std::string&) const;
%ignore DataTree::get(UMAT &, const std::string&) const;
%ignore DataTree::get(UCUBE &, const std::string&) const;
%ignore DataTree::get(Multi<arma::mat> &, const std::string&) const;
%ignore DataTree::get(Multi<arma::vec> &, const std::string&) const;
%ignore DataTree::get(Multi<arma::vec> &, const std::string&) const;
%ignore DataTree::get(Multi<arma::cube> &, const std::string&) const;
%ignore DataTree::get(Multi<IVEC > &, const std::string&) const;
%ignore DataTree::get(Multi<IMAT > &, const std::string&) const;
%ignore DataTree::get(Multi<ICUBE > &, const std::string&) const;
%ignore DataTree::get(Multi<UVEC > &, const std::string&) const;
%ignore DataTree::get(Multi<UMAT > &, const std::string&) const;
%ignore DataTree::get(Multi<UCUBE > &, const std::string&) const;

%ignore DataTree::get(arma::mat &, const std::string&, bool) const;
%ignore DataTree::get(arma::vec &, const std::string&, bool) const;
%ignore DataTree::get(arma::vec &, const std::string&, bool) const;
%ignore DataTree::get(arma::cube &, const std::string&, bool) const;
%ignore DataTree::get(IVEC &, const std::string&, bool) const;
%ignore DataTree::get(IMAT &, const std::string&, bool) const;
%ignore DataTree::get(ICUBE &, const std::string&, bool) const;
%ignore DataTree::get(UVEC &, const std::string&, bool) const;
%ignore DataTree::get(UMAT &, const std::string&, bool) const;
%ignore DataTree::get(UCUBE &, const std::string&, bool) const;
%ignore DataTree::get(Multi<arma::mat> &, const std::string&, bool) const;
%ignore DataTree::get(Multi<arma::vec> &, const std::string&, bool) const;
%ignore DataTree::get(Multi<arma::vec> &, const std::string&, bool) const;
%ignore DataTree::get(Multi<arma::cube> &, const std::string&, bool) const;
%ignore DataTree::get(Multi<IVEC > &, const std::string&, bool) const;
%ignore DataTree::get(Multi<IMAT > &, const std::string&, bool) const;
%ignore DataTree::get(Multi<ICUBE > &, const std::string&, bool) const;
%ignore DataTree::get(Multi<UVEC > &, const std::string&, bool) const;
%ignore DataTree::get(Multi<UMAT > &, const std::string&, bool) const;
%ignore DataTree::get(Multi<UCUBE > &, const std::string&, bool) const;

%ignore DataTree::set(const std::string&, const std::string&);
%ignore DataTree::set(const std::string&, const char*);
%ignore DataTree::set(const std::string&, INT);
%ignore DataTree::set(const std::string&, double);
%ignore DataTree::set(const std::string&, const arma::vec&);
%ignore DataTree::set(const std::string&, const arma::mat&);
%ignore DataTree::set(const std::string&, const arma::cube&);
%ignore DataTree::set(const std::string&, const IVEC&);
%ignore DataTree::set(const std::string&, const IMAT&);
%ignore DataTree::set(const std::string&, const ICUBE&);
%ignore DataTree::set(const std::string&, const UVEC&);
%ignore DataTree::set(const std::string&, const UMAT&);
%ignore DataTree::set(const std::string&, const UCUBE&);
%ignore DataTree::set(const std::string&, const Multi<UVEC>&);
%ignore DataTree::set(const std::string&, const Multi<UMAT>&);
%ignore DataTree::set(const std::string&, const Multi<UCUBE>&);
%ignore DataTree::set(const std::string&, const Multi<IVEC>&);
%ignore DataTree::set(const std::string&, const Multi<IMAT>&);
%ignore DataTree::set(const std::string&, const Multi<ICUBE>&);
%ignore DataTree::set(const std::string&, const Multi<arma::vec>&);
%ignore DataTree::set(const std::string&, const Multi<arma::mat>&);
%ignore DataTree::set(const std::string&, const Multi<arma::cube>&);
%ignore DataTree::set(const std::string&, const bool&);
%ignore DataTree::set(const std::string&);

%ignore Action::BlockingTry;

%ignore operator*(const std::string&, const INT&);

%ignore Interaction::operator()(INT) const;
%ignore operator "" _h(const char*, size_t const);

%ignore operator==(const UVEC &, const UVEC &);
%ignore operator==(const IVEC &, const IVEC &);
%ignore operator==(const arma::vec &, const arma::vec &);
%ignore operator==(const UMAT &, const UMAT &);
%ignore operator==(const IMAT &, const IMAT &);
%ignore operator==(const arma::mat &, const arma::mat &);
%ignore operator==(const UCUBE &, const UCUBE &);
%ignore operator==(const ICUBE &, const ICUBE &);
%ignore operator==(const arma::cube &, const arma::cube &);
%ignore operator!=(const UVEC &, const UVEC &);
%ignore operator!=(const IVEC &, const IVEC &);
%ignore operator!=(const arma::vec &, const arma::vec &);
%ignore operator!=(const UMAT &, const UMAT &);
%ignore operator!=(const IMAT &, const IMAT &);
%ignore operator!=(const arma::mat &, const arma::mat &);
%ignore operator!=(const UCUBE &, const UCUBE &);
%ignore operator!=(const ICUBE &, const ICUBE &);
%ignore operator!=(const arma::cube &, const arma::cube &);

%ignore operator==(const States &, const States &);
%ignore operator==(const State &, const State &);
%ignore operator==(const IndivState &, const IndivState &);
%ignore operator!=(const States &, const States &);
%ignore operator!=(const State &, const State &);
%ignore operator!=(const IndivState &, const IndivState &);

%ignore Tools::armaMD5(const arma::vec &);
%ignore Tools::armaMD5(const arma::mat &);
%ignore Tools::armaMD5(const arma::cube &);

%ignore Tools::infoStr(const std::vector<std::string> &);
%ignore Tools::infoStr(const INT &, bool);
%ignore Tools::infoStr(const UINT &, bool);
%ignore Tools::infoStr(const double &, bool);
%ignore Tools::infoStr(const std::string &, bool);
%ignore Tools::infoStr(const arma::vec &, bool);
%ignore Tools::infoStr(const arma::mat &, bool);
%ignore Tools::infoStr(const arma::cube &, bool);
%ignore Tools::infoStr(const IVEC &, bool);
%ignore Tools::infoStr(const UVEC &, bool);
%ignore Tools::infoStr(const IMAT &, bool);
%ignore Tools::infoStr(const UMAT &, bool);
%ignore Tools::infoStr(const ICUBE &, bool);
%ignore Tools::infoStr(const UCUBE &, bool);
%ignore Tools::infoStr(const States &, bool);
%ignore Tools::infoStr(const bool &);

%ignore Tools::infoStr(const INT &);
%ignore Tools::infoStr(const UINT &);
%ignore Tools::infoStr(const double &);
%ignore Tools::infoStr(const std::string &);
%ignore Tools::infoStr(const arma::vec &);
%ignore Tools::infoStr(const arma::mat &);
%ignore Tools::infoStr(const arma::cube &);
%ignore Tools::infoStr(const IVEC &);
%ignore Tools::infoStr(const UVEC &);
%ignore Tools::infoStr(const IMAT &);
%ignore Tools::infoStr(const UMAT &);
%ignore Tools::infoStr(const ICUBE &);
%ignore Tools::infoStr(const UCUBE &);
%ignore Tools::infoStr(const States &);
%ignore Tools::infoStr(const bool &);

%ignore Tools::info(std::string const &, arma::mat const &             );
%ignore Tools::info(std::string const &, arma::cube const &            );

%ignore Tools::info(std::string const &, IMAT const &                  );
%ignore Tools::info(std::string const &, ICUBE const &                 );

%ignore Tools::info(std::string const &, arma::mat const & , bool      );
%ignore Tools::info(std::string const &, arma::cube const &, bool      );

%ignore Tools::info(std::string const &, IMAT const &,       bool      );
%ignore Tools::info(std::string const &, ICUBE const &,      bool      );

//==============================================================================
//==============================================================================
//==============================================================================

// renamed methods / properties.

%rename(plotMap)  Plot::map(const std::string &, const arma::mat &, const Mesh &);
%rename(plotMap)  Plot::map(const std::string &, const arma::mat &, const Mesh &, const std::string &);
%rename(plotMap)  Plot::map(const std::string &, const arma::mat &, const Mesh &, const std::string &, const std::string &);
%rename(lbd)  Constraint::lambda;

%rename("del_key") del;
%rename("to_stream") operator<<;

//==============================================================================
//==============================================================================
//==============================================================================

%define PY_INFO_TO_REPR(TypeName)
%extend TypeName {
  %pythoncode %{
    def __repr__(self):
      return self.info()
    def show(self):
      pass
  %}
}
%enddef

//==============================================================================

PY_INFO_TO_REPR(Action)
PY_INFO_TO_REPR(Axis)
PY_INFO_TO_REPR(CliParser)
PY_INFO_TO_REPR(Constraint)
PY_INFO_TO_REPR(DataTree)
PY_INFO_TO_REPR(Discrete)
PY_INFO_TO_REPR(Field)
PY_INFO_TO_REPR(Fragments)
PY_INFO_TO_REPR(Generic)
PY_INFO_TO_REPR(GeometricalOperators)
PY_INFO_TO_REPR(Geometry)
PY_INFO_TO_REPR(States)
PY_INFO_TO_REPR(Interaction)
PY_INFO_TO_REPR(IOamedee)
PY_INFO_TO_REPR(IOberger)
PY_INFO_TO_REPR(IOjson)
PY_INFO_TO_REPR(IOmsgp)
PY_INFO_TO_REPR(Mesh)
PY_INFO_TO_REPR(Mixing)
PY_INFO_TO_REPR(MultipoleOperators)
PY_INFO_TO_REPR(Plot)
PY_INFO_TO_REPR(Qnumbers)
PY_INFO_TO_REPR(State)
PY_INFO_TO_REPR(SolverBasis)
PY_INFO_TO_REPR(SolverHFBBroyden)
PY_INFO_TO_REPR(SolverHFBGradient)
PY_INFO_TO_REPR(System)
PY_INFO_TO_REPR(Tools)


//==============================================================================
//==============================================================================
//==============================================================================

%feature("director") Callback;

//==============================================================================
//==============================================================================
//==============================================================================

// %template(armaIV) arma::Col<int>;

%template(vector_int) std::vector<int>;
%template(vector_lint) std::vector<long long int>;
%template(vector_uint) std::vector<long long unsigned int>;
%template(vector_vector_int) std::vector<std::vector<int> >;
%template(vector_vector_lint) std::vector<std::vector<long long int> >;
%template(vector_vector_uint) std::vector<std::vector<long long unsigned int> >;

%template(umap_array_int) std::unordered_map   <std::vector<int>, int, myHasher>;
%template(umap_array_double) std::unordered_map<std::vector<int>, double, myHasher>;
%template(umap_array_string) std::unordered_map<std::vector<int>, std::string, myHasher>;
%template(umap_array_vec) std::unordered_map   <std::vector<int>, arma::vec, myHasher>;
%template(umap_array_mat) std::unordered_map   <std::vector<int>, arma::mat, myHasher>;
%template(umap_array_cube) std::unordered_map  <std::vector<int>, arma::cube, myHasher>;
%template(umap_array_states) std::unordered_map<std::vector<int>, States, myHasher>;

%template(umap_array_ivec) std::unordered_map  <std::vector<int>, arma::Col<int>, myHasher>;
%template(umap_array_uvec) std::unordered_map  <std::vector<int>, arma::Col<long long unsigned int>, myHasher>;
%template(umap_array_imat) std::unordered_map  <std::vector<int>, arma::Mat<int>, myHasher>;
%template(umap_array_umat) std::unordered_map  <std::vector<int>, arma::Mat<long long unsigned int>, myHasher>;
%template(umap_array_icube) std::unordered_map <std::vector<int>, arma::Cube<int>, myHasher>;
%template(umap_array_ucube) std::unordered_map <std::vector<int>, arma::Cube<long long unsigned int>, myHasher>;

%template(vector_states) std::vector<States>;
%template(vector_geometry) std::vector<Geometry>;

//==============================================================================
//==============================================================================
//==============================================================================

%include "generic.h"
%include "action.h"
%include "axis.h"
%include "base64.h"
%include "basis.h"
%include "cliparser.h"
%include "constraint.h"
%include "datatree.h"
%include "discrete.h"
%include "field.h"
%include "field_cdm2.h"
%include "field_central.h"
%include "field_coulomb_slater.h"
%include "field_coulomb_exact.h"
%include "field_density.h"
%include "field_density_fr.h"
%include "field_kinetic.h"
%include "field_rearrangement.h"
%include "field_rearrangement_fr.h"
%include "field_spin_orbit.h"
%include "field_ws.h"
%include "fragments.h"
%include "gradientwalk.h"
%include "geometrical_operators.h"
%include "geometry.h"
%include "global.h"
%include "states.h"
%include "interaction.h"
%include "io_amedee.h"
%include "io_berger.h"
%include "io_hfb3.h"
%include "io_json.h"
%include "io_msgp.h"
%include "logger.h"
%include "mesh.h"
%include "mixing.h"
%include "multi.h"
%include "multipole_operators.h"
%include "plot.h"
%include "qnumbers.h"
%include "state.h"
%include "solver.h"
%include "solver_hfb_broyden.h"
%include "solver_hfb_gradient.h"
%include "solver_basis.h"
%include "solver_ws.h"
%include "system.h"
%include "tools.h"
%include "wspot.h"


//==============================================================================
//==============================================================================
//==============================================================================

%define PY_ADD_CALL_METHOD(TypeName)
%extend TypeName {
  %pythoncode %{
    def __call__(self, *vals):
      return self.__getitem__([*vals])
  %}
}
%enddef

PY_ADD_CALL_METHOD(Multi)

//==============================================================================
//==============================================================================
//==============================================================================

%rename(ArmaIV) arma::Col<int>;

%template(MultiI     ) Multi<int   >;
%template(MultiD     ) Multi<double>;
%template(MultiS     ) Multi<string>;
%template(MultiV     ) Multi<arma::vec>;
%template(MultiIV    ) Multi<IVEC  >;
%template(MultiUV    ) Multi<UVEC  >;
%template(MultiM     ) Multi<arma::mat>;
%template(MultiIM    ) Multi<IMAT  >;
%template(MultiUM    ) Multi<UMAT  >;
%template(MultiC     ) Multi<arma::cube>;
%template(MultiIC    ) Multi<ICUBE >;
%template(MultiUC    ) Multi<UCUBE >;
%template(MultiStates) Multi<States>;

%template(VectorState) std::vector<State>;

// %template(MultiGeometry1) Multi<Geometry>;
// %template(MultiQnumbers1) Multi<Qnumbers>;

//==============================================================================
//==============================================================================
//==============================================================================

%pythoncode %{
class PyCallback(Callback):
  def __init__(self):
    Callback.__init__(self)
  def run(self, str):
    #print(str, end = " ", flush = True)
    #print(str)
    print(str, flush = True)

cvar.useColors = True
cvar.msgToOut = [MSG_ERROR, MSG_MAIN, MSG_INFO, MSG_WARNING]
cvar.exitOnError = False

cvar.logger.setCallback(PyCallback().__disown__())

#cvar.logger.log(Tools.version())
#cvar.logger.log("Python logger callback activated")
%}

//==============================================================================
//==============================================================================
//==============================================================================

%pythoncode %{
import numpy as np

def dictToDataTree(ldict):

  result = DataTree()
  for v in ldict.items():
    if type(v[1]) is int:
      result.setI(v[0], v[1])
    elif type(v[1]) is float:
      result.setD(v[0], v[1])
    elif type(v[1]) is str:
      result.setS(v[0], v[1])
    elif type(v[1]) is np.ndarray:
      if v[1].dtype == np.float64:
        if v[1].ndim == 1:
          result.setV(v[0], v[1])
        if v[1].ndim == 2:
          result.setM(v[0], v[1])
        if v[1].ndim == 3:
          result.setC(v[0], v[1])
      if v[1].dtype == np.int32:
        if v[1].ndim == 1:
          result.setIV(v[0], v[1])
        if v[1].ndim == 2:
          result.setIM(v[0], v[1])
        if v[1].ndim == 3:
          result.setIC(v[0], v[1])
      if v[1].dtype == np.uint64:
        if v[1].ndim == 1:
          result.setUV(v[0], v[1])
        if v[1].ndim == 2:
          result.setUM(v[0], v[1])
        if v[1].ndim == 3:
          result.setUC(v[0], v[1])
    elif v[1] is None:
      result.setE(v[0])
  return result
%}

//==============================================================================
//==============================================================================
//==============================================================================

%{
DataTree bytesToDataTree(char *bytes, int length)
{
  DataTree result = DataTree::fromContent(std::string(bytes, length));
  return result;
}

PyObject* dataTreeToBytes(const DataTree &dataTree)
{
  std::string serialized = IOmsgp::serializeDataTree(dataTree);
  return PyBytes_FromStringAndSize(serialized.c_str(), serialized.size());
}
%}

//==============================================================================
//==============================================================================
//==============================================================================

%typemap(in) (char *bytes, int length) {
    Py_ssize_t len;
    PyBytes_AsStringAndSize($input, &$1, &len);
    $2 = (int)len;
}

//==============================================================================
//==============================================================================
//==============================================================================

DataTree bytesToDataTree(char *bytes, int length);
PyObject* dataTreeToBytes(const DataTree &dataTree);

