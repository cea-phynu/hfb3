#!/usr/bin/env python3

##==============================================================================
## HFB3
## Copyright CEA, DAM F-91297 Arpajon, France
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================

import os
import numpy
from setuptools import setup, Extension

#===============================================================================
#===============================================================================
#===============================================================================

hfb3_dir = "../.."
hfb3_build_dir = hfb3_dir + "/lib"
hfb3_include_dir = hfb3_dir + "/src"

#===============================================================================
#===============================================================================
#===============================================================================

# fix an annoying bug with setuptools
# https://stackoverflow.com/questions/29477298/setup-py-run-build-ext-before-anything-else/48942866#48942866

from setuptools.command.build_py import build_py as _build_py

class build_py(_build_py):
    def run(self):
        self.run_command("build_ext")
        return super().run()

#===============================================================================
#===============================================================================
#===============================================================================

def get_env(key):
  try:
    v = os.environ[key]
  except:
    v = ''
  return v

def get_skill_opts():
  skill = get_env("SKILL")
  if not skill:
    skill = "SKILL3"

  if (skill == "SKILL0"):
    return ["-O0", "-ggdb3", "-fno-inline-functions", "-D_GLIBCXX_DEBUG", "-DCHECKBOUNDS_FMULTI", "-DCHECK_ACCU"]
  if (skill == "SKILL1"):
    return ["-O1", "-g", "-march=native", "-DCHECKBOUNDS_FMULTI", "-DCHECK_ACCU"]
  if (skill == "SKILL2"):
    return ["-O2", "-g", "-march=native", "-DCHECKBOUNDS_FMULTI"]
  if (skill == "SKILL3"):
    return ["-O3", "-g", "-march=native", "-DARMA_NO_DEBUG"]
  if (skill == "SKILL4"):
    return ["-O3", "-g", "-march=native", "-DARMA_NO_DEBUG", "-DNO_DBG_STACK"]

  raise NotImplementedError


include_dirs = [hfb3_include_dir,
              get_env('PYTHON3_ROOT'),
              '../deps/msgpack-c/include',
              './include/armanpy/',
              get_env('ARMADILLO_INCDIR'),
              get_env('BOOST_INCDIR'),
              # get_env('BLIS_AOCL_INCDIR'),
              numpy.get_include()
             ]

include_dirs = [d for d in include_dirs if d != ""]


libraries = ['m', 'z', 'armadillo' ]
# libraries = ['m', 'z', 'lapack', 'blas', 'armadillo']

hfb3 = Extension('_hfb3', ['hfb3.i'],
                 include_dirs = include_dirs,
                 libraries = libraries,
                 library_dirs = [
                                 get_env('ARMADILLO_LIBDIR'),
                                 get_env('BOOST_LIBDIR'),
                                 # get_env('LAPACK_LIBDIR'),
                                 # get_env('BLIS_AOCL_LIBDIR'),
                                ],
                 extra_compile_args = ['-std=c++14'] + get_skill_opts(),

                 extra_objects = [hfb3_build_dir + '/libhfb3.a', ],
                 # depends = [hfb3_build_dir + '/libhfb3.a', ],
                 swig_opts = ["-Wall", "-c++", "-I../../src", "-I./include/armanpy/", "-I$PYTHON3_ROOT/include"],
                )

setup(
    name = 'hfb3',
    version = '1.0.0',
    install_requires= ['numpy', 'msgpack'],
    ext_modules = [ hfb3 ],
    cmdclass = {'build_py' : build_py}, # <= fix for the bug
    py_modules = [ 'hfb3' ],
)
