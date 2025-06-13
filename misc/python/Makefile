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

ROOT = ../..
PWD = "misc/python/"

include ../../Makefile.in

LIBHFB3 = $(ROOT)/lib/libhfb3.a

# the following flag may be needed for old versions of pip
# PIP_OPTS = --use-feature=in-tree-build

.PHONY: clean install

all:
	@$(call warning,$(PWD),"Do not launch this Makefile. Use the top Makefile instead.")

$(LIBHFB3):
	@$(MAKE) -C $(ROOT)

install:
	@$(call install,$(PWD),"Python3 bindings")
ifeq ($(QUIET),TRUE)
	@python3 -m pip $(PIP_OPTS) install . > /dev/null
else
	@python3 -m pip $(PIP_OPTS) install .
endif

clean:
	@$(call clean,$(PWD))
	@rm -f hfb3.py
	@rm -f hfb3_wrap.cpp
	@rm -f hfb3_wrap.h
	@rm -rf build
	@rm -rf dist
	@rm -rf .tox
	@rm -rf test_venv
	@rm -rf __pycache__
	@rm -rf .cache
	@rm -rf .pycache
	@rm -rf hfb3.egg-info
	@rm -rf hfb3*whl
	@rm -rf _hfb3*so
	@rm -rf .eggs/
	@find . -name '*.pyc' -delete
