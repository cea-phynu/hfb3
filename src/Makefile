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

ROOT = ..
PWD = "src/"

VERSION=$(shell cat $(ROOT)/VERSION)
include $(ROOT)/Makefile.in

LIB_TARGET=$(ROOT)/lib/libhfb3.a
MAIN_TARGET=$(ROOT)/bin/hfb3

OBJ_DIR=.objs

MAIN_OBJECT=$(OBJ_DIR)/main.o

LIB_OBJECTS=
LIB_OBJECTS+=$(OBJ_DIR)/action.o
LIB_OBJECTS+=$(OBJ_DIR)/solver_alternator.o
LIB_OBJECTS+=$(OBJ_DIR)/axis.o
LIB_OBJECTS+=$(OBJ_DIR)/base64.o
LIB_OBJECTS+=$(OBJ_DIR)/basis.o
LIB_OBJECTS+=$(OBJ_DIR)/state.o
LIB_OBJECTS+=$(OBJ_DIR)/cliparser.o
LIB_OBJECTS+=$(OBJ_DIR)/constraint.o
LIB_OBJECTS+=$(OBJ_DIR)/datatree.o
LIB_OBJECTS+=$(OBJ_DIR)/discrete.o
LIB_OBJECTS+=$(OBJ_DIR)/field.o
LIB_OBJECTS+=$(OBJ_DIR)/field_kinetic.o
LIB_OBJECTS+=$(OBJ_DIR)/field_cdm2.o
LIB_OBJECTS+=$(OBJ_DIR)/field_central.o
LIB_OBJECTS+=$(OBJ_DIR)/field_coulomb_slater.o
LIB_OBJECTS+=$(OBJ_DIR)/field_coulomb_exact.o
LIB_OBJECTS+=$(OBJ_DIR)/field_density.o
LIB_OBJECTS+=$(OBJ_DIR)/field_density_fr.o
LIB_OBJECTS+=$(OBJ_DIR)/field_rearrangement.o
LIB_OBJECTS+=$(OBJ_DIR)/field_rearrangement_fr.o
LIB_OBJECTS+=$(OBJ_DIR)/field_spin_orbit.o
LIB_OBJECTS+=$(OBJ_DIR)/field_ws.o
LIB_OBJECTS+=$(OBJ_DIR)/fragments.o
LIB_OBJECTS+=$(OBJ_DIR)/general.o
LIB_OBJECTS+=$(OBJ_DIR)/generic.o
LIB_OBJECTS+=$(OBJ_DIR)/geometrical_operators.o
LIB_OBJECTS+=$(OBJ_DIR)/geometry.o
LIB_OBJECTS+=$(OBJ_DIR)/global.o
LIB_OBJECTS+=$(OBJ_DIR)/gradientwalk.o
LIB_OBJECTS+=$(OBJ_DIR)/gzstream.o
LIB_OBJECTS+=$(OBJ_DIR)/interaction.o
LIB_OBJECTS+=$(OBJ_DIR)/io_amedee.o
LIB_OBJECTS+=$(OBJ_DIR)/io_berger.o
LIB_OBJECTS+=$(OBJ_DIR)/io_hfb3.o
LIB_OBJECTS+=$(OBJ_DIR)/io_json.o
LIB_OBJECTS+=$(OBJ_DIR)/io_msgp.o
LIB_OBJECTS+=$(OBJ_DIR)/logger.o
LIB_OBJECTS+=$(OBJ_DIR)/md5.o
LIB_OBJECTS+=$(OBJ_DIR)/mesh.o
LIB_OBJECTS+=$(OBJ_DIR)/mixing.o
LIB_OBJECTS+=$(OBJ_DIR)/multipole_operators.o
LIB_OBJECTS+=$(OBJ_DIR)/multi.o
LIB_OBJECTS+=$(OBJ_DIR)/plot.o
LIB_OBJECTS+=$(OBJ_DIR)/qnumbers.o
LIB_OBJECTS+=$(OBJ_DIR)/quadratures.o
LIB_OBJECTS+=$(OBJ_DIR)/solver.o
LIB_OBJECTS+=$(OBJ_DIR)/solver_basis.o
LIB_OBJECTS+=$(OBJ_DIR)/solver_hfb_broyden.o
LIB_OBJECTS+=$(OBJ_DIR)/solver_hfb_gradient.o
LIB_OBJECTS+=$(OBJ_DIR)/solver_ws.o
LIB_OBJECTS+=$(OBJ_DIR)/system.o
LIB_OBJECTS+=$(OBJ_DIR)/states.o
LIB_OBJECTS+=$(OBJ_DIR)/wspot.o
LIB_OBJECTS+=$(OBJ_DIR)/tools.o

DEP = $(LIB_OBJECTS:%.o=%.d) $(MAIN_OBJECT:%.o=%.d)


COMMON_OPTS = $(GLOBAL_CXXFLAGS) ${ARMADILLO_CFLAGS} -I$(ROOT)/$(MSGPACK_INCDIR)

COPTS = $(SKILLOPTS) $(COMMON_OPTS) -D$(SKILL) $(GIT_VERSION_FLAG) $(VERSION_FLAG)

LDFLAGS = ${LAPACK_LDFLAGS} $(GLOBAL_LDFLAGS) ${BLIS_LDFLAGS}
AR=ar

.PHONY: clean

all:
	@$(call warning,"Do not launch this Makefile. Use the top Makefile instead.")

-include $(DEP)

$(OBJ_DIR):
	@$(call directory,$(PWD),$@)
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp
	@$(call compile,$(PWD),$<)
ifeq ($(QUIET),TRUE)
	@$(CXX) $(COPTS) -MMD -c $< -o $@
else
	$(CXX) $(COPTS) -MMD -c $< -o $@
endif

$(LIB_OBJECTS): | $(OBJ_DIR)

$(LIB_TARGET): $(LIB_OBJECTS)
	@$(call generate,$(PWD),$@)
	@$(AR) rcs $@ $(LIB_OBJECTS)

$(MAIN_TARGET): $(LIB_TARGET) $(MAIN_OBJECT)
	@$(call link,$(PWD),$@,"")
ifeq ($(QUIET),TRUE)
	@$(CXX) $(COPTS) -o $@ $(MAIN_OBJECT) $(LIB_TARGET) $(LDFLAGS)
else
	$(CXX) $(COPTS) -o $@ $(MAIN_OBJECT) $(LIB_TARGET) $(LDFLAGS)
endif

clean:
	@$(call clean,$(PWD))
	@rm -f $(MAIN_OBJECT)
	@rm -f $(LIB_OBJECTS)
	@rm -f $(LIB_TARGET)
	@rm -f $(MAIN_TARGET)
	@rm -rf $(OBJ_DIR)
	@rm -f tags
