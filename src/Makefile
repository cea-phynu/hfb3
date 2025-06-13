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

include $(ROOT)/Makefile.in

LIB_TARGET=$(ROOT)/lib/libhfb3.a
MAIN_TARGET=$(ROOT)/bin/hfb3

OBJ_DIR=.objs

MAIN_OBJECT=$(OBJ_DIR)/main.o

LIB_OBJECTS=$(OBJ_DIR)/action.o                       \
            $(OBJ_DIR)/axis.o                         \
            $(OBJ_DIR)/base64.o                       \
            $(OBJ_DIR)/basis.o                        \
            $(OBJ_DIR)/state.o                        \
            $(OBJ_DIR)/cliparser.o                    \
            $(OBJ_DIR)/constraint.o                   \
            $(OBJ_DIR)/datatree.o                     \
            $(OBJ_DIR)/discrete.o                     \
            $(OBJ_DIR)/field.o                        \
            $(OBJ_DIR)/field_kinetic.o                \
            $(OBJ_DIR)/field_cdm2.o                   \
            $(OBJ_DIR)/field_central.o                \
            $(OBJ_DIR)/field_coulomb_slater.o         \
            $(OBJ_DIR)/field_coulomb_exact.o          \
            $(OBJ_DIR)/field_density.o                \
            $(OBJ_DIR)/field_density_fr.o             \
            $(OBJ_DIR)/field_rearrangement.o          \
            $(OBJ_DIR)/field_rearrangement_fr.o       \
            $(OBJ_DIR)/field_spin_orbit.o             \
            $(OBJ_DIR)/field_ws.o                     \
            $(OBJ_DIR)/fragments.o                    \
            $(OBJ_DIR)/general.o                      \
            $(OBJ_DIR)/generic.o                      \
            $(OBJ_DIR)/geometrical_operators.o        \
            $(OBJ_DIR)/geometry.o                     \
            $(OBJ_DIR)/global.o                       \
            $(OBJ_DIR)/gradientwalk.o                 \
            $(OBJ_DIR)/gzstream.o                     \
            $(OBJ_DIR)/interaction.o                  \
            $(OBJ_DIR)/io_amedee.o                    \
            $(OBJ_DIR)/io_berger.o                    \
            $(OBJ_DIR)/io_hfb3.o                      \
            $(OBJ_DIR)/io_json.o                      \
            $(OBJ_DIR)/io_msgp.o                      \
            $(OBJ_DIR)/logger.o                       \
            $(OBJ_DIR)/md5.o                          \
            $(OBJ_DIR)/mesh.o                         \
            $(OBJ_DIR)/mixing.o                       \
            $(OBJ_DIR)/multipole_operators.o          \
            $(OBJ_DIR)/multi.o                        \
            $(OBJ_DIR)/plot.o                         \
            $(OBJ_DIR)/qnumbers.o                     \
            $(OBJ_DIR)/quadratures.o                  \
            $(OBJ_DIR)/solver.o                       \
            $(OBJ_DIR)/solver_basis.o                 \
            $(OBJ_DIR)/solver_hfb_broyden.o           \
            $(OBJ_DIR)/solver_hfb_gradient.o          \
            $(OBJ_DIR)/solver_ws.o                    \
            $(OBJ_DIR)/system.o                       \
            $(OBJ_DIR)/states.o                       \
            $(OBJ_DIR)/wspot.o                        \
            $(OBJ_DIR)/tools.o                        \

DEP = $(LIB_OBJECTS:%.o=%.d) $(MAIN_OBJECT:%.o=%.d)


COMMON_OPTS = $(GLOBAL_CXXFLAGS) ${ARMADILLO_CFLAGS} -I$(ROOT)/$(MSGPACK_INCDIR)

COPTS = $(SKILLOPTS) $(COMMON_OPTS) -D$(SKILL) $(GIT_VERSION_FLAG)

LDFLAGS = ${ARMADILLO_LDFLAGS} $(GLOBAL_LDFLAGS) ${BLIS_LDFLAGS}
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
