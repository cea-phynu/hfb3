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

include $(ROOT)/Makefile.in

TEST_TARGET = $(ROOT)/bin/tests
LIBHFB3 = $(ROOT)/lib/libhfb3.a

OBJ_DIR=.objs

TEST_MAIN_OBJECT=$(OBJ_DIR)/tests.o

TEST_OBJECTS=$(OBJ_DIR)/test_action.o                      \
             $(OBJ_DIR)/test_axis.o                        \
             $(OBJ_DIR)/test_base64.o                      \
             $(OBJ_DIR)/test_basis.o                       \
             $(OBJ_DIR)/test_state.o                       \
             $(OBJ_DIR)/test_cliparser.o                   \
             $(OBJ_DIR)/test_constraint.o                  \
             $(OBJ_DIR)/test_datatree.o                    \
             $(OBJ_DIR)/test_discrete.o                    \
             $(OBJ_DIR)/test_field_ws.o                    \
             $(OBJ_DIR)/test_field.o                       \
             $(OBJ_DIR)/test_fieldkinetic.o                \
             $(OBJ_DIR)/test_fieldspinorbit.o              \
             $(OBJ_DIR)/test_fragments.o                   \
             $(OBJ_DIR)/test_gradientwalk.o                \
             $(OBJ_DIR)/test_geometry.o                    \
             $(OBJ_DIR)/test_global.o                      \
             $(OBJ_DIR)/test_ioamedee.o                    \
             $(OBJ_DIR)/test_ioberger.o                    \
             $(OBJ_DIR)/test_iohfb3.o                      \
             $(OBJ_DIR)/test_iojson.o                      \
             $(OBJ_DIR)/test_iomsgp.o                      \
             $(OBJ_DIR)/test_mesh.o                        \
             $(OBJ_DIR)/test_meta.o                        \
             $(OBJ_DIR)/test_mixing.o                      \
             $(OBJ_DIR)/test_multipole_operators.o         \
             $(OBJ_DIR)/test_multi.o                       \
             $(OBJ_DIR)/test_plot.o                        \
             $(OBJ_DIR)/test_qnumbers.o                    \
             $(OBJ_DIR)/test_solver_ws.o                   \
             $(OBJ_DIR)/test_solverhfb.o                   \
             $(OBJ_DIR)/test_system.o                      \
             $(OBJ_DIR)/test_wspot.o                       \
             $(OBJ_DIR)/test_tools.o

COMMON_OPTS = $(GLOBAL_CXXFLAGS) ${ARMADILLO_CFLAGS} -I$(ROOT)/src -I$(ROOT)/$(GTEST_INCDIR)

COPTS = $(SKILLOPTS) $(COMMON_OPTS) -D$(SKILL)

LDFLAGS = $(LIBHFB3) ${ARMADILLO_LDFLAGS} $(ROOT)/$(GTEST_LIB) $(GLOBAL_LDFLAGS) ${BLIS_LDFLAGS}

.PHONY: clean data target
PWD="misc/tests/"

all:
	@$(call warning,$(PWD),"Do not launch this Makefile. Use the top Makefile instead.")

target: $(TEST_TARGET)

$(OBJ_DIR):
	@$(call directory,$(PWD),$@)
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: %.cpp
	@$(call compile,$(PWD),$<)
ifeq ($(QUIET),TRUE)
	@$(CXX) $(COPTS) -c $< -o $@
else
	$(CXX) $(COPTS) -c $< -o $@
endif

$(TEST_TARGET): $(TEST_MAIN_OBJECT) $(TEST_OBJECTS) $(LIBHFB3)
	@$(call link,$(PWD),$@,$(GTEST_MESSAGE))
ifeq ($(QUIET),TRUE)
	@$(CXX) $(COPTS) -o $@ $(TEST_MAIN_OBJECT) $(TEST_OBJECTS) $(LDFLAGS)
else
	$(CXX) $(COPTS) -o $@ $(TEST_MAIN_OBJECT) $(TEST_OBJECTS) $(LDFLAGS)
endif

$(TEST_OBJECTS): | $(OBJ_DIR)
$(TEST_MAIN_OBJECT): | $(OBJ_DIR)

clean:
	@$(call clean,$(PWD))
	@rm -f $(TEST_OBJECTS)
	@rm -f $(TEST_TARGET)
	@rm -f $(TEST_MAIN_OBJECT)
	@rm -rf $(OBJ_DIR)
	@rm -rf __pycache__/
