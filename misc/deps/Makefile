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

include ../../Makefile.in

.PHONY: clean googletest
PWD="misc/deps/"

all:
	@$(call warning,"Do not launch this Makefile. Use the top Makefile instead.")

googletest: googletest.tgz
	@$(call generate,$(PWD),$@)
	@tar xfz $<
	@touch $@
	@mkdir -p $@/build
ifeq ($(QUIET),TRUE)
	@cd $@/build && cmake .. &> /dev/null && $(MAKE) &> /dev/null
else
	@cd $@/build && cmake .. && $(MAKE)
endif

msgpack-c: msgpack-c.tgz
	@$(call generate,$(PWD),$@)
	@tar xfz $<
	@touch $@

clean:
	@$(call clean,$(PWD))
	@rm -rf googletest
	@rm -rf msgpack-c
