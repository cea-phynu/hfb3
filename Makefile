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

include Makefile.in

PWD="./"

TARGETLIB=lib/libhfb3.a
TARGETBIN=bin/hfb3

all: info $(TARGETBIN)

.PHONY: doc \
				clean \
				msgpack \
				gtest \
				semiclean \
				python_install \
				tests \
				$(TARGETBIN) \
				$(TARGETLIB)

$(TARGETLIB): info msgpack
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C src ../$(TARGETLIB)

$(TARGETBIN): info $(TARGETLIB)
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C src ../$(TARGETBIN)

python_install: $(TARGETLIB)
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C misc/python install

tests: $(TARGETLIB) gtest
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C misc/tests target

doc:
	@$(call generate,$(PWD),"documentation pages")
	@doxygen misc/doxygen/Doxyfile

msgpack:
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C misc/deps msgpack-c

gtest:
	@$(MAKE) $(MAKE_GLOBAL_FLAGS) -C misc/deps googletest

clean:
	@$(MAKE) clean $(MAKE_GLOBAL_FLAGS) -C src
	@$(MAKE) clean $(MAKE_GLOBAL_FLAGS) -C misc/python
	@$(MAKE) clean $(MAKE_GLOBAL_FLAGS) -C misc/tests
	@$(MAKE) clean $(MAKE_GLOBAL_FLAGS) -C misc/deps
	@rm -rf doc


VT = "\e[7;30;45m"
BT = "\e[7;30;44m"
GT = "\e[7;30;42m"

VG = "\e[0;35;42m"
GV = "\e[0;32;45m"

BG = "\e[0;34;42m"
GB = "\e[0;32;44m"

VB = "\e[0;34;45m"
BV = "\e[0;35;44m"

info:
ifeq ($(LOGO),ANSI8B)
	# @echo -e
	# @echo -e "\x1b[0;32;49m\ue0ba██\x1b[0;35;42m\ue0ba\x1b[0;32;49m  \ue0ba██\x1b[0;35;42m\ue0ba\x1b[0;32;49m\ue0ba████████\x1b[0;35;42m\ue0ba\x1b[0;32;49m\ue0ba█████\x1b[0;35;42m\ue0ba\x1b[0;32;49m   \x1b[0;34;49m \x1b[0;32;49m\ue0ba████████\x1b[0;35;42m\ue0ba\x1b[0m"
	# @echo -e "\x1b[0;34;49m███\x1b[0;35;49m█  \x1b[0;34;49m███\x1b[0;35;49m█\x1b[0;34;49m█████████\x1b[0;35;49m█\x1b[0;34;49m██████\x1b[0;35;42m\ue0bc\x1b[0;32;49m██\x1b[0;35;42m\ue0ba\x1b[0;34;49m █████████\x1b[0;35;49m█\x1b[0m"
	# @echo -e "\x1b[0;34;49m███\x1b[0;35;42m\ue0bc\x1b[0;32;49m██\x1b[0;34;49m███\x1b[0;35;49m█\x1b[0;34;49m███\x1b[0;32;44m▁▁▁\x1b[0;34;49m███\x1b[0;35;49m\ue0bc\x1b[0;34;49m███\x1b[0;32;44m▁▁▁\x1b[0;34;49m███\x1b[0;35;49m█\x1b[0;34;49m ███\x1b[0;32;44m▁▁▁\x1b[0;34;49m███\x1b[0;35;49m█\x1b[0m"
	# @echo -e "\x1b[0;34;49m█████████\x1b[0;35;49m█\x1b[0;34;49m██████\x1b[0;35;49m█   \x1b[0;34;49m█████████\x1b[0;35;49m█\x1b[0;32;49m    \x1b[0;34;49m██████\x1b[0;35;49m█\x1b[0m"
	# @echo -e "\x1b[0;34;49m█████████\x1b[0;35;49m█\x1b[0;34;49m██████\x1b[0;35;49m\ue0bc   \x1b[0;34;49m███\x1b[0;32;44m▁▁▁\x1b[0;34;49m███\x1b[0;35;49m█\x1b[0;34;49m \x1b[0;32;49m\ue0ba██\x1b[0;32;44m▁▁▁\x1b[0;34;49m███\x1b[0;35;49m█\x1b[0m"
	# @echo -e "\x1b[0;34;49m███\x1b[0;35;49m█  \x1b[0;34;49m███\x1b[0;35;49m█\x1b[0;34;49m███\x1b[0;35;49m█      \x1b[0;34;49m█████████\x1b[0;35;49m\ue0bc\x1b[0;34;49m █████████\x1b[0;35;49m█\x1b[0m"
	# @echo -e "\x1b[0;34;49m███\x1b[0;35;49m\ue0bc  \x1b[0;34;49m███\x1b[0;35;49m\ue0bc\x1b[0;34;49m███\x1b[0;35;49m\ue0bc      \x1b[0;34;49m██████\x1b[0;35;49m\ue0bc\x1b[0;32;49m   \x1b[0;34;49m █████████\x1b[0;35;49m\ue0bc\x1b[0m"
	# @echo -e
endif

ifeq ($(LOGO),ANSI24B)
	@echo -e
	@echo -e "\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m  \ue0ba██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba████████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m\ue0ba█████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;128;179;255m   \x1b[38;2;42;127;255m \x1b[38;2;229;128;255m\ue0ba████████\x1b[48;2;229;128;255m\x1b[38;2;136;0;170m\ue0ba\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0bc\x1b[0;49m\x1b[38;2;128;179;255m██\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0ba\x1b[0;49m\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m█\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[48;2;128;179;255m\x1b[38;2;0;68;170m\ue0bc\x1b[0;49m\x1b[38;2;128;179;255m██\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m███\x1b[48;2;212;42;255m\x1b[38;2;229;128;255m▁▁▁\x1b[0;49m\x1b[38;2;212;42;255m███\x1b[38;2;136;0;170m█\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m█   \x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;128;179;255m \x1b[38;2;229;128;255m   \x1b[38;2;212;42;255m██████\x1b[38;2;136;0;170m█\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m\ue0bc   \x1b[38;2;42;127;255m███\x1b[48;2;42;127;255m\x1b[38;2;128;179;255m▁▁▁\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m \x1b[38;2;229;128;255m\ue0ba██\x1b[48;2;212;42;255m▁▁▁\x1b[0;49m\x1b[38;2;212;42;255m███\x1b[38;2;136;0;170m█\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m█      \x1b[38;2;42;127;255m█████████\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m█\x1b[0m"
	@echo -e "\x1b[0;49m\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc  \x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;42;127;255m███\x1b[38;2;0;68;170m\ue0bc      \x1b[38;2;42;127;255m██████\x1b[38;2;0;68;170m\ue0bc\x1b[38;2;128;179;255m   \x1b[38;2;42;127;255m \x1b[38;2;212;42;255m█████████\x1b[38;2;136;0;170m\ue0bc\x1b[0m"
	@echo -e
endif

ifeq ($(LOGO),ASCII)
	@echo -e
	@echo "  +---+   +---+ +-----------+ +-------+       +-----------+"
	@echo " /   /|  /   /|/           /|/       /+---+  /           /|"
	@echo "+---+ | +---+ +-----------+ +-------+/   /| +-----------+ |"
	@echo "|   | +-|   | |           | |       +---+ | |           | |"
	@echo "|   |/  |   | |           |/|           | | |           | |"
	@echo "|   +---+   | |   ----+---+ +   -----   | | +---+----   | |"
	@echo "|           | |       | +   |           | |   +-|       | |"
	@echo "|           | |       |/    |           | |  /  |       | |"
	@echo "|   +---+   | |   +---+     |   -----   | + +---+----   | |"
	@echo "|   | + |   | |   | +       |           |/  |           | +"
	@echo "|   |/  |   |/|   |/        |       +---+   |           |/ "
	@echo "+---+   +---+ +---+         +-------+       +-----------+  "
	@echo -e
endif

	@echo -e $(GREEN)  "PLATFORM  "$(YELLOW)"?|"$(NORM) $(shell uname)
	@echo -e $(GREEN)  "ARCH      "$(YELLOW)"?|"$(NORM) $(shell uname -p)
	@echo -e $(GREEN)  "CXX       "$(YELLOW)"?|"$(NORM) $(CXX) $(shell $(CXX) -dumpfullversion)
	@echo -e $(GREEN)  "VERSION   "$(YELLOW)"?|"$(NORM) $(GIT_VERSION_FLAG)
	@echo -e $(GREEN)  "THREADS   "$(BLUE)":|"$(NORM) $(shell getconf _NPROCESSORS_ONLN)
	@echo -e $(GREEN)  "SKILL     "$(BLUE)":|"$(NORM) $(SKILLOPTS)
	@echo -e $(GREEN)  "REQFLAGS  "$(RED)"!|"$(NORM) $(REQFLAGS)
	@echo -e $(GREEN)  "WARNFLAGS "$(BLUE)":|"$(NORM) $(WARNFLAGS)
	@echo -e $(GREEN)  "LDFLAGS   "$(BLUE)":|"$(NORM) $(GLOBAL_LDFLAGS)
	@echo -e


