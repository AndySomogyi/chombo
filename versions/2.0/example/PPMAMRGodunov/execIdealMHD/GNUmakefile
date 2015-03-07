# -*- Mode: Makefile; Modified: "Fri 19 Aug 2005 15:04:58 by dbs"; -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_PPMAMRGodunov_execIdealMHD

# the base name(s) of the application(s) in this directory
ebase = amrGodunov

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = AMRTimeDependent AMRElliptic AMRTools BoxTools 

# relative paths to source code directories
base_dir = .
src_dirs = ../srcIdealMHD

# input file for 'run' target
INPUT = wave.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# Note: tell Make.rules to not use -q for 'run' target
_space = # # NOTE: keep the space after the "=" or bad things happen
VERBOSE = $(_space)

# application-specific targets
doxygen:
	mkdir -p doc/doxygen
	sed s@../../Chombo/lib@$(CHOMBO_HOME)@g doxygen.config > /tmp/doxygen.config.tmp.$$$$ ; \
	$(DOXYGEN) /tmp/doxygen.config.tmp.$$$$ | \grep -v '^\(Preprocessing\|Parsing\|Generating code\|Generating docs\)' ; \
	$(RM) /tmp/doxygen.config.tmp.$$$$
	@echo "" ; echo "point browser at doc/doxygen/html/index.html" ; echo ""