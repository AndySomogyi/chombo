# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_AMRNodeElliptic_OpQuadratic

# the base name(s) of the application(s) in this directory
ebase := testOpQuadratic

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames := AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src

# input file for 'run' target
INPUT = inputs.3level

include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
