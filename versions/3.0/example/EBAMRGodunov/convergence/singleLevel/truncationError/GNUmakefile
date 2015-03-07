# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

## Note: requires USE_EB=TRUE

# trace the chain of included makefiles
makefiles += example_EBAMRGodunov_convergence_singleLevel_truncationError

# the base name(s) of the application(s) in this directory
ebase = convergenceTest

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../../../src

# input file for 'run' target
INPUT = trunc.smallrun.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

XTRACPPFLAGS += -DMODIANO_PROBLEM

# application-specific targets
