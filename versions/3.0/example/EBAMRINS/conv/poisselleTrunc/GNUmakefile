# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

## Note: requires USE_EB=TRUE

# trace the chain of included makefiles
makefiles += example_EBAMRINS_conv_truncation

# the base name(s) of the application(s) in this directory
ebase = truncationError

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRTimeDependent EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../../src

# input file for 'run' target
INPUT = regression.inputs

include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
