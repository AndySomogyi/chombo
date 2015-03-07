# -*- Mode: Makefile -*- 

## This makefile builds one or more subdirectories.

## NOTE: requires USE_EB=TRUE

## Define the variables needed by Make.example.subdirs

# trace the chain of included makefiles
makefiles += example_EBAdvect

# relative path to the Chombo "lib" directory
CHOMBO_HOME = ../../lib

# relative paths to the subdirectories to build
exec_dirs = convergence/singleLevel/anaDivTest convergence/singleLevel/coveredFaceExtrapTest \
    convergence/singleLevel/boundaryAutopsy convergence/singleLevel/progressionTest \
    convergence/singleLevel/autopsy convergence/singleLevel/solutionError \
    convergence/singleLevel/divergenceTrunc convergence/singleLevel/truncationError \
    convergence/singleLevel/fluxInterp exec/freeStream exec/explosion exec/planarShock \
    rzexec/explosion rzexec/nozzle rzexec/planarShock convergence/singleLevel/divDiffTest

# application-specific variables go here ...

# the shared code for building in example subdirectories
include $(CHOMBO_HOME)/mk/Make.example.subdirs

# application-specific targets go here ...
