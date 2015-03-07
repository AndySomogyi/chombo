# -*- Mode: Makefile -*- 

## This makefile builds one or more subdirectories.

## NOTE: requires USE_EB=TRUE

## Define the variables needed by Make.example.subdirs

# trace the chain of included makefiles
makefiles += example_EBAdvect

# relative path to the Chombo "lib" directory
CHOMBO_HOME = ../../lib

# relative paths to the subdirectories to build
exec_dirs = conv/advectiveDerivative conv/solution conv/simpleFilter conv/extrapScalar \
    conv/advectiveVel conv/udelsBoth conv/divFilteredUConv conv/gradDivFilterConv \
    conv/truncation conv/noflowtrunc conv/noflowsolution exec/viscInflowOutflow \
    conv/divUUConservation \
    exec/eulerInflowOutflow exec/eulerNoFlow

# application-specific variables go here ...

# the shared code for building in example subdirectories
include $(CHOMBO_HOME)/mk/Make.example.subdirs

# application-specific targets go here ...
