# -*- Mode: Makefile -*- 

## This makefile builds one or more subdirectories.

## NOTE: requires USE_EB=TRUE

## Define the variables needed by Make.example.subdirs

# trace the chain of included makefiles
makefiles += example_EBAMRPoisson

# relative path to the Chombo "lib" directory
CHOMBO_HOME = ../../lib

# relative paths to the subdirectories to build
exec_dirs = exec \
    tests/amr/solutionError \
    tests/amr/solutionVariableCoefError \
    tests/amr/truncationVariableCoefError \
    tests/amr/truncationError \
    tests/amr/heatEqnError \
    tests/projections/compositeDivU \
    tests/projections/compositeGrad \
    tests/projections/compositeMAC \
    tests/projections/cellProject \
    tests/projections/compositeRhoMAC \
    tests/projections/compositeCell \
    tests/projections/compositeCellWithFilter \
    tests/projections/macProject \
    tests/projections/compositeRhoCell \
    tests/projections/gradientAverage \
    tests/projections/paperProjTest


# application-specific variables go here ...

# the shared code for building in example subdirectories
include $(CHOMBO_HOME)/mk/Make.example.subdirs

# application-specific targets go here ...
