# -*- Mode: Makefile -*- 

### This is the top-level GNUMakefile.  It requires GNU make. It builds some
### or all of the example directories, depending on the targets specified.

### The way this makefile works is the target specified by the user
### on the command line is treated as the `action' to perform on all
### the `example' subdirectories.  In the rules, the subdirectories are
### treated like targets and the action is specified by a variable.
### The makefile builds each directory target by invoking make
### recursively on the directory with the action as the recursive target.
### This assumes that the makefiles in the subdirectories implement all
### the actions in the $(_all_actions) variable.

export makefiles += example

## This variable specifies where we are.  Lower-level make files will need an
## absolute path so we use 'pwd' rather than "."
export CHOMBO_HOME := $(shell pwd)/../lib

## Get some configuration variables that might be needed
include $(CHOMBO_HOME)/mk/Make.defs.defaults

## This variable defines the names of the subdirectories to be built by default.
## One or more can be used as a target on the command line.
#[NOTE: GNUmakefile.release doesn't have some of these examples]

ifneq ($(USE_ALL),TRUE)
  dirs = OldAMRGodunov/Unsplit/execAdvectDiffuse OldAMRGodunov/Unsplit/execPolytropic \
      OldAMRGodunov/Unsplit/convergencePolytropic/truncation AMRNodeElliptic/execDirichlet \
      AMRNodeElliptic/execDirichletAlt AMRNodeElliptic/execOpExponential AMRNodeElliptic/execOpQuadratic \
      AMRNodeElliptic/execOpTrig AMRNodeElliptic/execRunPeriodic AMRNodeElliptic/execRunPoly \
      AMRPoisson/exec AMRWaveEqn/exec IO/AMR IO/singleLevel AMRGodunov/execPolytropic \
      MappedGridPPMAMRGodunov/execPolytropic  AMRPoisson/designExec AMRPoisson/nodeExec
  ifeq ($(USE_EB),TRUE)
    dirs += EBAMRGodunov/exec/planarShock EBAMRGodunov/exec/freeStream \
      EBAMRGodunov/convergence/singleLevel/anaDivTest EBAMRGodunov/convergence/singleLevel/autopsy \
      EBAMRGodunov/convergence/singleLevel/divergenceTrunc EBAMRGodunov/convergence/singleLevel/fluxInterp \
      EBAMRGodunov/convergence/singleLevel/solutionError EBAMRGodunov/convergence/singleLevel/truncationError \
      EBMenagerie/exec
  endif
  # old examples that are not part of the release
  ifeq ($(USE_OLD),TRUE)
    dirs += OldAMRGodunov/Split/exec RichardsonAMRPoisson/exec expHeat/parVersion expHeat/sgVersion expHeat/slickVersion
  endif

else
  # special case for USE_ALL=TRUE on the command line (takes a long time to build)
  export USE_MF = TRUE
  export USE_EB = TRUE
  dirs = OldAMRGodunov/Unsplit/execAdvectDiffuse OldAMRGodunov/Unsplit/convergencePolytropic/truncation \
    OldAMRGodunov/Unsplit/execPolytropic \
   AMRNodeElliptic/execDirichletAlt AMRNodeElliptic/execRunPoly AMRNodeElliptic/execDirichlet \
    AMRNodeElliptic/execOpQuadratic AMRNodeElliptic/execOpExponential AMRNodeElliptic/execRunPeriodic \
    AMRNodeElliptic/execOpTrig \
   AMRPoisson/exec \
   AMRSelfGravity/exec \
   AMRUpwind/exec \
   AMRWaveEqn/exec \
   EBAMRGodunov/convergence/singleLevel/anaDivTest EBAMRGodunov/convergence/singleLevel/coveredFaceExtrapTest \
    EBAMRGodunov/convergence/singleLevel/boundaryAutopsy \
    EBAMRGodunov/convergence/singleLevel/progressionTest EBAMRGodunov/convergence/singleLevel/autopsy \
    EBAMRGodunov/convergence/singleLevel/solutionError EBAMRGodunov/convergence/singleLevel/divergenceTrunc \
    EBAMRGodunov/convergence/singleLevel/truncationError EBAMRGodunov/convergence/singleLevel/fluxInterp \
    EBAMRGodunov/exec/freeStream EBAMRGodunov/exec/explosion EBAMRGodunov/exec/planarShock EBAMRGodunov/rzexec/explosion \
    EBAMRGodunov/rzexec/nozzle EBAMRGodunov/rzexec/planarShock EBAMRGodunov/convergence/singleLevel/divDiffTest \
   EBAMRINS/conv/advectiveDerivative EBAMRINS/conv/solution EBAMRINS/conv/simpleFilter EBAMRINS/conv/extrapScalar \
    EBAMRINS/conv/advectiveVel EBAMRINS/conv/udelsBoth EBAMRINS/conv/divFilteredUConv EBAMRINS/conv/gradDivFilterConv \
    EBAMRINS/conv/truncation EBAMRINS/conv/noflowtrunc EBAMRINS/conv/noflowsolution EBAMRINS/exec/viscInflowOutflow \
    EBAMRINS/exec/eulerInflowOutflow EBAMRINS/exec/eulerNoFlow \
   EBAMRNS/execns\
   EBAMRPoisson/exec EBAMRPoisson/allRegExec \
      EBAMRPoisson/tests/amr/heatEqnError EBAMRPoisson/tests/amr/solutionError \
      EBAMRPoisson/tests/amr/solutionVariableCoefError EBAMRPoisson/tests/amr/truncationError \
      EBAMRPoisson/tests/amr/truncationVariableCoefError \
      EBAMRPoisson/tests/projections/cellProject EBAMRPoisson/tests/projections/compositeCell \
      EBAMRPoisson/tests/projections/compositeDivU EBAMRPoisson/tests/projections/compositeGrad \
      EBAMRPoisson/tests/projections/compositeMAC EBAMRPoisson/tests/projections/compositeRhoCell \
      EBAMRPoisson/tests/projections/compositeRhoMAC EBAMRPoisson/tests/projections/macProject \
   EBAdvect/convergence/truncationError EBAdvect/convergence/advectiveDerivative \
    EBAdvect/convergence/solutionError EBAdvect/exec \
   EBMenagerie/exec \
   IO/AMR IO/singleLevel \
   AMRGodunov/execIdealMHD AMRGodunov/execPolytropic AMRGodunov/execResistiveMHD
  ifeq ($(DIM),2)
    dirs += AMRClaw/applications/acoustics/2d/varying/interface AMRClaw/applications/acoustics/2d/WaveEqn \
     AMRClaw/applications/advection/2d/constant AMRClaw/applications/advection/2d/swirl \
     AMRClaw/applications/euler/2d/cylinder
  else
    dirs += AMRClaw/applications/acoustics/3d/varying/interface \
     AMRClaw/applications/advection/3d/constant AMRClaw/applications/advection/3d/swirl \
     AMRClaw/applications/euler/3d/shockbubble AMRClaw/applications/euler/3d/simplewave
  endif
  # old examples that are not part of the release
  dirs += OldAMRGodunov/Split/exec expHeat/sgVersion expHeat/parVersion expHeat/slickVersion RichardsonAMRPoisson/exec
  # these are not working
  broken_dirs = OldAMRGodunov/Unsplit/execIdealMHD AMROpacity/exec EBAMRNS/exec #EBAMRNS/execns
endif

##
## Targets implemented in this makefile are:
##   all           build all programs and the libraries they depend on
##   example       <same as above>
## Targets implemented as actions in the recursive makefiles:
##   example-only  builds all the executables, but doesn't try to rebuild libs
##   run           runs all the executables
##   clean         deletes files for this configuration
##   realclean     deletes files for all configurations
##   vars          prints values of various make variables

_actions = example-only run run-only clean

.PHONY: all example $(_actions) $(dirs) usage vars

##
## Default target
##
usage:
	@echo "This makefile supports targets:"
	@echo "     $(_actions) $(dirs) realclean usage vars"

##
## Targets implmented here
##
all example  :
	+$(ECHO)$(MAKE) --no-print-directory --directory ../lib lib
	+$(ECHO)$(MAKE) --no-print-directory example-only

##
## Targets implemented in subdirectories
##
## For each action, set the action variable to tell the recursive
## makefile what to do.  There should be a separate dependency for
## every value in $(_actions).
##[NOTE: the `NODEPENDS=TRUE' variable disables automatic dependency
##       analysis in the recursive makefile, which is useless for
##       actions that delete files. <dbs>]
##
_action = example-only

example-only : _action = example-only

run          : _action = run       NODEPENDS=TRUE

run-only     : _action = run-only  NODEPENDS=TRUE

clean        : _action = clean     NODEPENDS=TRUE


# this dependency ensures the action is done for every directory.
$(_actions) : $(dirs)

# this rule performs the action
$(dirs):
	+$(MAKE) --no-print-directory --directory $@ $(_action)

# This target removes all non-essential files from all subdirectories for all configurations
realclean:
	-find . -type d \( -name d -o -name p -o -name f -o -name o \) -depth -exec $(RM) -r {} \;
	-find . \( -name \*.ex -o -name \*.a -o -name \*_F.H -o -name work.pc\* -o -name \*.ii -o -name \*.stb \) -exec $(RM) {} \;
	-find . \( -name gmon.out -o -name core\* -o -name \*.stackdump -o -name leak.out\* -o -name \*PROF.prof -o -name \*.tmp\? \) -exec $(RM) {} \;
	-find . \( -name pout.\[0-9-\]\* -o -name \*.hdf5 -o -name \*.h5 -o -name \*~ -o -name .\#\* \) -exec $(RM) {} \;

# This target prints out various 'make' variables to the screen for debugging
vars:
	@$(MAKE) --no-print-directory -f $(CHOMBO_HOME)/mk/Make.printVariables $@
