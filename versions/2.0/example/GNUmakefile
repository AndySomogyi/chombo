# -*- Mode: Makefile; Modified: "Fri 13 Jan 2006 13:15:50 by dbs"; -*- 

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

dirs = AMRPoisson/exec  PPMAMRGodunov/execPolytropic PPMAMRGodunov/execIdealMHD \
        PPMAMRGodunov/execResistiveMHD AMRUpwind/exec \
       AMRPoisson/designExec AMRPoisson/nodeExec expHeat/slickVersion

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
