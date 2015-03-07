# -*- Mode: Makefile -*- 

### This is the top-level GNUMakefile.  It requires GNU make. It builds some
### or all of the example directories, depending on the targets specified.

### User-definable variables are in UPPER case.
### Internal variables are in lower case.
###
### Usage:
###   target:    purpose:
###   -------    ----------
###   all        build all libraries and tests for this configuration
###   clean      delete all built files for this configuration
###   defs       print values of user-customizable variables
###   dist       build a tar file of the whole dir. tree (does a "realclean" first)
###   dist-noclean same as "dist" except doesn't do a "realclean" first
###   doc        build the documentation directories (config vars are ignored)
###   include    copies all the Chombo .H files to $(CHOMBO_HOME)/include
###   lib        build only the libraries for the specified configuration
###   realclean  delete everything that wasn't in the original distribution
###   run        build and run the test programs
###   run-tests  only runs the test programs, doesn't rebuild anything
###   setup      creates the Make.defs.local file from the template
###   tags       build TAGS file for emacs
###   test       build only the test programs (will also build necessary libraries)
###   usage      list the available targets and variables
###   vars       print variable values for all variables for debugging make files
###
### Notes:
###   Some of the above targets wont work in this directory.
###
### Notes:
###
##################################################################################

export makefiles+=example

## This variable specifies where we are.  Lower-level make files will need an
## absolute path so we use 'pwd' rather than "."
export CHOMBO_HOME := $(shell pwd)/..

## This variable defines the names of the subdirectories to be built by default.
## One or more can be used as a target on the command line.
dirs =           \
       trunc     \
       multigrid \

##
## These are the actions (targets for the recursive make) supported by this
## makefile.  It assumes the makefiles in the subdirectories support them too.
##   all        builds all the executables, including the libraries they depend on
##   example    same as 'all'
##   run        runs all the executables
##   clean      deletes files for this configuration
##   realclean  deletes files for all configurations
##
_all_actions = all example run clean realclean

.PHONY: $(_all_actions) $(dirs) usage

##
## Default target
##
usage:
	@echo "This makefile supports targets:"
	@echo "     $(_all_actions) $(dirs) usage vars"

##
## For each action, set the action variable to tell the recursive
## makefile what to do.  There should be a separate dependency for
## every value in $(_all_actions).
##[NOTE: the `NODEPENDS=TRUE' variable disables automatic
##       dependency analysis in the recursive makefile, which is only
##       needed for actions that compile or link code.]
##
_action = all
all       : _action = all
example   : _action = all
run       : _action = run       NODEPENDS=TRUE
clean     : _action = clean     NODEPENDS=TRUE

# this dependency ensures the action is done for every directory.
$(_all_actions) : $(dirs)

# this rule performs the action
$(dirs):
	$(MAKE) --no-print-directory --directory $@ $(_action)


realclean: _action = realclean	NODEPENDS=TRUE
