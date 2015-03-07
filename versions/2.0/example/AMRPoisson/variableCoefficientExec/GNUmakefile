# -*- Mode: Makefile; Modified: "Fri 19 Aug 2005 16:00:38 by dbs"; -*- 

## Define the variables needed by Make.example

# trace the chain of included makefiles
makefiles += example_VCAMRPoisson

# the base name(s) of the application(s) in this directory
ebase = VCPoissonSolve

# the location of Chombo lib dir
CHOMBO_HOME = ../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EllipticDesign AMRTools BoxTools

# relative paths to source code directories
base_dir = .

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
