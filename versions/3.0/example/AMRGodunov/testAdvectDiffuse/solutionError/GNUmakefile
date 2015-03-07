# -*- Mode: Makefile -*- 

# the base name(s) of the application(s) in this directory
ebase = amrSolution

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../lib

# names of Chombo libraries needed by this program, in order of search.
LibNames =  AMRElliptic AMRTimeDependent AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../../srcAdvectDiffuse

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

