# -*- Mode: Makefile; Modified: "Tue 11 May 2004 09:58:43 by dbs"; -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_EllipticDesign

ebase := testAMRPoissonOp testBiCGStab testMultiGrid

LibNames := EllipticDesign AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test
