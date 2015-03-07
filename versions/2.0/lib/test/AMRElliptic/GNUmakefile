# -*- Mode: Makefile; Modified: "Thu 24 Feb 2005 12:58:27 by dbs"; -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_AMRElliptic

ebase := testPoisson

LibNames := EllipticDesign AMRElliptic AMRTools  BoxTools

INPUT = inputs.testPoisson.$(DIM)d

include $(CHOMBO_HOME)/mk/Make.test
