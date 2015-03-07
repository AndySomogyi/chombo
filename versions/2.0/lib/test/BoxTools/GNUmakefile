# -*- Mode: Makefile; Modified: "Tue 14 Mar 2006 16:38:13 by dbs"; -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_BoxTools

ebase =arrayTest broadcastTest copyTest  domainSplitTest\
  fabTest gatherTest HDF5attributes testCH_Attach HDF5boxIO HDF5data \
  newIVSTest testBox testRefCountedPtr testIntVectSet testParmParse \
  testLoadBalance testMeshRefine testPeriodic ivsfabTest testRealVect \
  test_complex test_parstream 

LibNames = BoxTools

include $(CHOMBO_HOME)/mk/Make.test
