# -*- Mode: Makefile -*- 

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# included makefiles need an absolute path
CHOMBO_HOME := $(shell pwd)/../..

makefiles+=lib_test_EBAMRTools

ebase = restart newIntRedistTest halfQuadTest allRegFluxRegTest   \
        aveConserveTest averageTest coarsenTest averageFluxTest   \
        pwlinterpTest pwQuadInterpTest fpExactTest fpqExactTest   \
        levelRedistTest fluxRegTest  fluxRegTestRZ fullRedistTest \
        quadCFITestEBCross

#ebase = newIntRedistTest fullRedistTest

LibNames = Workshop EBAMRTools EBTools AMRTools BoxTools

include $(CHOMBO_HOME)/mk/Make.test

CPPFLAGS += $(subst FALSE,,$(subst TRUE,-DCH_PROFILE,$(PROFILE)))
