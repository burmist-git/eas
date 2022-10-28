RELVERSION  = $(shell cat .release)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

MakefileFullPath = $(abspath $(lastword $(MAKEFILE_LIST)))
MakefileDirFullPath = $(shell dirname $(MakefileFullPath))
INSTALLDIR = $(MakefileDirFullPath)/install.$(RELVERSION)/

CXX  = g++
CXX += -I./

CXXFLAGS  = -g -Wall -fPIC -Wno-deprecated
CXXFLAGS += $(ROOTCFLAGS)
CXXFLAGS += $(ROOTLIBS)
CXXFLAGS += $(ROOTGLIBS)
CXXFLAGS += -std=c++14
CXXFLAGS += -fconcepts

#----------------------------------------------------#

all: eas

.PHONY: printmakehelp_and_reminder
printmakehelp_and_reminder: README Makefile
	$(info  /******************************************************************************/)
	$(info  * task --> printmakehelp_and_reminder: README Makefile          *)
	$(info  * $$@ ----> $@                                         *)
	$(info  * $$< --------------------------------> $<                   *)
	$(info  * $$^ --------------------------------> $^          *)
	$(info  /******************************************************************************/)

eas: eas.C
	$(CXX) -o $@ $^ $(CXXFLAGS)

clean:
	rm -f eas
	rm -f *~

