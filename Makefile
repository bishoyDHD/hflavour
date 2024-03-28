SHELL=/usr/bin/env bash

#include config
-include Makefile.ini

PYTHIA=$(PREFIX_LIB)/libpythia8$(LIB_SUFFIX)

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS += -O -g -Wall -fPIC $(ROOTCFLAGS) -I$(PYTHIA8_INCLUDE_DIR) -I$(ROOTSYS)/include
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

EXEC = eicLFgen
SRCS = flavorGen.cc
DEPS = disvars.h storage.h
OBJS = $(subst .cc,.o,$(SRCS))

LDLIBS =$(ROOTLIBS) -L$(PYTHIA8_LIBRARY) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl


all: $(EXEC)

$(EXEC): $(OBJS) $(DEPS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(OBJS) $(LDLIBS)

#$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: clean

clean:
	-rm $(OBJS) $(EXEC)

