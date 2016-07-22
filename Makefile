SHELL = /bin/bash
COMPILER = g++ -g
CFLAGS = $(shell root-config --cflags)
LIBS   = $(shell root-config --glibs --libs) -lPhysics -lThread -lMinuit -lHtml -lVMC -lEG -lGeom

.SUFFIXES: .cc .o

all :  Inclusive-Jet FitResolutionPtvsY FitResolutionPtvsPt

Inclusive-Jet : Inclusive-Jet.cc       						     
	${COMPILER} Inclusive-Jet.cc -o $@ ${CFLAGS} ${LIBS}; rm -rf Inclusive-Jet.dSYM

FitResolutionPtvsY : FitResolutionPtvsY.cc       						       
	${COMPILER} FitResolutionPtvsY.cc -o $@ ${CFLAGS} ${LIBS}; rm -rf FitResolutionPtvsY.dSYM

FitResolutionPtvsPt : FitResolutionPtvsPt.cc       						 
	${COMPILER} FitResolutionPtvsPt.cc -o $@ ${CFLAGS} ${LIBS}; rm -rf FitResolutionPtvsPt.dSYM

clean :
	rm Inclusive-Jet FitResolutionPtvsY FitResolutionPtvsPt

