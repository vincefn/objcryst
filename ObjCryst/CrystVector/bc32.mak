!include ..\rules.mak

libCrystVector.lib : CrystVector.obj
	tlib "libCrystVector.lib" -+CrystVector

lib: libCrystVector.lib

default: lib

all: lib
