!include ..\rules.mak

libCrystVector.lib : CrystVector.obj
	del libCrystVector.lib
	tlib "libCrystVector.lib" +CrystVector

lib: libCrystVector.lib

default: lib

all: lib