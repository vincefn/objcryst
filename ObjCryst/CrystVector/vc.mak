!include ..\rules.mak

libCrystVector.lib : CrystVector.obj
	lib -OUT:libCrystVector.lib CrystVector.obj

lib: libCrystVector.lib

default: lib

all: lib
