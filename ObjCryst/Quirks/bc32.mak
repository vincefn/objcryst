!include ../rules.mak
   
libQuirks.lib : VFNStreamFormat.obj VFNDebug.obj
	tlib "libQuirks.lib" +VFNStreamFormat +VFNDebug.obj

lib: libQuirks.lib

default: lib

all: lib