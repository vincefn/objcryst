!include ../rules.mak
   
libQuirks.lib : VFNStreamFormat.obj VFNDebug.obj
	lib -OUT:libQuirks.lib VFNStreamFormat.obj VFNDebug.obj

lib: libQuirks.lib

default: lib

all: lib
