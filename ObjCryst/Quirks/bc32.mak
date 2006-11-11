!include ../rules.mak
   
libQuirks.lib : ci_string.obj VFNStreamFormat.obj VFNDebug.obj
	tlib "libQuirks.lib" -+ci_string.obj -+VFNStreamFormat -+VFNDebug.obj

lib: libQuirks.lib

default: lib

all: lib
