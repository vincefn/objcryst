# Base ObjCryst directory
DIR_CRYST = .

#external libraries directories
DIR_ATOMINFO = ${DIR_CRYST}/../atominfo
DIR_BLITZ = ${DIR_CRYST}/../blitz
DIR_NEWMAT = ${DIR_CRYST}/../newmat
DIR_SGLITE = ${DIR_CRYST}/../sglite
DIR_TAU = ${DIR_CRYST}/../tau
# Just to find wxWindows headers. Libs are in /usr/lib or in /usr/local/lib
DIR_WXWINDOWS = /usr/lib/wx

#Internal directories
DIR_CRYSTVECTOR = ${DIR_CRYST}/CrystVector
DIR_EXAMPLE = ${DIR_CRYST}/example
DIR_LIBCRYST = ${DIR_CRYST}/ObjCryst
DIR_REFOBJ = ${DIR_CRYST}/RefinableObj
DIR_VFNQUIRKS = ${DIR_CRYST}/Quirks
DIR_WXWCRYST = ${DIR_CRYST}/wxCryst
DIR_DOC := ${DIR_CRYST}/doc

### Rules for Linux & GCC
# C compiler
CC     := gcc
CFLAGS  = ${DEPENDFLAGS}
# C++ compiler
CXX      := c++
CXXFLAGS  = ${DEPENDFLAGS} ${PROFILEFLAGS}
# FORTRAN compiler
FC     := f77
FFLAGS  = 
# linker
LINKER    := c++
LDFLAGS   = -L/usr/lib -L/usr/local/lib -L$(DIR_ATOMINFO) -L$(DIR_CRYSTVECTOR) -L$(DIR_LIBCRYST) -L$(DIR_NEWMAT) -L$(DIR_BLITZ)/lib -L$(DIR_REFOBJ) -L$(DIR_SGLITE) -L$(DIR_VFNQUIRKS) -L$(DIR_WXWCRYST) -L$(DIR_TAU)/i386_linux/lib

#to automatically generate dependencies
MAKEDEPEND = gcc -MM ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} $< > $*.dep

# header files
SEARCHDIRS = -I- -I${DIR_CRYST}/.. -I./ -I$(DIR_BLITZ)  -I$(DIR_TAU)/include -I$(DIR_NEWMAT) -I${DIR_CRYST}

#wxWindows flags
ifeq ($(wxcryst),1)
   WXCRYSTFLAGS = -D__WX__CRYST__ `/usr/local/bin/wx-config --cxxflags`
   WX_LDFLAGS = -L/usr/X11R6/lib -lwxcryst `/usr/local/bin/wx-config --libs` $(GL_WX_LIB)
else
   WXCRYSTFLAGS :=
   WX_LDFLAGS :=
endif

#activate profiling using TAU package
ifeq ($(profile),1)
   PROFILEFLAGS := -DPROFILING_ON -DTAU_STDCXXLIB
   PROFILELIB := -ltau
else
   PROFILEFLAGS :=
   PROFILELIB :=
endif

#Using OpenGL ?
ifeq ($(opengl),1)
GL_WX_LIB = `/usr/local/bin/wx-config --gl-libs`
GL_FLAGS := -DOBJCRYST_GL -IGL
else
GL_WX_LIB :=
GL_FLAGS :=
endif

#Set DEBUG options
#for Blitz++: -ftemplate-depth-30 
ifeq ($(debug),1)
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS}
   CPPFLAGS = -g -Wall -D__DEBUG__ 
   LOADLIBES = -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lsglite -latominfo ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS}
else
# do not use -fomit-frame-pointer, or throw() catch() does not work !! GCC BUG ?
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS}
   CPPFLAGS = -O3 -w -ffast-math
   LOADLIBES = -s -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lsglite -latominfo ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS}
endif

# Add to statically link: -nodefaultlibs -lgcc /usr/lib/libstdc++.a


######################################################################
#####################      LIBRAIRIES         ########################
######################################################################

#LibCryst++
libCryst:
	$(MAKE) -f gnu.mak -C ${DIR_LIBCRYST} lib

libcryst: libCryst

#wxCryst++
ifeq ($(wxcryst),1)
libwxCryst:
	$(MAKE) -f gnu.mak -C ${DIR_WXWCRYST} lib
else
libwxCryst:
endif

#Vector computation library
libCrystVector:
	$(MAKE) -f gnu.mak -C ${DIR_CRYSTVECTOR} lib

#Quirks, including a (crude) library to display float, vectors, matrices, strings with some formatting..
libQuirks:
	$(MAKE) -f gnu.mak -C ${DIR_VFNQUIRKS} lib

#Library to take care of refinable parameters, plus Global optimization and Least Squares refinements
libRefinableObj:
	$(MAKE) -f gnu.mak -C ${DIR_REFOBJ} lib

#Newmat Matrix Algebra library (used for SVD)
libnewmat:
	$(MAKE) -f nm_gnu.mak -C ${DIR_NEWMAT} libnewmat.a
     
#SgLite -Spacegroup Lib
libsglite:
	$(MAKE) -f gnu.mak -C ${DIR_SGLITE} lib
	ranlib ${DIR_SGLITE}/libsglite.a
     
#AtomInfo
libatominfo:
	$(MAKE) -f gnu.mak -C ${DIR_ATOMINFO} lib
