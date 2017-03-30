# Base ObjCryst directory
DIR_CRYST = .

#Other directories
DIR_ATOMINFO = ${DIR_CRYST}/../atominfo
DIR_BLITZ = ${DIR_CRYST}/../blitz
DIR_CRYSTVECTOR = ${DIR_CRYST}/CrystVector
DIR_EXAMPLE = ${DIR_CRYST}/example
DIR_LIBCRYST = ${DIR_CRYST}/ObjCryst
DIR_NEWMAT = ${DIR_CRYST}/../newmat
DIR_REFOBJ = ${DIR_CRYST}/RefinableObj
DIR_SGINFO = ${DIR_CRYST}/../sginfo
DIR_SGLITE = ${DIR_CRYST}/../sglite
DIR_TAU = ${DIR_CRYST}/../tau
DIR_VFNQUIRKS = ${DIR_CRYST}/Quirks
DIR_WXWCRYST = ${DIR_CRYST}/wxCryst
DIR_WXWINDOWS = ${DIR_CRYST}/../wxWindows

DIR_DOC := ${DIR_CRYST}/doc

### Rules for Linux & GCC
# C compiler
CC     := gcc
CFLAGS  = ${DEPENDFLAGS}
# C++ compiler
CXX      := g++
CXXFLAGS  = ${DEPENDFLAGS} ${PROFILEFLAGS}
# FORTRAN compiler
FC     := f77
FFLAGS  =
# linker
LINKER    := g++
LDFLAGS   = -L$(DIR_ATOMINFO) -L$(DIR_CRYSTVECTOR) -L$(DIR_LIBCRYST) -L$(DIR_NEWMAT) -L$(DIR_BLITZ)/lib -L$(DIR_REFOBJ) -L$(DIR_SGLITE) -L$(DIR_VFNQUIRKS) -L$(DIR_TAU)/i386_linux/lib -L$(DIR_WXWCRYST)

#to automatically generate dependencies
MAKEDEPEND = gcc -MM ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} -o $*.dep $<

# header files
SEARCHDIRS = -I- -I${DIR_CRYST}/.. -I./ -I$(DIR_BLITZ)  -I$(DIR_TAU)/include -I$(DIR_NEWMAT) -I$(DIR_WXWINDOWS)/include/ -I${DIR_CRYST}

#wxWindows flags
ifeq ($(wxcryst),1)
   WXCRYSTFLAGS = -D__WX__CRYST__ --pipe -fvtable-thunks -c -D_X86_=1 -D_WIN32 -DWINVER=0x0400 -D__WIN95__ -D__GNUWIN32__ -DSTRICT  -D__WXMSW__ -D__WINDOWS__
   WX_LDFLAGS =  -L$(DIR_WXWINDOWS)/lib -L/usr/lib -L/usr/X11R6/lib -lwxcryst -lwx -lpng -ljpeg -lzlib -lxpm -ltiff -lstdc++ -lgcc -lcomctl32 -lole32 -luuid -lopengl32 -lglu32 -lglut32 -lgdi32 -lcomdlg32
	#-mwindows -Wl,--subsystem,windows -mno-cygwin -lodbc32 -lctl3d32 -loleaut32 -lwsock32 -lwinspool -lwinmm -lshell32 -ladvapi32
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
GL_DIR   = ${DIR_CRYST}/../OpenGL
GL_FLAGS :=
GL_LIB   :=
GL_WX_LIB =
GL_FLAGS := -DOBJCRYST_GL
GL_OBJ   :=
else
GL_DIR   :=
GL_FLAGS :=
GL_LIB   :=
GL_WX_LIB :=
GL_FLAGS :=
GL_OBJ   :=
endif

#Set DEBUG options
#for Blitz++: -ftemplate-depth-30
ifeq ($(debug),1)
   DEPENDFLAGS = -g -Wall -D__DEBUG__ ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS}
   BIN_DIR  := obj_debug
else
# do not use -fomit-frame-pointer, or throw() catch() does not work !! GCC BUG ?
   DEPENDFLAGS = -O2 -w -ffast-math ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS}
   BIN_DIR  := obj_optimized
endif
LOADLIBES = -s -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lsglite -latominfo ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS}

CPPFLAGS =

######################################################################
#####################      LIBRAIRIES         ########################
######################################################################

#LibCryst++
libCryst:
	$(MAKE) -f linux.mak -C ${DIR_LIBCRYST} lib

libcryst: libCryst

#wxCryst++
ifeq ($(wxcryst),1)
libwxCryst:
	$(MAKE) -f linux.mak -C ${DIR_WXWCRYST} lib
else
libwxCryst:
endif

#Vector computation library
libCrystVector:
	$(MAKE) -f linux.mak -C ${DIR_CRYSTVECTOR} lib

#Quirks, including a (crude) library to display float, vectors, matrices, strings with some formatting..
libQuirks:
	$(MAKE) -f linux.mak -C ${DIR_VFNQUIRKS} lib

#Library to take care of refinable parameters, plus Global optimization and Least Squares refinements
libRefinableObj:
	$(MAKE) -f linux.mak -C ${DIR_REFOBJ} lib

#Newmat Matrix Algebra library (used for SVD)
libnewmat:
	$(MAKE) -f nm_gnu.mak -C ${DIR_NEWMAT} libnewmat.a

#SgLite -Spacegroup Lib
libsglite:
	$(MAKE) -f linux.mak -C ${DIR_SGLITE} lib

#AtomInfo
libatominfo:
	$(MAKE) -f linux.mak -C ${DIR_ATOMINFO} lib
