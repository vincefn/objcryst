# Base ObjCryst directory
DIR_CRYST=C:\Dev\Fox\ObjCryst
# Borland Directory
DIR_BORLAND=C:\Borland\BCC55
# wxWindows Directory
DIR_WXWINDOWS = C:\Dev\wxWidgets

#You should not need to modify anything below
###################################################################
#use wxCryst ?
wxcryst=1
#debugging ? (only activates debugging messages)
debug=0

#Other directories
DIR_CCTBX = ${DIR_CRYST}\..\cctbx
DIR_CRYSTVECTOR = ${DIR_CRYST}\CrystVector
DIR_EXAMPLE = ${DIR_CRYST}\example
DIR_LIBCRYST = ${DIR_CRYST}\ObjCryst
DIR_NEWMAT = ${DIR_CRYST}\..\newmat
DIR_REFOBJ = ${DIR_CRYST}\RefinableObj
DIR_TAU = ${DIR_CRYST}\..\tau
DIR_VFNQUIRKS = ${DIR_CRYST}\Quirks
DIR_WXWCRYST = ${DIR_CRYST}\wxCryst

DIR_DOC := ${DIR_CRYST}/doc

### Rules for Borland C++ 5.5
###################################################################
#
#       Borland specific directives ---
#
.SWAP
.AUTODEPEND

###################################################################
#
#       set default values:

!ifndef ENVIRON
ENVIRON = WIN32
!endif

!ifndef BINDING
BINDING = STATIC
!endif

!ifndef THREAD
THREAD = MULTI
!endif

!ifndef BMODE
BMODE = RELEASE
!endif

###################################################################
#
# Flag illegal options:
#

!if $(ENVIRON) != WIN32
! error Illegal value for ENVIRON option
!endif

!if $(BINDING) != DLL && $(BINDING) != STATIC
!  error Illegal value for BINDING option
!endif

!if $(THREAD) != SINGLE && $(THREAD) != MULTI
!  error Illegal value for THREAD option
!endif

!if $(BMODE) != RELEASE && $(BMODE) != DEBUG
!  error Illegal value for BMODE option
!endif

###################################################################
#
# Set tool and version names:

!if $(ENVIRON) == WIN32
CPP        = bcc32
CPP32      = cpp32
LIBRARIAN  = tlib /P128
LINKER     = ilink32
RC         = brc32
ENVNAME    =
!endif

###################################################################
#
# Set the various flags:

!if $(BMODE) == DEBUG
DBGOPT= -v -N -x -xp
CCLINKOPT = -lGn
!else
CCLINKOPT = -lGn
!endif

!if $(THREAD) == MULTI
CCLINKOPT = $(CCLINKOPT) -tWM
LIBSUF=mt
!else
CCLINKOPT = $(CCLINKOPT) -tWM-
LIBSUF=
!endif

###################################################################
#
# Set any relevant defines (-Dxxx)

DEFOPTS =

!if $(BINDING) == DLL
DEFOPTS=$(DEFOPTS) -tWCR
TARGSUF=R
LIBSUF=$(LIBSUF)i
!else
DEFOPTS = $(DEFOPTS) -tWC
LIBSUF=$(LIBSUF)
TARGSUF=
!endif



###################################################################
#
# compiler options

PCHROOT=stl_pch
CCOPTS = -w- -jb -j3 -O2 -6 -ff -OS
#-Hc -H=$(PCHROOT).csm

#Compile flags:
CPPFLAGS= $(CCOPTS) $(DBGOPT)  $(ENVOPTS) $(DEFOPTS) $(THROPTS) $(CCLINKOPT) $(SEARCHDIRS) -DOBJCRYST_GL $(WXCRYSTFLAGS) $(CPPDEBUGFLAGS)
#LINKFLAGS= -Gn -Gi -Tpd -aa -L$(MAKEDIR)\..\lib -x
LINKFLAGS= ${WX_LDFLAGS} -L$(DIR_CCTBX) -L$(DIR_CRYSTVECTOR) -L$(DIR_LIBCRYST) -L$(DIR_NEWMAT) -L$(DIR_BLITZ)lib -L$(DIR_REFOBJ) -L$(DIR_VFNQUIRKS) -L$(DIR_WXWCRYST) -L$(DIR_WXWINDOWS)\lib \
	-L$(DIR_BORLAND)\lib;$(DIR_BORLAND)\lib\psdk

LINKSTARTUP= c0d32.obj
LINKLIBS= import32.lib cw32$(LIBSUF).lib libcctbx.lib libcrystvector.lib libquirks.lib librefinableobj.lib libcryst.lib $(WX_LIBS)
RCFLAGS= -r -i$(MAKEDIR)\..\include;$(MAKEDIR)\..\include\windows


# header files
SEARCHDIRS = -I${DIR_CRYST} -I. -I.. -I..\.. -I..\..\.. -I$(DIR_NEWMAT) -I${DIR_CRYST} -I$(DIR_CCTBX) -I$(DIR_BORLAND)\include -I$(DIR_WXWINDOWS)\include

#debug ?
!if $(debug)==1
   CPPDEBUGFLAGS = -D__DEBUG__
!else
   CPPDEBUGFLAGS :=
!endif

#wxWindows flags
!if $(wxcryst)==1
   WXCRYSTFLAGS = -D__WX__CRYST__ @$(DIR_WXWINDOWS)\src\msw\wxw32.cfg
#-aa : windows application -ap: console application
   WX_LDFLAGS = -ap
   WX_LIBS = libwxcryst.lib wx24s_bcc.lib jpeg.lib tiff.lib winpng.lib zlib.lib
!else
   WXCRYSTFLAGS :=
   WX_LDFLAGS :=
   WX_LIBS :=
!endif

#activate profiling using TAU package
!if $(profile)==1
   PROFILEFLAGS := -DPROFILING_ON -DTAU_STDCXXLIB
   PROFILELIB := -ltau
!else
   PROFILEFLAGS :=
   PROFILELIB :=
!endif

#Using OpenGL ?
!if $(opengl)==1
   GL_DIR   = ${DIR_CRYST}/../OpenGL
   GL_FLAGS := -I$(GL_DIR)/include
   GL_LIB   := -L$(GL_DIR)/lib
   GL_WX_LIB =
   GL_FLAGS := -DOBJCRYST_GL
   GL_OBJ   :=
!else
   GL_DIR   :=
   GL_FLAGS :=
   GL_LIB   :=
   GL_WX_LIB :=
   GL_FLAGS :=
   GL_OBJ   :=
!endif

############## TARGETS #####
#clean up
clean:
	del *.obj
	del *.lib
	del *.res

# build the ressource part of an application
.rc.res:
	brc32 -r -I$(DIR_BORLAND)\include -I$(DIR_WXWINDOWS)\include {$< }


######################################################################
#####################      LIBRAIRIES         ########################
######################################################################

#LibCryst++
libcryst:
	cd ${DIR_LIBCRYST}
	$(MAKE) -f bc32.mak lib

#wxCryst++
!if $(wxcryst)==1
libwxCryst:
	cd ${DIR_WXWCRYST}
	$(MAKE) -f bc32.mak  lib
!else
libwxCryst:
!endif

#Vector computation library
libCrystVector:
	cd ${DIR_CRYSTVECTOR}
	$(MAKE) -f bc32.mak lib

#Quirks, including a (crude) library to display float, vectors, matrices, strings with some formatting..
libQuirks:
	cd ${DIR_VFNQUIRKS}
	$(MAKE) -f bc32.mak lib

#Library to take care of refinable parameters, plus Global optimization and Least Squares refinements
libRefinableObj:
	cd ${DIR_REFOBJ}
	$(MAKE) -f bc32.mak lib

#Newmat Matrix Algebra library (used for SVD)
libnewmat:
	cd ${DIR_NEWMAT}
	$(MAKE) -f nm_b55.mak newmat.lib

#cctbx
libcctbx:
	@cd ${DIR_CCTBX}
	$(MAKE) -f bc32.mak lib
