# Base ObjCryst directory
DIR_CRYST=C:\Dev\Fox\ObjCryst
# wxWindows Directory
DIR_WXWINDOWS = C:\Dev\wxWindows

#You should not need to modify anything below
###################################################################
#use wxCryst ?
wxcryst=1
#debugging ? (only activates debugging messages)
debug=0

#Other directories
DIR_ATOMINFO = $(DIR_CRYST)\..\atominfo
DIR_BLITZ = $(DIR_CRYST)\..\blitz
DIR_CRYSTVECTOR = $(DIR_CRYST)\CrystVector
DIR_EXAMPLE = $(DIR_CRYST)\example
DIR_LIBCRYST = $(DIR_CRYST)\ObjCryst
DIR_NEWMAT = $(DIR_CRYST)\..\newmat
DIR_REFOBJ = $(DIR_CRYST)\RefinableObj
DIR_SGINFO = $(DIR_CRYST)\..\sginfo
DIR_SGLITE = $(DIR_CRYST)\..\sglite
DIR_TAU = $(DIR_CRYST)\..\tau
DIR_VFNQUIRKS = $(DIR_CRYST)\Quirks
DIR_WXWCRYST = $(DIR_CRYST)\wxCryst

DIR_DOC = $(DIR_CRYST)/doc

###################################################################
#
# compiler options

#CCOPTS = -nologo -O2 -MT -GX -I"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include"
#CCOPTS = -nologo -O2 -G7 -arch:SSE -GX -Ox -GR  -MT -I"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include"
CCOPTS = -nologo -O2 -G6 -GX -Ox -GR -MT -I"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include"


# header files
SEARCHDIRS = -I$(DIR_CRYST) -I. -I.. -I..\.. -I..\..\.. -I$(DIR_BLITZ) -I$(DIR_NEWMAT) -I$(DIR_CRYST) -I$(DIR_SGLITE) -I$(DIR_ATOMINFO) -I$(DIR_WXWINDOWS)\include -I"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include" -I"E:\Program Files\Microsoft SDK\include"

#debug ?
!if ($(debug)==1)
CPPDEBUGFLAGS = -D__DEBUG__
!else
CPPDEBUGFLAGS = 
!endif

#wxWindows flags
!if $(wxcryst)==1
#@$(DIR_WXWINDOWS)\src\msw\wxw32.cfg
WXCRYSTFLAGS = -D__WX__CRYST__ -DOBJCRYST_GL
WX_LDFLAGS = -LIBPATH:$(DIR_WXWCRYST) -LIBPATH:$(DIR_WXWINDOWS)\lib
WINLIBS=  kernel32.lib user32.lib gdi32.lib comdlg32.lib winspool.lib winmm.lib shell32.lib oldnames.lib comctl32.lib odbc32.lib ole32.lib oleaut32.lib uuid.lib rpcrt4.lib advapi32.lib wsock32.lib 
WX_LIBS = libwxcryst.lib wxmsw.lib jpeg.lib tiff.lib png.lib zlib.lib opengl32.lib glu32.lib $(WINLIBS)
!else
WXCRYSTFLAGS = 
WX_LDFLAGS = 
WX_LIBS = 
!endif


#Compile flags:
CFLAGS= $(CCOPTS) $(DBGOPT)  $(ENVOPTS) $(DEFOPTS) $(THROPTS) $(CCLINKOPT) $(SEARCHDIRS) -DOBJCRYST_GL $(WXCRYSTFLAGS) $(CPPDEBUGFLAGS)
CPPFLAGS= $(CCOPTS) $(DBGOPT)  $(ENVOPTS) $(DEFOPTS) $(THROPTS) $(CCLINKOPT) $(SEARCHDIRS) -DOBJCRYST_GL $(WXCRYSTFLAGS) $(CPPDEBUGFLAGS)


#Link
#-machine:i386 -subsystem:windows,4
#-entry:WinMainCRTStartup
#LINKFLAGS= /INCREMENTAL:NO /NOLOGO -subsystem:windows $(WX_LDFLAGS) -LIBPATH:$(DIR_ATOMINFO) -LIBPATH:$(DIR_CRYSTVECTOR) -LIBPATH:$(DIR_LIBCRYST) -LIBPATH:$(DIR_NEWMAT) -LIBPATH:$(DIR_BLITZ)lib -LIBPATH:$(DIR_REFOBJ) -LIBPATH:$(DIR_SGLITE) -LIBPATH:$(DIR_VFNQUIRKS) -LIBPATH:"C:\Program Files\Microsoft Visual C++ Toolkit 2003\Lib" -LIBPATH:"E:\Program Files\Microsoft SDK\Lib"
LINKFLAGS= -INCREMENTAL:NO -nologo -subsystem:console -entry:WinMainCRTStartup $(WX_LDFLAGS) -LIBPATH:$(DIR_ATOMINFO) -LIBPATH:$(DIR_CRYSTVECTOR) -LIBPATH:$(DIR_LIBCRYST) -LIBPATH:$(DIR_NEWMAT) -LIBPATH:$(DIR_BLITZ)lib -LIBPATH:$(DIR_REFOBJ) -LIBPATH:$(DIR_SGLITE) -LIBPATH:$(DIR_VFNQUIRKS) -LIBPATH:"C:\Program Files\Microsoft Visual C++ Toolkit 2003\Lib" -LIBPATH:"E:\Program Files\Microsoft SDK\Lib"

LINKLIBS= libatominfo.lib libsglite.lib libcrystvector.lib libquirks.lib librefinableobj.lib libcryst.lib $(WX_LIBS)
#RCFLAGS= -r -i$(MAKEDIR)\..\include;$(MAKEDIR)\..\include\windows



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
	cd $(DIR_LIBCRYST)
	$(MAKE) -f vc.mak lib

#wxCryst++
!if $(wxcryst)==1
libwxCryst:
   cd $(DIR_WXWCRYST)
	$(MAKE) -f vc.mak  lib
!else
libwxCryst:
!endif

#Vector computation library
libCrystVector:
	cd $(DIR_CRYSTVECTOR)
	$(MAKE) -f vc.mak lib

#Quirks, including a (crude) library to display float, vectors, matrices, strings with some formatting..
libQuirks:
	cd $(DIR_VFNQUIRKS)
	$(MAKE) -f vc.mak lib

#Library to take care of refinable parameters, plus Global optimization and Least Squares refinements
libRefinableObj:
	cd $(DIR_REFOBJ)
	$(MAKE) -f vc.mak lib

#Newmat Matrix Algebra library (used for SVD)
libnewmat:
	cd $(DIR_NEWMAT)
	$(MAKE) -f nm_m6.mak newmat.lib

#SgLite -Spacegroup Lib
libsglite:
	@cd $(DIR_SGLITE)
	$(MAKE) -f vc.mak lib

#AtomInfo
libatominfo:
	cd $(DIR_ATOMINFO)
	$(MAKE) -f vc.mak lib
