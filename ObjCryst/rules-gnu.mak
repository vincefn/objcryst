# Base building directory - This must be defined in the Makefile including this one
#ROOT_DIR = ${CURDIR}
# Base ObjCryst directory
DIR_CRYST = $(BUILD_DIR)/ObjCryst

#Libraries to be statically linked are installed in $(DIR_STATIC_LIBS)/lib,
#with their headers in DIR_STATIC_LIBS)/include 
DIR_STATIC_LIBS = $(BUILD_DIR)/static-libs

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
CXX      := g++
CXXFLAGS  = ${DEPENDFLAGS} ${PROFILEFLAGS}
# FORTRAN compiler
FC     := f77
FFLAGS  = 
# linker
LINKER    := g++
LDFLAGS   = -L/usr/lib -L/usr/local/lib -L$(DIR_CRYSTVECTOR) -L$(DIR_LIBCRYST) -L$(DIR_REFOBJ) -L$(DIR_STATIC_LIBS)/lib -L$(DIR_VFNQUIRKS) -L$(DIR_WXWCRYST) -L$(DIR_TAU)/i386_linux/lib

#to automatically generate dependencies
MAKEDEPEND = gcc -MM ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} $< > $*.dep

# header files
SEARCHDIRS = -I${DIR_CRYST}/.. -I./ -I$(DIR_TAU)/include -I${DIR_CRYST} -I$(DIR_STATIC_LIBS)/include

#wxWindows flags
ifeq ($(wxcryst),1)
   WXCRYSTFLAGS = -D__WX__CRYST__ `$(WXCONFIG) --cxxflags`
   WX_LDFLAGS = -L/usr/X11R6/lib -lwxcryst `$(WXCONFIG) --libs` $(GL_WX_LIB)
else
   WXCRYSTFLAGS :=
   WX_LDFLAGS :=
endif

#Profiling
ifeq ($(profile),1) #activate profiling using TAU package
   PROFILEFLAGS := -DPROFILING_ON -DTAU_STDCXXLIB -I$(DIR_TAU)/include
   PROFILELIB := -ltau
else
   ifeq ($(profile),2) # *generate* profiling using gcc
      PROFILEFLAGS := -fprofile-generate
      PROFILELIB := -fprofile-generate
   else
      ifeq ($(profile),3) # *use* profiling using gcc
         PROFILEFLAGS := -fprofile-use
         PROFILELIB := -fprofile-use
      else
         PROFILEFLAGS :=
         PROFILELIB :=
      endif
   endif
endif

#Use static linking to wx and freeglut libraries ?
ifneq ($(shared),1)
WXCONFIG= $(DIR_CRYST)/../static-libs/bin/wx-config
else
WXCONFIG= wx-config
endif
# If using glut (freeglut)
GLUT_FLAGS= -DHAVE_GLUT
GLUT_LIB= -lglut

#Using OpenGL ?
ifeq ($(opengl),1)
GL_WX_LIB = `$(WXCONFIG) --gl-libs` -lGL -lGLU $(GLUT_LIB)
GL_FLAGS = -DOBJCRYST_GL -I/usr/X11R6/include -IGL $(GLUT_FLAGS)
else
GL_WX_LIB :=
GL_FLAGS :=
endif

#Using fftw
ifneq ($(fftw),0)
FFTW_LIB = -lfftw3f
FFTW_FLAGS = -DHAVE_FFTW
else
FFTW_LIB :=
FFTW_FLAGS :=
endif

#Set DEBUG options
# $(DIR_CRYST)/../static-libs/lib/libfftw3f.a
ifeq ($(debug),1)
   ifdef RPM_OPT_FLAGS
      # we are building a RPM !
      CPPFLAGS = ${RPM_OPT_FLAGS} 
   else
      CPPFLAGS = -g -Wall -D__DEBUG__ 
   endif
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS} ${FFTW_FLAGS}
   LOADLIBES = -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lcctbx -lnewmat ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS} ${FFTW_LIB}
else
# -march=athlon,pentiumpro
   ifdef RPM_OPT_FLAGS
      # we are building a RPM !
      CPPFLAGS = ${RPM_OPT_FLAGS} 
   else
      # Athlon XP, with auto-vectorization
      #CPPFLAGS = -O3 -w -ffast-math -march=athlon-xp -mmmx -msse -m3dnow -mfpmath=sse -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=0
      # AMD64 Opteron , with auto-vectorization
      #CPPFLAGS = -O3 -w -ffast-math -march=opteron -mmmx -msse -msse2 -m3dnow -mfpmath=sse -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=0
      #default flags
      CPPFLAGS = -O3 -w -ffast-math -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops
   endif
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS} ${FFTW_FLAGS}
   LOADLIBES = -s -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lcctbx -lnewmat ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS} ${FFTW_LIB}
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
$(DIR_STATIC_LIBS)/lib/libnewmat.a:
	cd $(BUILD_DIR) && tar -xjf newmat.tar.bz2
	$(MAKE) -f nm_gnu.mak -C $(BUILD_DIR)/newmat libnewmat.a
	mkdir -p $(DIR_STATIC_LIBS)/lib/
	cp $(BUILD_DIR)/newmat/libnewmat.a $(DIR_STATIC_LIBS)/lib/
	mkdir -p $(DIR_STATIC_LIBS)/include/newmat
	cp $(BUILD_DIR)/newmat/*.h $(DIR_STATIC_LIBS)/include/newmat/
	#rm -Rf $(BUILD_DIR)/newmat

libnewmat: $(DIR_STATIC_LIBS)/lib/libnewmat.a
     
#cctbx
$(DIR_STATIC_LIBS)/lib/libcctbx.a:
	cd $(BUILD_DIR) && tar -xjf cctbx.tar.bz2        
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/cctbx install
	#rm -Rf $(BUILD_DIR)/cctbx

libcctbx: $(DIR_STATIC_LIBS)/lib/libcctbx.a

libfftw: $(DIR_STATIC_LIBS)/lib/libfftw3f.a
	cd $(BUILD_DIR) && tar -xjf fftw.tar.bz2
	cd $(BUILD_DIR)/fftw && ./configure --enable-single --prefix $(DIR_STATIC_LIBS) && make install
	rm -Rf $(BUILD_DIR)/fftw

