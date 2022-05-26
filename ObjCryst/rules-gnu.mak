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

# DO we want to use shared libraries for wxGTK, freeglut, fftw & newmat ?
# User can also use shared libraries only for some libraries, by using
# "make shared-wxgtk=1" instead of "make shared=1"
ifeq ($(shared),1)
shared-newmat=1
shared-wxgtk=1
shared-fftw=1
shared-glut=1
endif
### Rules for Linux & GCC
# C compiler
#CC     := gcc
CFLAGS  = ${DEPENDFLAGS}
# C++ compiler
#CXX      := g++
CXXFLAGS  = ${DEPENDFLAGS} ${PROFILEFLAGS}
# FORTRAN compiler
FC     := f77
FFLAGS  =
# linker
LINKER    := ${CXX}
CRYST_LDFLAGS   = ${LDFLAGS} -L/usr/lib -L/usr/local/lib -L$(DIR_CRYSTVECTOR) -L$(DIR_LIBCRYST) -L$(DIR_REFOBJ) -L$(DIR_STATIC_LIBS)/lib -L$(DIR_VFNQUIRKS) -L$(DIR_WXWCRYST) -L$(DIR_TAU)/x86_64/lib

#to automatically generate dependencies
MAKEDEPEND = gcc -MM ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} $< > $*.dep

# header files
SEARCHDIRS = -I$(DIR_TAU)/include -I${DIR_CRYST}/.. -I$(DIR_STATIC_LIBS)/include

#wxWindows flags
ifeq ($(wxcryst),1)
   WXCRYSTFLAGS = -D__WX__CRYST__ `$(WXCONFIG) --cxxflags`
   WX_LDFLAGS = -L/usr/X11R6/lib -lwxcryst `$(WXCONFIG) --libs adv,core,base,net` $(GL_WX_LIB)
else
   WXCRYSTFLAGS :=
   WX_LDFLAGS :=
endif

#Profiling
ifeq ($(profile),1) #activate profiling using TAU package
   DIR_TAU=$(BUILD_DIR)/../../utils/tau
   PROFILEFLAGS := -DPROFILING_ON -DTAU_GNU -DTAU_DOT_H_LESS_HEADERS -DTAU_LINUX_TIMERS -DTAU_LARGEFILE -D_LARGEFILE64_SOURCE -DHAVE_TR1_HASH_MAP -UTAU_MPI -I$(DIR_TAU)/include
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

#Use static linking to wx and freeglut libraries ? Unicode
ifneq ($(shared-wxgtk),1)
WXCONFIG= $(DIR_CRYST)/../static-libs/bin/wx-config --unicode=yes
else
WXCONFIG= wx-config --unicode=yes
endif

# Which command to use for download ?
WGET=$(shell which wget 2>/dev/null)
ifneq ($(WGET),)
DOWNLOAD_COMMAND=wget
else
DOWNLOAD_COMMAND=curl -L -O
endif

# If using glut (freeglut)
GLUT_FLAGS= -DHAVE_GLUT
GLUT_LIB= -lglut

#Specifiy REAL= float or double
ifeq ($(real_double),1)
REAL_FLAG= -DREAL=double
else
REAL_FLAG= -DREAL=float
endif
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

ifneq ($(sse),0)
SSE_FLAGS = -DHAVE_SSE_MATHFUN -DUSE_SSE2 -march=native
else
SSE_FLAGS =
endif

ifneq ($(shared-newmat),1)
LDNEWMAT := $(DIR_STATIC_LIBS)/lib/libnewmat.a
else
LDNEWMAT := -lnewmat
endif

#Using COD, with MySQL
ifneq ($(cod),0)
 COD_FLAGS = -D__FOX_COD__
else
 COD_FLAGS :=
endif

ifeq ($(shared_libcryst),1)
 CPPFLAGS = -O3 -w -fPIC -g -ffast-math -fstrict-aliasing -pipe -funroll-loops ${SSE_FLAGS}
 DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS} ${FFTW_FLAGS} ${REAL_FLAG}
else
 ifeq ($(debug),1)
 #Set DEBUG options
   ifdef RPM_OPT_FLAGS
      # we are building a RPM !
      CPPFLAGS = ${RPM_OPT_FLAGS}
   else
      CPPFLAGS = -g -Wall -D__DEBUG__ ${SSE_FLAGS} ${COD_FLAGS}
   endif
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS} ${FFTW_FLAGS} ${REAL_FLAG}
   LOADLIBES = -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lcctbx ${LDNEWMAT} ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS} ${FFTW_LIB}
 else
   ifdef RPM_OPT_FLAGS
      # we are building a RPM !
      CPPFLAGS = ${RPM_OPT_FLAGS}
   else
      #default flags - use "sse=1" to enable SSE optimizations
      CPPFLAGS = -O3 -w -ffast-math -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops -ftree-vectorize ${SSE_FLAGS} ${COD_FLAGS}
   endif
   DEPENDFLAGS = ${SEARCHDIRS} ${GL_FLAGS} ${WXCRYSTFLAGS} ${FFTW_FLAGS} ${REAL_FLAG}
   LOADLIBES = -lm -lcryst -lCrystVector -lQuirks -lRefinableObj -lcctbx ${LDNEWMAT} ${PROFILELIB} ${GL_LIB} ${WX_LDFLAGS} ${FFTW_LIB}
 endif
endif
# Add to statically link: -nodefaultlibs -lgcc /usr/lib/libstdc++.a

######################################################################
#####################      LIBRAIRIES         ########################
######################################################################
#Newmat Matrix Algebra library (used for SVD)

$(BUILD_DIR)/newmat.tar.bz2:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) https://github.com/vincefn/objcryst/releases/download/v2021-3rdPartyLibs/newmat.tar.bz2

$(DIR_STATIC_LIBS)/lib/libnewmat.a: $(BUILD_DIR)/newmat.tar.bz2
	cd $(BUILD_DIR) && tar -xjf newmat.tar.bz2
	$(MAKE) -f nm_gnu.mak -C $(BUILD_DIR)/newmat libnewmat.a
	mkdir -p $(DIR_STATIC_LIBS)/lib/
	cp $(BUILD_DIR)/newmat/libnewmat.a $(DIR_STATIC_LIBS)/lib/
	mkdir -p $(DIR_STATIC_LIBS)/include/newmat
	cp $(BUILD_DIR)/newmat/*.h $(DIR_STATIC_LIBS)/include/newmat/
	#rm -Rf $(BUILD_DIR)/newmat

ifneq ($(shared-newmat),1)
libnewmat= $(DIR_STATIC_LIBS)/lib/libnewmat.a
else
libnewmat=
endif

$(BUILD_DIR)/freeglut.tar.bz2:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) https://github.com/vincefn/objcryst/releases/download/v2021-3rdPartyLibs/freeglut.tar.bz2

$(BUILD_DIR)/static-libs/lib/libglut.a: $(BUILD_DIR)/freeglut.tar.bz2
	cd $(BUILD_DIR) && tar -xjf freeglut.tar.bz2
	cd $(BUILD_DIR)/freeglut && ./configure --prefix=$(BUILD_DIR)/static-libs --disable-shared --disable-warnings --x-includes=/usr/X11R6/include/ && $(MAKE) install
	rm -Rf freeglut

ifeq ($(opengl),1)
ifneq ($(shared-glut),1)
libfreeglut= $(DIR_STATIC_LIBS)/lib/libglut.a
else
libfreeglut=
endif
else
libfreeglut=
endif

$(BUILD_DIR)/wxWidgets-3.1.6.tar.bz2:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) https://github.com/wxWidgets/wxWidgets/releases/download/v3.1.6/wxWidgets-3.1.6.tar.bz2

$(BUILD_DIR)/static-libs/include/wx-3.1/wx/wx.h: $(BUILD_DIR)/wxWidgets-3.1.6.tar.bz2
	cd $(BUILD_DIR) && rm -Rf wxWidgets-3.1.6 && tar -xjf wxWidgets-3.1.6.tar.bz2
	cd $(BUILD_DIR)/wxWidgets-3.1.6 && ./configure --with-gtk --with-opengl --disable-glcanvasegl --prefix=$(BUILD_DIR)/static-libs --enable-unicode  --enable-optimise --disable-shared --x-includes=/usr/X11R6/include/ && $(MAKE) install
	rm -Rf $(BUILD_DIR)/wxWidgets-3.1.6

ifneq ($(wxcryst),0)
ifneq ($(shared-wxgtk),1)
libwx = $(BUILD_DIR)/static-libs/include/wx-3.1/wx/wx.h
else
libwx=
endif
else
libwx=
endif

#cctbx
$(BUILD_DIR)/cctbx.tar.bz2:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) https://github.com/vincefn/objcryst/releases/download/v2021-3rdPartyLibs/cctbx.tar.bz2

$(DIR_STATIC_LIBS)/lib/libcctbx.a: $(BUILD_DIR)/cctbx.tar.bz2
	mkdir -p $(DIR_STATIC_LIBS)/lib/ $(DIR_STATIC_LIBS)/include/
	cd $(BUILD_DIR) && tar -xjf cctbx.tar.bz2
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/cctbx install
	#ln -sf $(BUILD_DIR)/boost $(DIR_STATIC_LIBS)/include/
	#rm -Rf $(BUILD_DIR)/cctbx

libcctbx: $(DIR_STATIC_LIBS)/lib/libcctbx.a

$(BUILD_DIR)/fftw-3.3.9.tar.gz:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) http://fftw.org/fftw-3.3.9.tar.gz

$(DIR_STATIC_LIBS)/lib/libfftw3f.a: $(BUILD_DIR)/fftw-3.3.9.tar.gz
	cd $(BUILD_DIR) && tar -xzf fftw-3.3.9.tar.gz
	cd $(BUILD_DIR)/fftw-3.3.9 && ./configure --enable-single --prefix $(DIR_STATIC_LIBS) && $(MAKE) install
	rm -Rf $(BUILD_DIR)/fftw-3.3.9

ifneq ($(fftw),0)
ifneq ($(shared-fftw),1)
libfftw = $(DIR_STATIC_LIBS)/lib/libfftw3f.a
else
libfftw=
endif
else
libfftw=
endif

# Boost library - used for cctbx and Fox. 
# TODO: clarify install to see if we use boost from cctbx or independently installed...
$(BUILD_DIR)/boost_1_68_0.tar.bz2:
	cd $(BUILD_DIR) && $(DOWNLOAD_COMMAND) https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.bz2

libboost:$(BUILD_DIR)/boost_1_68_0.tar.bz2
	cd $(BUILD_DIR) && tar -xjf boost_1_68_0.tar.bz2
	cd $(BUILD_DIR)/boost_1_68_0 && rsync -ar boost --exclude="accumulators/" --exclude="asio/" --exclude="asign/" --exclude="bimap/" --exclude="circular_buffer/" --exclude="concept_check/" --exclude="dynamic_bitset/" --exclude="filesystem/" --exclude="flyweight/" --exclude="function_types/" --exclude="gil/" --exclude="graph/" --exclude="interprocess/" --exclude="intrusive/" --exclude="iostreams/" --exclude="lambda/" --exclude="logic/" --exclude="mpi/" --exclude="parameter/" --exclude="pending/" --exclude="pool/" --exclude="program_options/" --exclude="proto/" --exclude="ptr_container/" --exclude="python/" --exclude="random/" --exclude="regex/" --exclude="signals/" --exclude="spirit/" --exclude="statechart/" --exclude="system/" --exclude="test/" --exclude="thread/" --exclude="units/" --exclude="unordered/" --exclude="wave/" --exclude="xpressive/" --filter="+ */" --filter="+ *.hpp"  --filter="+ *.ipp" $(DIR_STATIC_LIBS)/include/
	rm -Rf $(BUILD_DIR)/boost_1_68_0

#ObjCryst++
libCryst: $(libwx) libcctbx
	$(MAKE) -f gnu.mak -C ${DIR_LIBCRYST} lib

libcryst: libCryst

#wxCryst++
libwxCryst: $(libwx) $(libfreeglut) $(libfftw) libcctbx
	$(MAKE) -f gnu.mak -C ${DIR_WXWCRYST} lib

#Vector computation library
libCrystVector: $(libwx)
	$(MAKE) -f gnu.mak -C ${DIR_CRYSTVECTOR} lib

#Quirks, including a (crude) library to display float, vectors, matrices, strings with some formatting..
libQuirks: $(libwx)
	$(MAKE) -f gnu.mak -C ${DIR_VFNQUIRKS} lib

#Library to take care of refinable parameters, plus Global optimization and Least Squares refinements
libRefinableObj:$(libnewmat) $(libwx) libcctbx
	$(MAKE) -f gnu.mak -C ${DIR_REFOBJ}/ lib
