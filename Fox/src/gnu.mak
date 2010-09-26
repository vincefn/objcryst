BUILD_DIR = $(CURDIR)/../..
include $(BUILD_DIR)/ObjCryst/rules.mak

ifeq ($(profile),2)
%.o : %.c
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CC} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
else
%.o : %.c
	@$(MAKEDEPEND)
	${CC} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
endif

ifeq ($(profile),2)
%.o :  %.cpp libwx libnewmat libcctbx
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -D'__FOXVERSION__="$(foxversion)"' -c $< -o $@
else
%.o :  %.cpp libwx libnewmat libcctbx
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -D'__FOXVERSION__="$(foxversion)"' -c $< -o $@
endif

%.o : %.rc
	windres -i $< -o $@ --include-dir ${DIR_WXWINDOWS}/include
	
-include Fox.dep

#Main Application
Fox: Fox.o FoxJob.o GridResult.o IOSocket.o FoxServerThread.o FoxServer.o FoxClient.o WXFoxClient.o WXFoxServer.o WXGridWindow.o libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

Fox-nogui: libnewmat libcctbx libCrystVector libQuirks libRefinableObj libCryst Fox.o 
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

FoxJob.o: FoxJob.cpp
	gcc $(CFLAGS) -c FoxJob.cpp

GridResult.o: GridResult.cpp
	gcc $(CFLAGS) -c GridResult.cpp

IOSocket.o: IOSocket.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c IOSocket.cpp

FoxServerThread.o: FoxServerThread.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c FoxServerThread.cpp

FoxServer.o: FoxServer.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c FoxServer.cpp

FoxClient.o: FoxClient.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c FoxClient.cpp

WXFoxClient.o: WXFoxClient.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c WXFoxClient.cpp

WXFoxServer.o: WXFoxServer.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c WXFoxServer.cpp

WXGridWindow.o: WXGridWindow.cpp libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	gcc $(CFLAGS) -c WXGridWindow.cpp

fox: Fox


# target for making everything
.PHONY : all
all: Fox

# target for removing all object files
.PHONY : tidy
tidy::
	@${RM} core *.o *.dep

# target for removing all object files and libraries
.PHONY : clean
clean:: tidy
	@${RM} *.a Fox

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
