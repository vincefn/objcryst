BUILD_DIR = $(CURDIR)/../..
include $(BUILD_DIR)/ObjCryst/rules.mak

OBJ = Fox.o FoxJob.o GridResult.o IOSocket.o FoxServerThread.o FoxServer.o FoxClient.o WXFoxClient.o WXFoxServer.o WXGridWindow.o

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
%.o :  %.cpp $(libwx) $(libnewmat) libcctbx
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
else
%.o :  %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
endif

%.o : %.rc
	windres -i $< -o $@ --include-dir ${DIR_WXWINDOWS}/include

-include $(OBJ:.o=.dep)

libFox: $(OBJ)

#Main Application
Fox: $(libwx) $(libnewmat)  $(libfftw) $(COD_LIB) libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst
	$(MAKE) -f gnu.mak libFox
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} $(OBJ) ${LOADLIBES}

Fox-nogui: $(libnewmat) libcctbx libCrystVector libQuirks libRefinableObj libCryst Fox.o
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

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
