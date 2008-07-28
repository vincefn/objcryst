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
%.o : %.cpp libwx
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -D'__FOXVERSION__="$(foxversion)"' -c $< -o $@
else
%.o : %.cpp #libwx
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -D'__FOXVERSION__="$(foxversion)"' -c $< -o $@
endif

%.o : %.rc
	windres -i $< -o $@ --include-dir ${DIR_WXWINDOWS}/include
	
-include Fox.dep

#Main Application
Fox: Fox.o libwx libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst libwxCryst libfftw
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

Fox-nogui: Fox.o libnewmat libCrystVector libQuirks libRefinableObj libcctbx libCryst
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

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
