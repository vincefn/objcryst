include ../ObjCryst/rules.mak
DIR_CRYST := ../ObjCryst

%.o : %.c
	@$(MAKEDEPEND)
	${CC} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

%.o : %.rc
	windres -i $< -o $@ --include-dir ${DIR_WXWINDOWS}/include
	
-include Fox.dep

#wxCryst Application ( wxCrystApp_resource.o for cygwin)
Fox: Fox.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

Fox-nogui: Fox.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst
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
