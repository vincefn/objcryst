include ../rules.mak
DIR_CRYST := ..

OBJ= trackball.o wxCryst.o wxRefinableObj.o wxScatteringPower.o wxScatterer.o wxAtom.o wxCrystal.o wxZScatterer.o wxPowderPattern.o wxGlobalOptimObj.o

%.o : %.c
	@$(MAKEDEPEND)
	${CC} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

%.o : %.rc
	windres -i $< -o $@ --include-dir ${DIR_WXWINDOWS}/include
	
-include $(OBJ:.o=.dep) Fox.dep

#wxCryst librarry
libwxcryst.a : $(OBJ)
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}

#wxCryst Application ( wxCrystApp_resource.o for cygwin)
Fox: Fox.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so lib%, $^} ${LOADLIBES} 

fox: Fox

lib: libwxcryst.a

#install Fox in /usr/local/bin
install: Fox
	cp Fox /usr/local/bin/

# target for making everything
.PHONY : all
all: lib Fox

# target for removing all object files
.PHONY : tidy
tidy::
	@${RM} core *.o *.dep

# target for removing all object files and libraries
.PHONY : clean
clean:: tidy
	@${RM} *.a *.exe

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
