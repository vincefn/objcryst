include ../rules.mak
DIR_CRYST := ..

#wxGeneticAlgorithm.o
OBJ= wxScatteringPowerSphere.o trackball.o wxDiffractionSingleCrystal.o wxCryst.o wxRefinableObj.o wxScatteringPower.o wxScatterer.o wxAtom.o wxCrystal.o wxZScatterer.o wxPowderPattern.o wxGlobalOptimObj.o

%.o : %.c
	@$(MAKEDEPEND)
	${CC} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

#wxCryst librarry
libwxcryst.a : $(OBJ)
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}

lib: libwxcryst.a

# target for making everything
.PHONY : all
all: lib

# target for removing all object files
.PHONY : tidy
tidy::
	@${RM} core *.o *.dep

# target for removing all object files and libraries
.PHONY : clean
clean:: tidy
	@${RM} *.a

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
