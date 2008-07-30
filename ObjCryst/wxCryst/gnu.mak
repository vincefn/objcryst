BUILD_DIR=$(CURDIR)/../..
include $(BUILD_DIR)/ObjCryst/rules.mak

#wxGeneticAlgorithm.o
OBJ= wxLSQ.o wxTrackerGraph.o wxMultiGraph.o wxScatteringPowerSphere.o trackball.o wxDiffractionSingleCrystal.o wxCryst.o wxRefinableObj.o wxScatteringPower.o wxScatterer.o wxAtom.o wxMolecule.o wxCrystal.o wxZScatterer.o wxPowderPattern.o wxGlobalOptimObj.o mpVector.o MC.o

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
%.o : %.cpp
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
else
%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
endif

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
