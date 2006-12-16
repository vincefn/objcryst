include ../rules.mak
DIR_CRYST := ../

OBJ= Indexing.o CIF.o ReflectionProfile.o PowderPatternBackgroundBayesianMinimiser.o Polyhedron.o ScatteringCorr.o ZScatterer.o SpaceGroup.o Scatterer.o Atom.o Molecule.o ScatteringPower.o  ScatteringPowerSphere.o Crystal.o ScatteringData.o DiffractionDataSingleCrystal.o PowderPattern.o Exception.o geomStructFactor.o geomStructFactor_001.o geomStructFactor_002.o geomStructFactor_067.o geomStructFactor_097.o geomStructFactor_230.o geomStructFactor_centro.o IO.o UnitCell.o test.o ${GL_OBJ}

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

# libcryst : 
libcryst.a : ${OBJ}
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}

-include $(OBJ:.o=.dep)
   
lib: libcryst.a

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
