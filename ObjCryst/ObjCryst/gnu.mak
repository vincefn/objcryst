include ../rules.mak
DIR_CRYST := ../

OBJ= Molecule.o ScatteringCorr.o ZScatterer.o SpaceGroup.o Scatterer.o Atom.o ScatteringPower.o  ScatteringPowerSphere.o Crystal.o ScatteringData.o DiffractionDataSingleCrystal.o PowderPattern.o Exception.o geomStructFactor.o geomStructFactor_001.o geomStructFactor_002.o geomStructFactor_067.o geomStructFactor_097.o geomStructFactor_230.o geomStructFactor_centro.o IO.o UnitCell.o ${GL_OBJ}

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

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
