!include ..\rules.mak

libcryst.lib : Molecule.obj ScatteringPowerFullerene.obj ScatteringCorr.obj Spacegroup.obj Scatterer.obj Atom.obj ScatteringPower.obj ZScatterer.obj Crystal.obj ScatteringData.obj DiffractionDataSingleCrystal.obj PowderPattern.obj Exception.obj geomStructFactor.obj geomStructFactor_001.obj geomStructFactor_002.obj geomStructFactor_067.obj geomStructFactor_097.obj geomStructFactor_230.obj geomStructFactor_centro.obj IO.obj ${GL_OBJ}
	tlib "libcryst.lib" -+Molecule.obj -+ScatteringPowerFullerene.obj -+ScatteringCorr.obj -+Spacegroup.obj -+Scatterer.obj -+Atom.obj -+ScatteringPower.obj -+ZScatterer.obj -+Crystal.obj -+ScatteringData.obj -+DiffractionDataSingleCrystal.obj -+PowderPattern.obj -+Exception.obj -+geomStructFactor.obj -+geomStructFactor_001.obj -+geomStructFactor_002.obj -+geomStructFactor_067.obj -+geomStructFactor_097.obj -+geomStructFactor_230.obj -+geomStructFactor_centro.obj -+IO.obj -+${GL_OBJ}

lib: libcryst.lib
