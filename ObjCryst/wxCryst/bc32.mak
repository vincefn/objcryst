!include ../rules.mak

lib : trackball.obj wxScatteringPowerSphere.obj wxAtom.obj wxCryst.obj wxCrystal.obj wxGlobalOptimObj.obj wxPowderPattern.obj wxRefinableObj.obj wxScatterer.obj wxScatteringPower.obj wxZScatterer.obj wxDiffractionSingleCrystal.obj mpVector.obj MC.obj
	tlib "libwxCryst.lib" -+trackball.obj -+wxScatteringPowerSphere.obj -+wxAtom.obj -+wxCryst.obj -+wxCrystal.obj -+wxGlobalOptimObj.obj -+wxPowderPattern.obj -+wxRefinableObj.obj -+wxScatterer.obj -+wxScatteringPower.obj -+wxZScatterer.obj -+wxDiffractionSingleCrystal.obj -+mpVector.obj -+MC.obj

all: lib
