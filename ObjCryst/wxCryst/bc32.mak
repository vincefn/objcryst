!include ../rules.mak

lib : trackball.obj wxAtom.obj wxCryst.obj wxCrystal.obj wxGlobalOptimObj.obj wxPowderPattern.obj wxRefinableObj.obj wxScatterer.obj wxScatteringPower.obj wxZScatterer.obj wxDiffractionSingleCrystal.obj
	tlib "libwxCryst.lib" -+trackball.obj -+wxAtom.obj -+wxCryst.obj -+wxCrystal.obj -+wxGlobalOptimObj.obj -+wxPowderPattern.obj -+wxRefinableObj.obj -+wxScatterer.obj -+wxScatteringPower.obj -+wxZScatterer.obj -+wxDiffractionSingleCrystal.obj

all: lib
