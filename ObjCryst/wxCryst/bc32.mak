!include ../rules.mak

#!include $(DIR_WXWINDOWS)\src\makeprog.b32

lib : trackball.obj wxAtom.obj wxCryst.obj wxCrystal.obj wxGlobalOptimObj.obj wxPowderPattern.obj wxRefinableObj.obj wxScatterer.obj wxScatteringPower.obj wxZScatterer.obj wxDiffractionSingleCrystal.obj
	del libwxCryst.lib
	tlib "libwxCryst.lib" +trackball.obj +wxAtom.obj +wxCryst.obj +wxCrystal.obj +wxGlobalOptimObj.obj +wxPowderPattern.obj +wxRefinableObj.obj +wxScatterer.obj +wxScatteringPower.obj +wxZScatterer.obj +wxDiffractionSingleCrystal.obj

TARGET= Fox

OBJECTS = $(TARGET).obj

#build a wx application
$(TARGET).exe: $(OBJECTS) $(TARGET).res libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst libwxCryst
	ilink32 $(LINKFLAGS) @&&!
c0w32.obj $(OBJECTS)
$(TARGET)
nul
$(LINKLIBS)

$(TARGET).res
!

all: $(TARGET).exe
