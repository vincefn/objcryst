!include ..ObjCryst\rules.mak

TARGET= Fox

OBJECTS = $(TARGET).obj

#build a wx application
$(TARGET).exe: $(OBJECTS) $(TARGET).res libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst libwxCryst
	ilink32 $(LINKFLAGS) @&&!
c0w32.obj $(OBJECTS)
$(TARGET)

$(LINKLIBS)

$(TARGET).res
!

all: $(TARGET).exe
