!include ..\ObjCryst\rules.mak

TARGET= Fox

OBJECTS = $(TARGET).obj

#build a wx application
$(TARGET).exe: $(OBJECTS) $(TARGET).res libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst libwxCryst
	cd ${DIR_CRYST}\..\src
	ilink32 $(LINKFLAGS) -LC:${DIR_CRYST}\..\Fox\src @&&!
c0w32.obj $(OBJECTS)
$(TARGET)

$(LINKLIBS)

$(TARGET).res
!

all: $(TARGET).exe
