!include ..\ObjCryst\rules.mak

TARGET= Fox

OBJECTS = $(TARGET).obj

$(TARGET).exe: $(OBJECTS) libcctbx libnewmat libCrystVector libQuirks libRefinableObj libcryst libwxCryst
	cd $(DIR_CRYST)\..\src
	rc -r -i$(DIR_WXWINDOWS)\include -i$(DIR_WXWINDOWS)\contrib\include -i"C:\Program Files\Microsoft Visual C++ Toolkit 2003\include" -i"E:\Program Files\Microsoft SDK\include" -fo$(TARGET).res $(TARGET).rc
	link $(TARGET).obj $(TARGET).res $(LINKFLAGS) $(LINKLIBS) -out:$(TARGET).exe

all:$(TARGET).exe
