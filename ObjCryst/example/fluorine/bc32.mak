!include ../../rules.mak

xray-single: xray-single.obj
	$(CPP) $(LINKFLAGS) xray-single.obj $(LINKLIBS)

xray: xray.obj
	$(CPP) $(LINKFLAGS) xray.obj $(LINKLIBS)

neutron: neutron.obj
	$(CPP) $(LINKFLAGS) neutron.obj $(LINKLIBS)

all: xray-single xray neutron
