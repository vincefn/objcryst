!include ../../rules.mak

xray-single: xray-single.obj libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst
	cd ${DIR_CRYST}\example\pbso4
	ilink32 $(LINKFLAGS) -Tpe -ap c0x32 xray-single.obj,xray-single,, import32.lib cw32mt.lib $(LINKLIBS)

xray: xray.obj libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst
	cd ${DIR_CRYST}\example\pbso4
	ilink32 $(LINKFLAGS) -Tpe -ap c0x32 xray.obj,xray,, import32.lib cw32mt.lib $(LINKLIBS)

neutron: neutron.obj libsglite libatominfo libnewmat libCrystVector libQuirks libRefinableObj libcryst
	cd ${DIR_CRYST}\example\pbso4
	ilink32 $(LINKFLAGS) -Tpe -ap c0x32 neutron.obj,neutron,, import32.lib cw32mt.lib $(LINKLIBS)

all: xray-single xray neutron
