include ../../rules.mak
DIR_CRYST = ../..

OBJ= $@.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

xray-single: xray-single.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

xray: xray.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

neutron: neutron.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

# target for making everything
.PHONY : all
all: xray-single xray neutron

default: all

.PHONY : tidy
tidy::
	@${RM} core *.o *.obj *.dep

.PHONY : clean
clean:: tidy
	@${RM} core *.a *.exe

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
