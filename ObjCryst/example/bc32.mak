include ../rules.mak

${BIN_DIR}/%.o : %.cpp
	mkdir -p ${BIN_DIR}
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

#test-pbso4 on powder from RoundRobin
pbso4-neutron: ${BIN_DIR}/pbso4-neutron.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

#Same as above, but on X-Ray
pbso4-xray: ${BIN_DIR}/pbso4-xray.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

#Same as above, working on extracted from F's
pbso4-xray2: ${BIN_DIR}/pbso4-xray2.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

# target for making everything
.PHONY : all
all: pbso4-neutron pbso4-xray pbso4-xray2

dep:
	makedepend *.cpp -Y -I. -I../ -I../..

# target for removing all object files
.PHONY : tidy
tidy::
	@${RM} core obj_debug/*.o obj_optimized/*.o

# target for removing all object files and libraries
.PHONY : clean
clean:: tidy
