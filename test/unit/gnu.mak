ROOT_DIR ?= $(abspath $(CURDIR)/../..)
BUILD_DIR := $(ROOT_DIR)

RM ?= rm -f

ifneq ($(filter clean tidy,$(MAKECMDGOALS)),)
SKIP_BUILD_SETUP := 1
endif

ifndef SKIP_BUILD_SETUP
include $(BUILD_DIR)/ObjCryst/rules.mak
endif

BINARIES := unit_cell_smoke crystallography_workflow \
            api_spacegroup api_crystal api_molecule api_scattering \
            api_powderpattern api_cif api_optimization api_indexing \
            ground_truth

OBJ := $(addsuffix .o,$(BINARIES))

.PHONY: all build run clean tidy

all: build run

build: $(BINARIES)

run:
	./run_unit_tests.sh

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

unit_cell_smoke: unit_cell_smoke.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

crystallography_workflow: crystallography_workflow.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_spacegroup: api_spacegroup.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_crystal: api_crystal.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_molecule: api_molecule.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_scattering: api_scattering.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_powderpattern: api_powderpattern.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_cif: api_cif.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_optimization: api_optimization.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

api_indexing: api_indexing.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

ground_truth: ground_truth.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

tidy:
	@$(RM) core *.o *.dep

clean: tidy
	@$(RM) $(BINARIES)
