ROOT_DIR ?= $(abspath $(CURDIR)/../..)
BUILD_DIR := $(ROOT_DIR)

RM ?= rm -f

ifneq ($(filter clean tidy,$(MAKECMDGOALS)),)
SKIP_BUILD_SETUP := 1
endif

ifndef SKIP_BUILD_SETUP
include $(BUILD_DIR)/ObjCryst/rules.mak
endif

OBJ := unit_cell_smoke.o crystallography_workflow.o api_surface_tests.o

.PHONY: all build run clean tidy

all: build run

build: unit_cell_smoke crystallography_workflow api_surface_tests

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

api_surface_tests: api_surface_tests.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

tidy:
	@$(RM) core *.o *.dep

clean: tidy
	@$(RM) unit_cell_smoke crystallography_workflow api_surface_tests
