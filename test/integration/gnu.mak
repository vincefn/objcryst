ROOT_DIR ?= $(abspath $(CURDIR)/../..)
BUILD_DIR := $(ROOT_DIR)

RM ?= rm -f

ifneq ($(filter clean tidy,$(MAKECMDGOALS)),)
SKIP_BUILD_SETUP := 1
endif

ifndef SKIP_BUILD_SETUP
include $(BUILD_DIR)/ObjCryst/rules.mak
endif

BINARIES := tutorial_cimetidine tutorial_pbso4

OBJ := $(addsuffix .o,$(BINARIES))

.PHONY: all build run clean tidy

all: build run

build: $(BINARIES)

run:
	./fox_nogui_cli.sh $(ROOT_DIR)
	./run_integration_tests.sh

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

tutorial_cimetidine: tutorial_cimetidine.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

tutorial_pbso4: tutorial_pbso4.o libCrystVector libQuirks libRefinableObj $(libcctbx) libCryst
	${LINKER} ${CRYST_LDFLAGS} -o $@ ${filter-out %.h %.a %.so lib%, $^} ${LOADLIBES}

tidy:
	@$(RM) core *.o *.dep

clean: tidy
	@$(RM) $(BINARIES)
	@rm -rf .tmp
