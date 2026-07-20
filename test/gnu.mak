ROOT_DIR := $(abspath $(CURDIR)/..)
OBJCRYST_DIR := $(ROOT_DIR)/ObjCryst
LOG_DIR := $(CURDIR)/.tmp/build-logs

TEST_FLAGS := wxcryst=0 opengl=0 fftw=0 cod=0
FOX_BUILD_VERBOSE ?= 0

SUBDIRS := unit integration regression

.PHONY: all setup build-fox-nogui $(SUBDIRS) clean tidy

all: setup unit build-fox-nogui integration regression

setup:
	@ln -sf rules-gnu.mak $(OBJCRYST_DIR)/rules.mak

unit:
	@mkdir -p $(LOG_DIR)
	@$(MAKE) --no-print-directory -s -f gnu.mak -C unit build $(TEST_FLAGS) ROOT_DIR=$(ROOT_DIR) >$(LOG_DIR)/unit-build.log 2>&1 || (cat $(LOG_DIR)/unit-build.log && exit 1)
	@$(MAKE) --no-print-directory -s -f gnu.mak -C unit run ROOT_DIR=$(ROOT_DIR)

build-fox-nogui:
	@mkdir -p $(LOG_DIR)
ifeq ($(FOX_BUILD_VERBOSE),1)
	@echo "[test] Building Fox-nogui (verbose mode enabled)"
	@time $(MAKE) --no-print-directory -f gnu.mak -C $(ROOT_DIR)/Fox Fox-nogui debug=$(debug)
else
	@$(MAKE) --no-print-directory -s -f gnu.mak -C $(ROOT_DIR)/Fox Fox-nogui debug=$(debug) >$(LOG_DIR)/fox-build.log 2>&1 || (cat $(LOG_DIR)/fox-build.log && exit 1)
endif

integration:
	@mkdir -p $(LOG_DIR)
	@$(MAKE) --no-print-directory -s -f gnu.mak -C integration build $(TEST_FLAGS) ROOT_DIR=$(ROOT_DIR) >$(LOG_DIR)/integration-build.log 2>&1 || (cat $(LOG_DIR)/integration-build.log && exit 1)
	@$(MAKE) --no-print-directory -s -f gnu.mak -C integration run ROOT_DIR=$(ROOT_DIR)

regression:
	@$(MAKE) --no-print-directory -s -f gnu.mak -C regression ROOT_DIR=$(ROOT_DIR)

tidy:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -s -f gnu.mak -C $$d tidy ROOT_DIR=$(ROOT_DIR); done

clean:
	@for d in $(SUBDIRS); do $(MAKE) --no-print-directory -s -f gnu.mak -C $$d clean ROOT_DIR=$(ROOT_DIR); done
	@rm -rf $(CURDIR)/.tmp
