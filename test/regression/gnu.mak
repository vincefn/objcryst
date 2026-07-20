ROOT_DIR ?= $(abspath $(CURDIR)/../..)

.PHONY: all clean tidy

all:
	./fox_nogui_output_naming.sh $(ROOT_DIR)

tidy:
	@true

clean:
	@rm -rf .tmp
