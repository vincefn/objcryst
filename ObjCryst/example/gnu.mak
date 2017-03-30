SUBDIRS= pbso4

all:
	@for d in $(SUBDIRS); do (cd $$d && make -f gnu.mak all debug=$(debug) opengl=$(opengl) wxcryst=$(wxcryst)); done

default:all

# target for removing all object files
.PHONY : tidy
tidy::
	@for d in $(SUBDIRS); do (cd $$d && make -f gnu.mak tidy); done

# target for removing all object files and libraries
.PHONY : clean
clean::
	@for d in $(SUBDIRS); do (cd $$d && make -f gnu.mak clean); done
