include ../rules.mak
DIR_CRYST = ..

OBJ= VFNStreamFormat.o VFNDebug.o ci_string.o

ifeq ($(profile),2)
%.o : %.cpp
	@rm -f $(*F).gcda $(*F).gcno
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
else
%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@
endif

-include $(OBJ:.o=.dep)


# libVFNStreamFormat

libQuirks.a : ${OBJ}
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}

lib: libQuirks.a

default: lib

dep:
	makedepend *.cpp -Y -I. -I../ -I../..

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} *.o *.dep

# target for removing all object files and libraries

.PHONY : clean
clean:: tidy
	@${RM} *.a *.exe

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
