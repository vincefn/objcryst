include ../rules.mak
DIR_CRYST = ..

OBJ= VFNStreamFormat.o VFNDebug.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} -c $< -o $@

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
