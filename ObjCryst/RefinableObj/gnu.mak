include ../rules.mak
DIR_CRYST = ..

OBJ=RefinableObj.o LSQNumObj.o GlobalOptimObj.o IO.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} -c $< -o $@
   
libRefinableObj.a : ${OBJ}
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}


-include $(OBJ:.o=.dep)

lib: libRefinableObj.a

default: lib

all: lib

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core *.o *.dep

# target for removing all object files and libraries

.PHONY : clean
clean:: tidy
	@${RM} *.a *.exe

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
