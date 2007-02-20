BUILD_DIR=$(CURDIR)/../..
include $(BUILD_DIR)/ObjCryst/rules.mak

# GeneticAlgorithm.o Powell.o ConjugateGradient.o 
OBJ= Tracker.o Simplex.o RefinableObj.o GlobalOptimObj.o IO.o LSQNumObj.o 

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
