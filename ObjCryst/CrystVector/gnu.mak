include ../rules.mak

DIR_CRYST = ..

##Using blitz library ??
##Make it (in blitz dir)  using : ./configure --with-cxx=gcc --enable-optimize
##and then : make lib
#ifeq (${VFNVECTOR_USE_BLITZ},1)
#   C_BLITZFLAG := -D__VFN_VECTOR_USE_BLITZ__
#   AR_BLITZ := ${DIR_BLITZ}/lib/libblitz.a
#else
   C_BLITZFLAG :=
   AR_BLITZ := 
#endif

OBJ=CrystVector.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} ${C_BLITZFLAG} -c $< -o $@

-include $(OBJ:.o=.dep)
   
libCrystVector.a : ${OBJ}
	${AR} crs $@ ${AR_BLITZ} ${filter-out %.a %.so, $^}

lib: libCrystVector.a

dep:
	makedepend *.cpp -Y -I. -I../ 

default: lib

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
