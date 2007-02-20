BUILD_DIR=$(CURDIR)/../../..
include $(BUILD_DIR)/ObjCryst/rules.mak

OBJ= $@.o

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

-include $(OBJ:.o=.dep)

fluorine: fluorine.o libCrystVector libQuirks libRefinableObj libsglite libatominfo libCryst libwxCryst
	${LINKER} ${LDFLAGS} -o $@.exe ${filter-out %.a %.so lib%, $^} ${LOADLIBES}

# target for making everything
.PHONY : all
all: fluorine

default: all

.PHONY : tidy
tidy::
	@${RM} core *.o *.obj *.dep

.PHONY : clean
clean:: tidy
	@${RM} core *.a *.exe

cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
