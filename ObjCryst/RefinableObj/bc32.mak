!include ../rules.mak

libRefinableObj.lib : RefinableObj.obj LSQNumObj.obj GlobalOptimObj.obj IO.obj
	del libRefinableObj.lib
	tlib "libRefinableObj.lib" +RefinableObj.obj +LSQNumObj.obj +GlobalOptimObj.obj +IO.obj 

lib: libRefinableObj.lib

all: lib
