!include ../rules.mak

libRefinableObj.lib : RefinableObj.obj LSQNumObj.obj GlobalOptimObj.obj IO.obj
	lib -OUT:libRefinableObj.lib RefinableObj.obj LSQNumObj.obj GlobalOptimObj.obj IO.obj 

lib: libRefinableObj.lib

all: lib
