!include ../rules.mak

libRefinableObj.lib : Simplex.obj RefinableObj.obj LSQNumObj.obj GlobalOptimObj.obj IO.obj
	lib -OUT:libRefinableObj.lib Simplex.obj RefinableObj.obj LSQNumObj.obj GlobalOptimObj.obj IO.obj 

lib: libRefinableObj.lib

all: lib
