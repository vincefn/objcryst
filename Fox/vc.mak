default: all

all: Fox

Fox:
	cd src
	$(MAKE) -f vc.mak all
	del ..\Fox.exe
	move Fox.exe ..\Fox.exe

clean:
	cd src
	$(MAKE) -f vc.mak clean
	cd ..\ObjCryst
	$(MAKE) -f vc.mak cleanall

update:
	cvs -z3 update
	cd ObjCryst
	cvs -z3 update
