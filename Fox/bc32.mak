default: all

all: Fox

Fox:
	cd src
	$(MAKE) -f bc32.mak all
	del ..\Fox.exe
	move Fox.exe ..\Fox.exe

clean:
	cd src
	$(MAKE) -f bc32.mak clean
	cd ..\ObjCryst
	$(MAKE) -f bc32.mak cleanall

update:
	cd src
	cvs -z3 update
	cd ..\..\ObjCryst
	cvs -z3 update
