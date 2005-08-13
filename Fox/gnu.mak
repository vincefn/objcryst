DIR_CRYST := ObjCryst

default: all

all: Fox

Fox:
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 debug=$(debug) -C src all

Fox-nogui:
	$(MAKE) -f gnu.mak wxcryst=0 opengl=0 debug=$(debug) -C src all

doc:
	$(MAKE) -f gnu.mak -C src-doc all

clean:
	$(MAKE) -f gnu.mak -C src clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} clean

tidy:
	$(MAKE) -f gnu.mak -C src tidy
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} tidy

#install Fox in /usr/local/bin
install:
	install -m 755 src/Fox /usr/local/bin

update:
	cvs -z3 update
	cd ${DIR_CRYST} ; cvs -z3 update
   
cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
	$(MAKE) -f gnu.mak -C src cvsignore
	$(MAKE) -f gnu.mak -C src-doc cvsignore

dist:
	cd .. && tar -cjf Fox.tar.bz2 --exclude "*.o" --exclude "Fox-LastOptimizationStop.xml" --exclude ".#*" --exclude "*.a" --exclude "*.dep" --exclude "*.exe"  --exclude "Obj*.xml" --exclude "profile*" --exclude "Fox/src/Fox" Fox

rpm: dist
	mkdir ../rpmbuild
	mkdir ../rpmbuild/SOURCES
	mkdir ../rpmbuild/SPECS
	mkdir ../rpmbuild/SRPMS
	cp ObjCryst-Fox.spec ../rpmbuild/SPECS
	cd ../rpmbuild/SPECS && rpm -ba ObjCryst-Fox.spec

#Switch to ssh developer access
cvs-ext:
	perl -pi -e 's|pserver|ext|g' `find . -name Root`
	perl -pi -e 's|anonymous|vincefn|g' `find . -name Root`

#Switch to anonymous CVS
cvs-anon:
	perl -pi -e 's|ext|pserver|g' `find . -name Root`
	perl -pi -e 's|vincefn|anonymous|g' `find . -name Root`
