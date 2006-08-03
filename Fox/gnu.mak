DIR_CRYST := ObjCryst

all: Fox

default: all

ifneq ($(shared),1)
static-libs/lib/libglut.a:
	tar -xjf freeglut-2.4.0.tar.bz2
	cd freeglut-2.4.0 && ./configure --prefix=$(PWD)/static-libs --disable-shared --disable-warnings && make install
	rm -Rf freeglut-2.4.0
static-libs/lib/libwx_base-2.6.a:
	tar -xjf wxGTK-2.6.3.tar.bz2 # wxGtK source, with "demos" "samples" "contrib" removed
	cd wxGTK-2.6.3 && ./configure --with-gtk --with-opengl --prefix=$(PWD)/static-libs --enable-optimise --disable-shared --with-gtk=any --disable-clipboard && make install
	rm -Rf wxGTK-2.6.3
endif


Fox:
ifneq ($(shared),1)
	make static-libs/lib/libglut.a static-libs/lib/libwx_base-2.6.a
endif
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 debug=$(debug) shared=$(shared) -C src all

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
	svn update
	cd ${DIR_CRYST} ; svn update
   
dist:
	cd .. && tar -cjf Fox.tar.bz2 --exclude "*.o" --exclude "Fox-LastOptimizationStop.xml" --exclude ".#*" --exclude "*.a" --exclude "*.dep" --exclude "*.exe"  --exclude "Obj*.xml" --exclude "profile*" --exclude "Fox/src/Fox" --exclude "*~" --exclude "static-limbs" --exclude "doc" --exclude "*.bak" --exclude "*.pdf" Fox

rpm: dist
	cp ObjCryst-Fox.spec /usr/src/RPM/SPECS/
	cp ../Fox.tar.bz2 /usr/src/RPM/SOURCES/
	cd /usr/src/RPM/SPECS && rpm -ba ObjCryst-Fox.spec
