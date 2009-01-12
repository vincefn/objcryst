BUILD_DIR := $(CURDIR)/..

FOXVERSION:=1.7.6
FOXRELEASE:=SVN$(shell svnversion)_$(shell date "+%Y%m%d")

ifeq ($(shared),1)
shared-newmat=1
shared-wxgtk=1
shared-fftw=1
shared-glut=1
endif

all: Fox

default: all

Fox:
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 fftw=1 debug=$(debug) unicode=$(unicode) shared=$(shared) shared-wxgtk=$(shared-wxgtk) shared-glut=$(shared-glut) shared-newmat=$(shared-newmat) shared-fftw=$(shared-fftw) -C src Fox

Fox-nogui:
	$(MAKE) -f gnu.mak wxcryst=0 opengl=0 fftw=0 debug=$(debug) -C src Fox-nogui

doc:
	python wiki2pdf.py
	rm -Rf wikihtml
	mv FoxWiki.pdf ../

clean:
	$(MAKE) -f gnu.mak -C src clean
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/ObjCryst clean

tidy:
	$(MAKE) -f gnu.mak -C src tidy
	$(MAKE) -f gnu.mak -C ${BUILD_DIR}/ObjCryst tidy

#install Fox in /usr/local/bin
install:
	install -m 755 src/Fox /usr/local/bin

update:
	cd $(BUILD_DIR) && svn update

dist:
	cd ../.. && tar -cjf Fox.tar.bz2 --exclude "*.o" --exclude "Fox-LastOptimizationStop.xml" --exclude ".#*" --exclude "*.a" --exclude "*.dep" --exclude "*.exe"  --exclude "Obj*.xml" --exclude "profile*" --exclude "Fox/src/Fox" --exclude "*~" --exclude "cctbx" --exclude "newmat" --exclude "static-libs" --exclude "doc" --exclude "*.bak" --exclude "*.pdf" Fox

rpmmacros: ${HOME}/.rpmmacros
	#@if test -f ${HOME}/.rpmmacros ; then echo "Please remove existing ${HOME}/.rpmmacros first !" ; else mkdir -p ${HOME}/RPM/{BUILD,RPMS,SOURCES,SPECS,SRPMS,tmp} ; echo "%home %(echo ${HOME})" >> ${HOME}/.rpmmacros;       echo "%_topdir %{home}/RPM" >> ${HOME}/.rpmmacros ;echo "%_tmppath %{home}/RPM/tmp" >> ${HOME}/.rpmmacros; fi
	mkdir -p ${HOME}/RPM/{BUILD,RPMS,SOURCES,SPECS,SRPMS,tmp}
	echo "%home %(echo ${HOME})" >> ${HOME}/.rpmmacros
	echo "%_topdir %{home}/RPM" >> ${HOME}/.rpmmacros
	echo "%_tmppath %{home}/RPM/tmp" >> ${HOME}/.rpmmacros

rpm-src:
	@if test -f ../../Fox.tar.bz2 ; then echo "Using existing ../../Fox.tar.bz2" ; else make dist ;fi
	cp ObjCryst-Fox.spec ${HOME}/RPM/SPECS/
	perl -pi -e 's/FOXVERSION/$(FOXVERSION)/g' ${HOME}/RPM/SPECS/ObjCryst-Fox.spec
	perl -pi -e 's/FOXRELEASE/$(FOXRELEASE)/g' ${HOME}/RPM/SPECS/ObjCryst-Fox.spec
	cp ../../Fox.tar.bz2 ${HOME}/RPM/SOURCES/
	cd ${HOME}/RPM/SPECS && rpmbuild -bs ObjCryst-Fox.spec

rpm: rpm-src
	@if test -f ../../Fox.tar.bz2 ; then echo "Using existing ../../Fox.tar.bz2" ; else make dist ;fi
	cp ObjCryst-Fox.spec ${HOME}/RPM/SPECS/
	perl -pi -e 's/FOXVERSION/$(FOXVERSION)/g' ${HOME}/RPM/SPECS/ObjCryst-Fox.spec
	perl -pi -e 's/FOXRELEASE/$(FOXRELEASE)/g' ${HOME}/RPM/SPECS/ObjCryst-Fox.spec
	cp ../../Fox.tar.bz2 ${HOME}/RPM/SOURCES/
	cd ${HOME}/RPM/SPECS && rpmbuild -ba ObjCryst-Fox.spec
