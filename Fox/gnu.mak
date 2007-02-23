BUILD_DIR := $(CURDIR)/..

FOXVERSION:=1.7.0
FOXRELEASE:=SVN$(shell svnversion)_$(shell date "+%Y%m%d")

all: Fox

default: all

ifneq ($(shared),1)
$(BUILD_DIR)/static-libs/lib/libglut.a:
	cd $(BUILD_DIR) && tar -xjf freeglut.tar.bz2
	cd $(BUILD_DIR)/freeglut && ./configure --prefix=$(BUILD_DIR)/static-libs --disable-shared --disable-warnings --x-includes=/usr/X11R6/include/ && make install
	rm -Rf freeglut
$(BUILD_DIR)/static-libs/bin/wx-config:
	cd $(BUILD_DIR) && tar -xjf wxGTK.tar.bz2 # wxGtK source, with "demos" "samples" "contrib" removed
	cd $(BUILD_DIR)/wxGTK && ./configure --with-gtk --with-opengl --prefix=$(BUILD_DIR)/static-libs --disable-unicode --enable-optimise --disable-shared --with-gtk=any --disable-clipboard --x-includes=/usr/X11R6/include/ && make install
	rm -Rf wxGTK
endif

Fox:
ifneq ($(shared),1)
	make $(BUILD_DIR)/static-libs/lib/libglut.a $(BUILD_DIR)/static-libs/bin/wx-config
endif
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 debug=$(debug) shared=$(shared) foxversion="$(FOXVERSION)-$(FOXRELEASE)" -C src all

Fox-nogui:
	$(MAKE) -f gnu.mak wxcryst=0 opengl=0 debug=$(debug) foxversion="$(FOXVERSION)-$(FOXRELEASE)" -C src all

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
