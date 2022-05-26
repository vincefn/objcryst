BUILD_DIR := $(CURDIR)/..

DESTDIR ?= /usr/local/bin

ifeq ($(shared),1)
shared-newmat=1
shared-wxgtk=1
shared-fftw=1
shared-glut=1
shared-mysql=1
endif

all: Fox

default: all

Fox:
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 fftw=1 debug=$(debug) unicode=$(unicode) cod=$(cod) shared=$(shared) shared-wxgtk=$(shared-wxgtk) shared-glut=$(shared-glut) shared-newmat=$(shared-newmat) shared-fftw=$(shared-fftw) shared-mysql=$(shared-mysql) -C src Fox

Fox-nogui:
	$(MAKE) -f gnu.mak wxcryst=0 opengl=0 fftw=0 cod=0 debug=$(debug) -C src Fox-nogui

doc:
	python wiki2pdf.py
	rm -Rf wikihtml
	mv FoxWiki.pdf ../

clean: tidy
	$(MAKE) -f gnu.mak -C src clean
	$(MAKE) -f gnu.mak -C $(BUILD_DIR)/ObjCryst clean
	${RM} -Rf ${BUILD_DIR}/static-libs/*

tidy:
	$(MAKE) -f gnu.mak -C src tidy
	$(MAKE) -f gnu.mak -C ${BUILD_DIR}/ObjCryst tidy
	${RM} -Rf ${BUILD_DIR}/fftw-3.3.4 ${BUILD_DIR}/newmat ${BUILD_DIR}/wxWidgets-3.1.6

#install Fox in /usr/local/bin
install:
	install -Dm 755 src/Fox $(DESTDIR)

update:
	cd $(BUILD_DIR) && git pull

dist:
	cd ../.. && tar -cjf Fox.tar.bz2 --exclude "*.o" --exclude "Fox-LastOptimizationStop.xml" --exclude ".#*" --exclude "*.a" --exclude "*.dep" --exclude "*.exe"  --exclude "Obj*.xml" --exclude "profile*" --exclude "Fox/src/Fox" --exclude "*~" --exclude "cctbx" --exclude "newmat" --exclude "static-libs" --exclude "doc" --exclude "*.bak" --exclude "*.pdf" Fox
