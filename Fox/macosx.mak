# NOTE: beware of multiple installed libraries, when using homebrew or macports.
# For example iconv in /usr/lib and /opt/local/lib, leading to link errors...
# For a clean compile exluding non-system libraries, use "PATH=/usr/bin:/bin:/usr/sbin:/sbin make Fox"

# Which command to use for download ?
CURL=$(shell which curl 2>/dev/null)
ifneq ($(CURL),)
DOWNLOAD_COMMAND=curl -L -O
else
DOWNLOAD_COMMAND=wget
endif

../cctbx.tar.bz2:
	cd .. && $(DOWNLOAD_COMMAND)  https://github.com/vincefn/objcryst/releases/download/v2024.1/cctbx.tar.bz2

../cctbx: ../cctbx.tar.bz2
	cd .. && tar -xjf cctbx.tar.bz2

../newmat.tar.bz2:
	cd .. && $(DOWNLOAD_COMMAND)  https://github.com/vincefn/objcryst/releases/download/v2021-3rdPartyLibs/newmat.tar.bz2

../newmat: ../newmat.tar.bz2
	cd .. && tar -xjf newmat.tar.bz2

../fftw-3.3.10.tar.gz:
	cd .. &&  curl -O http://fftw.org/fftw-3.3.10.tar.gz

../static-libs/lib/libfftw3f.a: ../fftw-3.3.10.tar.gz
	rm -f $(PWD)/../static-libs/lib/*fftw*
	cd .. && tar -xzf fftw-3.3.10.tar.gz && mv fftw-3.3.10 fftw
	cd ../fftw && MACOSX_DEPLOYMENT_TARGET=11.0 ./configure --enable-single --prefix $(PWD)/../static-libs CFLAGS="-arch arm64 -arch x86_64 -mmacosx-version-min=11.0" && MACOSX_DEPLOYMENT_TARGET=11.0 make clean && MACOSX_DEPLOYMENT_TARGET=11.0 make -j4 install
	rm -Rf ../fftw

libfftw: ../static-libs/lib/libfftw3f.a


../wxWidgets-3.2.5.tar.bz2:
	cd .. && $(DOWNLOAD_COMMAND)  https://github.com/wxWidgets/wxWidgets/releases/download/v3.2.5/wxWidgets-3.2.5.tar.bz2

../static-libs/bin/wx-config: ../wxWidgets-3.2.5.tar.bz2
	cd .. && tar -xjf wxWidgets-3.2.5.tar.bz2
	cd ../wxWidgets-3.2.5 && MACOSX_DEPLOYMENT_TARGET=11.0 ./configure --with-opengl --disable-debug --disable-webviewwebkit --enable-optimise --disable-shared  --enable-monolithic --disable-mediactrl --without-libtiff --enable-cxx11 --enable-universal_binary=x86_64,arm64 --prefix=$(PWD)/../static-libs && make -j4 install
	rm -Rf ../wxWidgets-3.2.5

libwx: ../static-libs/bin/wx-config

default: Fox

Fox: libfftw libwx ../cctbx ../newmat
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Deployment

Fox-nogui: ../cctbx ../newmat
	make -f gnu.mak Fox-nogui

dist:Fox
	rm -Rf Fox-`date "+%Y-%m-%d"` Fox-`date "+%Y-%m-%d"`.dmg
	mkdir Fox-`date "+%Y-%m-%d"`
	cp -R build/Deployment/Fox.app example ../LICENSE.txt ../ChangeLog.txt ../README.rst README-Fox.txt Fox-`date "+%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-`date "+%Y-%m-%d"` Fox-`date "+%Y-%m-%d"`.dmg
	rm -Rf Fox-`date "+%Y-%m-%d"`

all: Fox Fox-nogui

tidy:
	make -f gnu.mak tidy

clean:
	rm -Rf build/Fox.build
	make -f gnu.mak clean
