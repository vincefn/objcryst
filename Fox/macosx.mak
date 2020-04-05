# Which command to use for download ?
CURL=$(shell which curl 2>/dev/null)
ifneq ($(CURL),)
DOWNLOAD_COMMAND=curl -L -O
else
DOWNLOAD_COMMAND=wget
endif

../cctbx:
	cd .. && tar -xjf cctbx.tar.bz2

../newmat:
	cd .. && tar -xjf newmat.tar.bz2

../fftw-3.3.4.tar.gz:
	cd .. &&  curl -O http://fftw.org/fftw-3.3.4.tar.gz

../static-libs/lib/libfftw3f.a: ../fftw-3.3.4.tar.gz
	rm -f $(PWD)/../static-libs/lib/*fftw*
	cd .. && tar -xzf fftw-3.3.4.tar.gz && mv fftw-3.3.4 fftw
	cd ../fftw && MACOSX_DEPLOYMENT_TARGET=10.5 ./configure --enable-single --prefix $(PWD)/../static-libs && MACOSX_DEPLOYMENT_TARGET=10.5 make clean && MACOSX_DEPLOYMENT_TARGET=10.5 make -j4 install
	rm -Rf ../fftw

libfftw: ../static-libs/lib/libfftw3f.a

# MySQL
../mysql-5.6.24.tar.gz:
	cd .. && $(DOWNLOAD_COMMAND) http://dev.mysql.com/get/Downloads/MySQL-5.6/mysql-5.6.24.tar.gz

#:TODO: find a way to only compile the static version of libmysqlclient ?
../static-libs/lib/libmysqlclient.a: ../mysql-5.6.24.tar.gz
	cd .. && tar -xzf mysql-5.6.24.tar.gz
	cd ../mysql-5.6.24 && MACOSX_DEPLOYMENT_TARGET=10.5 cmake -DCMAKE_INSTALL_PREFIX=$(PWD)/../static-libs && MACOSX_DEPLOYMENT_TARGET=10.5 $(MAKE) -j4 install
	rm -f $(PWD)/../static-libs/lib/libmysql*.dylib
	rm -Rf ../mysql-5.6.24

libmysql: ../static-libs/lib/libmysqlclient.a


../wxWidgets-3.0.4.tar.bz2:
	cd .. && $(DOWNLOAD_COMMAND)  https://github.com/wxWidgets/wxWidgets/releases/download/v3.0.4/wxWidgets-3.0.4.tar.bz2

../static-libs/bin/wx-config: ../wxWidgets-3.0.4.tar.bz2
	cd .. && tar -xjf wxWidgets-3.0.4.tar.bz2
	cd ../wxWidgets-3.0.4 && ./configure --with-opengl --disable-webviewwebkit --enable-optimise --disable-shared  --enable-monolithic --disable-mediactrl --prefix=$(PWD)/../static-libs && make -j4 install
	rm -Rf ../wxWidgets-3.0.4

libwx: ../static-libs/bin/wx-config

default: Fox

Fox: libfftw libwx ../cctbx ../newmat libmysql
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

clean:
	rm -Rf build/Fox.build
	make -f gnu.mak clean
