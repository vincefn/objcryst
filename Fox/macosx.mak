../static-libs/lib/libfftw3f.a:
	cd .. && tar -xjf fftw.tar.bz2
	cd ../fftw && ./configure --enable-single --prefix ../static-libs && make install
	rm -Rf ../fftw


libfftw: ../static-libs/lib/libfftw3f.a

../static-libs/bin/wx-config:
	cd .. && tar -xjf wxMac-2.8.7.tar.bz2
	cd ../wxMac-2.8.7 && ./configure --with-opengl --enable-optimise --disable-shared --enable-monolithic --prefix=$HOME/dev/Fox/static-libs && make install


libwx: ../static-libs/bin/wx-config

default: Fox

Fox: libfftw libwx
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Deployment

Fox-nogui:
	make -f gnu.mak Fox-nogui

dist:Fox Fox-nogui
	rm -Rf Fox-`arch`-`date "+%Y-%m-%d"` Fox-`arch`-`date "+%Y-%m-%d"`.dmg
	mkdir Fox-`arch`-`date "+%Y-%m-%d"`
	cp -R build/Deployment/Fox.app example Fox-`arch`-`date "+%Y-%m-%d"`/
	cp src/Fox Fox-`arch`-`date "+%Y-%m-%d"`/Fox-noGUI
	hdiutil create -srcfolder Fox-`arch`-`date "+%Y-%m-%d"` Fox-`arch`-`date "+%Y-%m-%d"`.dmg
	rm -Rf Fox-`arch`-`date "+%Y-%m-%d"`

all: Fox Fox-nogui

clean:
	rm -Rf build/Fox.build
	make -f gnu.mak clean

