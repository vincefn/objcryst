../cctbx:
	cd .. && tar -xjf cctbx.tar.bz2

../newmat:
	cd .. && tar -xjf newmat.tar.bz2

../static-libs/lib/libfftw3f.a:
	cd .. && tar -xjf fftw.tar.bz2
	cd ../fftw && ./configure --enable-single CFLAGS="-arch i386" --prefix $(PWD)/../static-libs && make clean && make install && mv .libs/libfftw3f.a ./libfftw3f-i386.a
	cd ../fftw && ./configure --enable-single CFLAGS="-arch ppc" --prefix $(PWD)/../static-libs && make clean && make && mv .libs/libfftw3f.a ./libfftw3f-ppc.a
	rm -f $(PWD)/../static-libs/lib/*fftw*
	cd ../fftw && lipo -create libfftw3f-i386.a libfftw3f-ppc.a -output $(PWD)/../static-libs/lib/libfftw3f.a
	rm -Rf ../fftw


libfftw: ../static-libs/lib/libfftw3f.a

../wxMac-2.8.7.tar.bz2:
	cd .. && curl -O http://switch.dl.sourceforge.net/sourceforge/wxwindows/wxMac-2.8.7.tar.bz2

../static-libs/bin/wx-config: ../wxMac-2.8.7.tar.bz2
	cd .. && tar -xjf wxMac-2.8.7.tar.bz2
	cd ../wxMac-2.8.7 && ./configure --with-opengl --enable-optimise --disable-shared --enable-monolithic --enable-universal_binary --prefix=$(PWD)/../static-libs && make install
	rm -Rf ../wxMac-2.8.7

libwx: ../static-libs/bin/wx-config

default: Fox

Fox: libfftw libwx ../cctbx ../newmat
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Deployment

Fox-nogui: ../cctbx ../newmat
	make -f gnu.mak Fox-nogui

dist:Fox
	rm -Rf Fox-`arch`-`date "+%Y-%m-%d"` Fox-`arch`-`date "+%Y-%m-%d"`.dmg
	mkdir Fox-`arch`-`date "+%Y-%m-%d"`
	cp -R build/Deployment/Fox.app example Fox-`arch`-`date "+%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-`arch`-`date "+%Y-%m-%d"` Fox-`arch`-`date "+%Y-%m-%d"`.dmg
	rm -Rf Fox-`arch`-`date "+%Y-%m-%d"`

all: Fox Fox-nogui

clean:
	rm -Rf build/Fox.build
	make -f gnu.mak clean

