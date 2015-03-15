../cctbx:
	cd .. && tar -xjf cctbx.tar.bz2

../newmat:
	cd .. && tar -xjf newmat.tar.bz2

../fftw-3.3.4.tar.gz:
	cd .. &&  curl -O http://fftw.org/fftw-3.3.4.tar.gz

../static-libs/lib/libfftw3f.a: ../fftw-3.3.4.tar.gz
	rm -f $(PWD)/../static-libs/lib/*fftw*
	cd .. && tar -xzf fftw-3.3.4.tar.gz && mv fftw-3.3.4 fftw
	cd ../fftw && ./configure --enable-single --prefix $(PWD)/../static-libs && make clean && make -j4 install
	#cd ../fftw && ./configure --enable-single CFLAGS="-arch i386" --prefix $(PWD)/../static-libs && make clean && make install && mv .libs/libfftw3f.a ./libfftw3f-i386.a
	#cd ../fftw && ./configure --enable-single CFLAGS="-arch ppc" --prefix $(PWD)/../static-libs && make clean && make && mv .libs/libfftw3f.a ./libfftw3f-ppc.a
	#cd ../fftw && lipo -create libfftw3f-i386.a libfftw3f-ppc.a -output $(PWD)/../static-libs/lib/libfftw3f.a
	rm -Rf ../fftw


libfftw: ../static-libs/lib/libfftw3f.a

../wxWidgets-3.0.2.tar.bz2:
	cd .. && curl -O  ftp://ftp.wxwidgets.org/pub/3.0.2/wxWidgets-3.0.2.tar.bz2

#../static-libs/bin/wx-config: ../wxMac-2.8.12.tar.bz2
#	cd .. && tar -xjf wxMac-2.8.12.tar.bz2
#	cd ../wxMac-2.8.12 && ./configure --with-opengl --enable-optimise --disable-shared --enable-monolithic --enable-universal_binary --prefix=$(PWD)/../static-libs && make install
#	rm -Rf ../wxMac-2.8.12

../static-libs/bin/wx-config: ../wxWidgets-3.0.2.tar.bz2
	cd .. && tar -xjf wxWidgets-3.0.2.tar.bz2
	cd ../wxWidgets-3.0.2 && ./configure --with-opengl --disable-webviewwebkit --enable-optimise --disable-shared --enable-debug --prefix=$(PWD)/../static-libs && make -j4 install
	rm -Rf ../wxWidgets-3.0.2

libwx: ../static-libs/bin/wx-config

default: Fox

Fox: libfftw libwx ../cctbx ../newmat
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Deployment

Fox-nogui: ../cctbx ../newmat
	make -f gnu.mak Fox-nogui

dist:Fox
	rm -Rf Fox-Universal-`date "+%Y-%m-%d"` Fox-Universal-`date "+%Y-%m-%d"`.dmg
	mkdir Fox-Universal-`date "+%Y-%m-%d"`
	cp -R build/Deployment/Fox.app example Fox-Universal-`date "+%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-Universal-`date "+%Y-%m-%d"` Fox-Universal-`date "+%Y-%m-%d"`.dmg
	rm -Rf Fox-Universal-`date "+%Y-%m-%d"`

all: Fox Fox-nogui

clean:
	rm -Rf build/Fox.build
	make -f gnu.mak clean

