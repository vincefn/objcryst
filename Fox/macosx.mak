default: all

all: Aqua Aqua-G4 NoGUI NoGUI-G4 

Aqua:
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Aqua

Aqua-G4:
	xcodebuild -project Fox.xcodeproj -target Fox -configuration Aqua-G4

Fox-nogui: NoGUI

NoGUI:
	xcodebuild -project Fox.xcodeproj -target NoGUI -configuration NoGUI

Aqua-dist:Aqua
	rm -Rf Fox`date "+-%Y-%m-%d"` Fox`date "+-%Y-%m-%d"`.dmg
	mkdir Fox`date "+-%Y-%m-%d"`
	cp -R build/Aqua/Fox.app example Fox`date "+-%Y-%m-%d"`/
	hdiutil create -srcfolder Fox`date "+-%Y-%m-%d"` Fox`date "+-%Y-%m-%d"`.dmg
	rm -Rf Fox`date "+-%Y-%m-%d"`

Aqua-G4-dist:Aqua-G4
	rm -Rf Fox-G4`date "+-%Y-%m-%d"` Fox-G4`date "+-%Y-%m-%d"`.dmg
	mkdir Fox-G4`date "+-%Y-%m-%d"`
	cp -R build/Aqua-G4/Fox.app example Fox-G4`date "+-%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-G4`date "+-%Y-%m-%d"` Fox-G4`date "+-%Y-%m-%d"`.dmg
	rm -Rf Fox-G4`date "+-%Y-%m-%d"`

clean:
	$(MAKE) -f macosx.mak -C src clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} clean

update:
	cvs -z3 update
	cd ${DIR_CRYST} ; cvs -z3 update
   
cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
	$(MAKE) -f gnu.mak -C src cvsignore
	$(MAKE) -f gnu.mak -C src-doc cvsignore
