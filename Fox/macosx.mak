default: all

all: Aqua Aqua-G4 NoGUI NoGUI-G4 

Aqua-G3:
	xcodebuild -project Fox.xcodeproj -target Fox -configuration G3 

Aqua-G4:
	xcodebuild -project Fox.xcodeproj -target Fox -configuration G4

Aqua-G5:
	xcodebuild -project Fox.xcodeproj -target Fox -configuration G5

Fox-nogui: NoGUI-G3

NoGUI-G3:
	xcodebuild -project Fox.xcodeproj -target NoGUI -configuration G3

NoGUI-G4:
	xcodebuild -project Fox.xcodeproj -target NoGUI -configuration G4

NoGUI-G5:
	xcodebuild -project Fox.xcodeproj -target NoGUI -configuration G5

G3: Aqua-G3 NoGUI-G3

G4: Aqua-G4 NoGUI-G4

G5: Aqua-G5 NoGUI-G5

Aqua-G3-dist:G3
	rm -Rf Fox-G3-`date "+-%Y-%m-%d"` Fox`date "+-%Y-%m-%d"`.dmg
	mkdir Fox-G3-`date "+-%Y-%m-%d"`
	cp -R build/G3/Fox.app build/G3/Fox-NoGUI example Fox-G3-`date "+-%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-G3-`date "+-%Y-%m-%d"` Fox-G3-`date "+-%Y-%m-%d"`.dmg
	rm -Rf Fox-G3-`date "+-%Y-%m-%d"`

Aqua-G4-dist:Aqua-G4
	rm -Rf Fox-G4-`date "+-%Y-%m-%d"` Fox-G4-`date "+-%Y-%m-%d"`.dmg
	mkdir Fox-G4-`date "+-%Y-%m-%d"`
	cp -R build/G4/Fox.app build/G4/Fox-NoGUI example Fox-G4-`date "+-%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-G4-`date "+-%Y-%m-%d"` Fox-G4-`date "+-%Y-%m-%d"`.dmg
	rm -Rf Fox-G4-`date "+-%Y-%m-%d"`

Aqua-G5-dist:Aqua-G5
	rm -Rf Fox-G5-`date "+-%Y-%m-%d"` Fox-G5-`date "+-%Y-%m-%d"`.dmg
	mkdir Fox-G5-`date "+-%Y-%m-%d"`
	cp -R build/G5/Fox.app build/G5/Fox-NoGUI example Fox-G5-`date "+-%Y-%m-%d"`/
	hdiutil create -srcfolder Fox-G5-`date "+-%Y-%m-%d"` Fox-G5-`date "+-%Y-%m-%d"`.dmg
	rm -Rf Fox-G5-`date "+-%Y-%m-%d"`

clean:
	rm -Rf build/Fox.build

update:
	cvs -z3 update
	cd ${DIR_CRYST} ; cvs -z3 update
   
cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
	$(MAKE) -f gnu.mak -C src cvsignore
	$(MAKE) -f gnu.mak -C src-doc cvsignore

#Switch to ssh developer access
cvs-ext:
	perl -pi -e 's|pserver|ext|g' `find . -name Root`
	perl -pi -e 's|anonymous|vincefn|g' `find . -name Root`

#Switch to anonymous CVS
cvs-anon:
	perl -pi -e 's|ext|pserver|g' `find . -name Root`
	perl -pi -e 's|vincefn|anonymous|g' `find . -name Root`

