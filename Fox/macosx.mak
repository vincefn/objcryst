DIR_CRYST := ObjCryst

default: all

all: Fox

Fox:
	$(MAKE) -f macosx.mak wxcryst=1 opengl=1 debug=$(debug) -C src all

doc:
	$(MAKE) -f gnu.mak -C src-doc all

clean:
	$(MAKE) -f macosx.mak -C src clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} clean

tidy:
	$(MAKE) -f macosx.mak -C src tidy
	$(MAKE) -f gnu.mak -C ${DIR_CRYST} tidy

#install Fox in /usr/local/bin
install:
	install -m 755 src/Fox /usr/local/bin

update:
	cvs -z3 update
	cd ${DIR_CRYST} ; cvs -z3 update
   
cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
	$(MAKE) -f gnu.mak -C src cvsignore
	$(MAKE) -f gnu.mak -C src-doc cvsignore
