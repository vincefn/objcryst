DIR_CRYST := ObjCryst

default: all

all: Fox doc

Fox:
	$(MAKE) -f gnu.mak wxcryst=1 opengl=1 debug=$(debug) -C src all

doc:
	$(MAKE) -f gnu.mak -C src-doc all

clean:
	$(MAKE) -f gnu.mak -C src clean

tidy:
	$(MAKE) -f gnu.mak -C src tidy

#install Fox in /usr/local/bin
install:
	install -m 755 Fox /usr/local/bin

update:
	cvs -z3 update
	cd ${DIR_CRYST} ; cvs -z3 update
   
cvsignore:
	cp -f ${DIR_CRYST}/.cvsignore ./
	$(MAKE) -f gnu.mak -C src cvsignore
	$(MAKE) -f gnu.mak -C src-doc cvsignore
