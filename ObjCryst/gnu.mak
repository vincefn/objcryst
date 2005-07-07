include rules.mak

#local directory is base directory
DIR_CRYST := ./

######################################################################
#####################      LIBRAIRIES         ########################
######################################################################

#All defined in rules.mak

######################################################################
#####################        PROGRAMS         ########################
######################################################################

#test-cristallo
test:
	  cd ${DIR_TEST} && make all

######################################################################
#####################      OTHER TARGETS      ########################
######################################################################

ctags:
	$(MAKE) -f gnu.mak -C ${DIR_DOC} ctags


# target for making both libraries

.PHONY : all
all: libcryst test doc

# target for removing all object files (does not affect blitz/newmat/sglite/atominfo)
.PHONY : tidy
tidy::
	$(MAKE) -f gnu.mak -C ${DIR_LIBCRYST} tidy
	$(MAKE) -f gnu.mak -C ${DIR_CRYSTVECTOR} tidy
	$(MAKE) -f gnu.mak -C ${DIR_VFNQUIRKS} tidy
	$(MAKE) -f gnu.mak -C ${DIR_REFOBJ} tidy
	$(MAKE) -f gnu.mak -C ${DIR_WXWCRYST} tidy
	$(MAKE) -f gnu.mak -C ${DIR_EXAMPLE} tidy

# target for removing all object files and libraries
# (does not affect blitz/newmat/sglite/atominfo)
.PHONY : clean
clean::
	$(MAKE) -f gnu.mak -C ${DIR_LIBCRYST} clean
	$(MAKE) -f gnu.mak -C ${DIR_CRYSTVECTOR} clean
	$(MAKE) -f gnu.mak -C ${DIR_VFNQUIRKS} clean
	$(MAKE) -f gnu.mak -C ${DIR_REFOBJ} clean
	$(MAKE) -f gnu.mak -C ${DIR_WXWCRYST} clean
	$(MAKE) -f gnu.mak -C ${DIR_EXAMPLE} clean


#target to make documentation (requires doxygen)
#also makes tags file, although it is not related to doxygen
doc:
	cd ${DIR_DOC}; $(MAKE) -f gnu.mak doc

#target to make distribution archive of libcryst++
dist:
	tar -czf ../archives/ObjCryst.tar.gz  -C .. --exclude='*.o' --exclude='.systemG.Desktop' --exclude='profile.0.0.0' --exclude='*.a' --exclude='*.pov' --exclude='latex' --exclude='*.exe' --exclude='*.out' --exclude='tags' --exclude='wxCryst/Fox' --exclude='Makefile' ObjCryst --dereference

#these are libraries/programs used by ObjCryst but developped by other people.
#These are needed to use ObjCryst, but not modified-so only get it once.
#
# ObjCryst/blitz
dist-libs:
	tar -czf ../archives/ObjCryst-libs.tar.gz  -C .. --exclude='*.o' --exclude='.systemG.Desktop'  --exclude='*.a' --exclude='*.exe' cctbx newmat --exclude='Makefile' --dereference

#target to make a complete archive of ObjCryst++
archive:
	tar -czf ../archives/ObjCryst-complete.tar.gz  -C .. --exclude='*.o' --exclude='.systemG.Desktop' --exclude='*.a' --exclude='*.exe' --exclude='*.out' --exclude='CVS' --exclude='ObjCryst/doc/html' --exclude='ObjCryst/doc/latex' --exclude='*.oxy' --exclude='profile.0.0.0' --exclude='Makefile' cctbx newmat AsymProfLarryFinger ObjCryst --dereference

#Copy the .cvsignore in all relevant directories
cvsignore:
	$(MAKE) -f gnu.mak -C ${DIR_LIBCRYST} cvsignore
	$(MAKE) -f gnu.mak -C ${DIR_CRYSTVECTOR} cvsignore
	$(MAKE) -f gnu.mak -C ${DIR_VFNQUIRKS} cvsignore
	$(MAKE) -f gnu.mak -C ${DIR_REFOBJ} cvsignore
	$(MAKE) -f gnu.mak -C ${DIR_WXWCRYST} cvsignore
	$(MAKE) -f gnu.mak -C ${DIR_EXAMPLE} cvsignore
	cp -f ${DIR_CRYST}/.cvsignore doc/
