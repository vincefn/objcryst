!include rules.mak

# target for removing all object files (does not affect blitz/newmat/sglite/atominfo)
.PHONY : tidyall
tidyall:
   cd ${DIR_LIBCRYST}
	$(MAKE) -f bc32.mak tidy
   cd ${DIR_CRYSTVECTOR}
	$(MAKE) -f bc32.mak tidy
   cd ${DIR_VFNQUIRKS}
	$(MAKE) -f bc32.mak tidy
   cd ${DIR_REFOBJ}
	$(MAKE) -f bc32.mak tidy
   cd ${DIR_WXWCRYST}
	$(MAKE) -f bc32.mak tidy

# target for removing all object files and libraries
# (does not affect blitz/newmat/sglite/atominfo)
.PHONY : cleanall
cleanall:
   cd ${DIR_LIBCRYST}
	$(MAKE) -f bc32.mak clean
   cd ${DIR_CRYSTVECTOR}
	$(MAKE) -f bc32.mak clean
   cd ${DIR_VFNQUIRKS}
	$(MAKE) -f bc32.mak clean
   cd ${DIR_REFOBJ}
	$(MAKE) -f bc32.mak clean
   cd ${DIR_WXWCRYST}
	$(MAKE) -f bc32.mak clean
