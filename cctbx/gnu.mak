# Which command to use for download ?
CURL=$(shell which curl 2>/dev/null)
ifneq ($(CURL),)
DOWNLOAD_COMMAND=curl -L -O
else
DOWNLOAD_COMMAND=wget
endif

../cctbx_bundle.tar.gz:
	cd .. && $(DOWNLOAD_COMMAND) http://cci.lbl.gov/cctbx_build/results/last_published/cctbx_bundle.tar.gz

#This target can be used to re-download cctbx from source
../cctbx.tar.bz2:../cctbx_bundle.tar.gz
	cd .. && tar -xjf cctbx_bundle.tar.gz && csh cctbx_install_script.csh 4
	rm -Rf cctbx scitbx include *.a
	rsync -ar ../cctbx_sources/cctbx --filter="+ */" --filter="- adp_restraints/" --filter="- adptbx/" --filter="- boost_python/" --filter="- command_line/" --filter="- crystal/" --filter="- development/" --filter="- dmtbx/" --filter="- examples/" --filter="- geometry_restraints/" --filter="- macro_mol/" --filter="- maptbx/" --filter="- neutron/" --filter="- reference/" --filter="- regression/" --filter="- source_generators/" --filter="- translation_search/" --filter="- web/" --filter="+ *.cpp" --filter="- *" ./
	cp ../cctbx_build/cctbx/eltbx/*.cpp cctbx/eltbx/
	mkdir include
	rsync -ar ../cctbx_sources/cctbx --filter="+ */" --filter="- adp_restraints/" --filter="- adptbx/" --filter="- boost_python/" --filter="- command_line/" --filter="- crystal/" --filter="- development/" --filter="- dmtbx/" --filter="- examples/" --filter="- geometry_restraints/" --filter="- macro_mol/" --filter="- maptbx/" --filter="- neutron/" --filter="- reference/" --filter="- regression/" --filter="- source_generators/" --filter="- translation_search/" --filter="- web/" --filter="+ *.h" --filter="- *" include/
	mkdir -p scitbx/include/
	rsync -ar ../cctbx_sources/scitbx --filter="+ */" --filter="- boost_python/" --filter="+ *.h" --filter="- *" include/
	cp -Rf ../cctbx_build/include/scitbx include/
	# Copy boost from the cctbx bundle
	rsync -ar ../cctbx_sources/boost/boost --exclude="accumulators/" --exclude="archive/" --exclude="asio/" --exclude="asign/" --exclude="bimap/" --exclude="bind/" --exclude="circular_buffer/" --exclude="concept_check/" --exclude="dynamic_bitset/" --exclude="filesystem/" --exclude="flyweight/" --exclude="function/" --exclude="function_types/"  --exclude="gil/" --exclude="graph/" --exclude="interprocess/" --exclude="intrusive/" --exclude="iostreams/" --exclude="lambda/" --exclude="logic/" --exclude="mpi/" --exclude="parameter/" --exclude="pending/" --exclude="pool/" --exclude="program_options/" --exclude="property_tree/" --exclude="proto/" --exclude="ptr_container/" --exclude="python/" --exclude="random/" --exclude="regex/" --exclude="serialization/" --exclude="signals/" --exclude="spirit/" --exclude="statechart/" --exclude="system/" --exclude="test/" --exclude="thread/" --exclude="units/" --exclude="unordered/" --exclude="wave/" --exclude="xpressive/" --filter="+ */" --filter="+ *.hpp"  --filter="+ *.ipp" include/
	cd .. && tar -cjf cctbx.tar.bz2 cctbx
	rm -Rf ../cctbx_install_script.csh ../cctbx_sources ../cctbx_build


#Boost library - if we want to use a more recent version from boost.org
../boost_1_58_0.tar.bz2:
	cd .. && $(DOWNLOAD_COMMAND)  http://downloads.sourceforge.net/project/boost/boost/1.58.0/boost_1_58_0.tar.bz2

libboost-orig:../boost_1_58_0.tar.bz2
	cd .. && tar -xjf boost_1_58_0.tar.bz2
	rm -Rf include/boost
	cd ../boost_1_58_0 && rsync -ar boost --exclude="accumulators/" --exclude="archive/" --exclude="asio/" --exclude="asign/" --exclude="bimap/" --exclude="bind/" --exclude="circular_buffer/" --exclude="concept_check/" --exclude="dynamic_bitset/" --exclude="filesystem/" --exclude="flyweight/" --exclude="function/" --exclude="function_types/"  --exclude="gil/" --exclude="graph/" --exclude="interprocess/" --exclude="intrusive/" --exclude="iostreams/" --exclude="lambda/" --exclude="logic/" --exclude="mpi/" --exclude="parameter/" --exclude="pending/" --exclude="pool/" --exclude="program_options/" --exclude="property_tree/" --exclude="proto/" --exclude="ptr_container/" --exclude="python/" --exclude="random/" --exclude="regex/" --exclude="serialization/" --exclude="signals/" --exclude="spirit/" --exclude="statechart/" --exclude="system/" --exclude="test/" --exclude="thread/" --exclude="units/" --exclude="unordered/" --exclude="wave/" --exclude="xpressive/" --filter="+ */" --filter="+ *.hpp"  --filter="+ *.ipp" $(PWD)/include/
	rm -Rf ../boost_1_58_0



CXXFLAGS = -O3 -Wall -pedantic -Iinclude -I../static-libs/include/ -Wno-long-long

OBJ_ELTBX= cctbx/eltbx/basic.o cctbx/eltbx/icsd_radii.o cctbx/eltbx/tiny_pse.o cctbx/eltbx/fp_fdp.o \
           cctbx/eltbx/neutron.o cctbx/eltbx/wavelengths.o \
           cctbx/eltbx/xray_scattering/it1992.o cctbx/eltbx/xray_scattering/n_gaussian.o \
           cctbx/eltbx/xray_scattering/n_gaussian_raw.o cctbx/eltbx/xray_scattering/wk1995.o \
           cctbx/eltbx/henke.o cctbx/eltbx/henke_tables_25_36.o  cctbx/eltbx/henke_tables_61_72.o  \
           cctbx/eltbx/henke_tables_01_12.o cctbx/eltbx/henke_tables_37_48.o cctbx/eltbx/henke_tables_73_84.o  \
           cctbx/eltbx/henke_tables_13_24.o cctbx/eltbx/henke_tables_49_60.o cctbx/eltbx/henke_tables_85_92.o \
           cctbx/eltbx/sasaki.o cctbx/eltbx/sasaki_tables_25_36.o cctbx/eltbx/sasaki_tables_61_72.o \
           cctbx/eltbx/sasaki_tables_01_12.o cctbx/eltbx/sasaki_tables_37_48.o cctbx/eltbx/sasaki_tables_73_82.o \
           cctbx/eltbx/sasaki_tables_13_24.o cctbx/eltbx/sasaki_tables_49_60.o


OBJ_SGTBX = cctbx/sgtbx/bricks.o cctbx/sgtbx/miller.o cctbx/sgtbx/select_generators.o \
            cctbx/sgtbx/tr_group.o cctbx/sgtbx/change_of_basis_op.o cctbx/sgtbx/reciprocal_space_asu.o \
            cctbx/sgtbx/seminvariant.o cctbx/sgtbx/tr_vec.o \
            cctbx/sgtbx/find_affine.o cctbx/sgtbx/reciprocal_space_ref_asu.o cctbx/sgtbx/site_symmetry.o cctbx/sgtbx/utils.o \
            cctbx/sgtbx/group_codes.o cctbx/sgtbx/rot_mx.o  cctbx/sgtbx/space_group.o cctbx/sgtbx/wyckoff.o \
            cctbx/sgtbx/hall_in.o cctbx/sgtbx/rot_mx_info.o cctbx/sgtbx/space_group_type.o \
            cctbx/sgtbx/lattice_symmetry.o cctbx/sgtbx/row_echelon_solve.o cctbx/sgtbx/symbols.o \
            cctbx/sgtbx/lattice_tr.o cctbx/sgtbx/rt_mx.o  cctbx/sgtbx/tensor_rank_2.o \
            cctbx/sgtbx/reference_settings/hall_symbol_table.o cctbx/sgtbx/reference_settings/normalizer.o \
            cctbx/sgtbx/reference_settings/matrix_group_code_table.o cctbx/sgtbx/reference_settings/wyckoff.o


OBJ_MILLER = cctbx/miller/asu.o cctbx/miller/index_generator.o cctbx/miller/match_bijvoet_mates.o \
             cctbx/miller/sym_equiv.o cctbx/miller/bins.o cctbx/miller/index_span.o cctbx/miller/match_indices.o

OBJ_UCTBX = cctbx/uctbx/uctbx.o cctbx/uctbx/spoil_optimization.o


OBJ = $(OBJ_ELTBX) $(OBJ_SGTBX) $(OBJ_MILLER) $(OBJ_UCTBX)

MAKEDEPEND = gcc -MM ${CXXFLAGS} $< > $*.dep

%.o : %.cpp
	@$(MAKEDEPEND)
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@

libcctbx.a : ${OBJ}
	@${RM} $@
	${AR} crs $@ ${filter-out %.a %.so, $^}

-include $(OBJ:.o=.dep)

lib: libcctbx.a

install:lib
	cp -f libcctbx.a ../static-libs/lib/
	rm -Rf ../static-libs/include/cctbx  ../static-libs/include/scitbx ../static-libs/include/boost
	cp -Rf include/* ../static-libs/include/

# target for removing all object files
.PHONY : tidy
tidy::
	@${RM} -f core ${OBJ} *.dep */*.dep */*/*.dep

# target for removing all object files and libraries
.PHONY : clean
clean:: tidy
	@${RM} *.a

#scitbx.constants.h
