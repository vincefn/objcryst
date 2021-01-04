#!/bin/bash

# This script will create the cctbx and scitbx sub-directory from an installed cctbx in /usr/local/cctbx/,
#including the required boost headers


rm -Rf cctbx scitbx include *.a

rsync -ar /usr/local/cctbx/cctbx_sources/cctbx --filter="+ */" --filter="- adp_restraints/" --filter="- adptbx/" --filter="- boost_python/" --filter="- command_line/" --filter="- crystal/" --filter="- development/" --filter="- dmtbx/" --filter="- examples/" --filter="- geometry_restraints/" --filter="- macro_mol/" --filter="- maptbx/" --filter="- neutron/" --filter="- reference/" --filter="- regression/" --filter="- source_generators/" --filter="- translation_search/" --filter="- web/" --filter="+ *.cpp" --filter="- *" ./

cp /usr/local/cctbx/cctbx_build/cctbx/eltbx/*.cpp cctbx/eltbx/

mkdir include
rsync -ar /usr/local/cctbx/cctbx_sources/cctbx --filter="+ */" --filter="- adp_restraints/" --filter="- adptbx/" --filter="- boost_python/" --filter="- command_line/" --filter="- crystal/" --filter="- development/" --filter="- dmtbx/" --filter="- examples/" --filter="- geometry_restraints/" --filter="- macro_mol/" --filter="- maptbx/" --filter="- neutron/" --filter="- reference/" --filter="- regression/" --filter="- source_generators/" --filter="- translation_search/" --filter="- web/" --filter="+ *.h" --filter="- *" include/

mkdir -p scitbx/include/
rsync -ar /usr/local/cctbx/cctbx_sources/scitbx --filter="+ */" --filter="- boost_python/" --filter="+ *.h" --filter="- *" include/
cp -Rf /usr/local/cctbx/cctbx_build/include/scitbx include/

cp /usr/local/cctbx/cctbx_sources/cctbx/LICENSE_2_0.txt CCTBX_LICENSE_2_0.txt


# Now copy the boost parts cctbx needs
rsync -ar /usr/local/cctbx/cctbx_sources/boost/boost --exclude="accumulators/" --exclude="asio/" --exclude="asign/" --exclude="bimap/" --exclude="circular_buffer/" --exclude="concept_check/" --exclude="dynamic_bitset/" --exclude="filesystem/" --exclude="flyweight/" --exclude="function_types/" --exclude="gil/" --exclude="graph/" --exclude="interprocess/" --exclude="intrusive/" --exclude="iostreams/" --exclude="lambda/" --exclude="logic/" --exclude="mpi/" --exclude="parameter/" --exclude="pending/" --exclude="pool/" --exclude="program_options/" --exclude="proto/" --exclude="ptr_container/" --exclude="python/" --exclude="random/" --exclude="regex/" --exclude="signals/" --exclude="spirit/" --exclude="statechart/" --exclude="system/" --exclude="test/" --exclude="thread/" --exclude="units/" --exclude="unordered/" --exclude="wave/" --exclude="xpressive/" --filter="+ */" --filter="+ *.hpp"  --filter="+ *.ipp" --filter="- *" include/

cp /usr/local/cctbx/cctbx_sources/boost/LICENSE_1_0.txt BOOST_LICENSE.txt

cd .. ; tar -cjf cctbx.tar.bz2 cctbx ; cd cctbx
