/*! \page page_install ObjCryst++ Download & Installation Notes
*
*\section download Download
*     You can download the ObjCryst library on the 
*     <a href="http://www.sourceforge.net/projects/objcryst/">ObjCryst project page</a>
*     on SourceForge. You will need to download the ObjCryst tar.gz, as well as
*     the 'external libraries' tar.gz archives (the latter include libraries which are
*     used by this project, but which are not part of it). Simply uncompress both archives
*     in the same directory. This will create the atominfo, sglite and newmat directories,
*     as well as the ObjCryst directory. This documentation is in ObjCryst/doc/html/
*
*\section compile Compiling
*     
*\par ANSI C++ Compliance
*     The ObjCryst++ library makes use of some (relatively) recent C++ features, such as
*     templates, exceptions, mutable members. This implies that you need a compiler with
*     reasonable C++ compliance in order to use this library.
*\par Makefiles
*     Under linux/unix, you have to edit the rules.mak file. Note that the rules in this
*     makefiles follow the GNU make format...Some other make programs may choke on that..
*     You may want to remove the 'profile' and 'debug' sections.
*
*\section platform Platforms
*\par Linux
*     The library is developped on a Linux computer, using the <a href="http://gcc.gnu.org">
*     GNU compiler gcc </a>, version 2.95.2 or more recent.
*\par Windows
*     The library has been tested under windows (NT) using <a href="http://www.cygwin.com">
*     cygwin's port of the GNU gcc compiler </a>. The free cygwin port creates a complete
*     UNIX environment within the windows OS. Other compiler should work, provided that
*     they support enough of the C++ standard (...).
*\par Macintosh
*     The library can also be compiled on MacOS, using <a href="http://www.metrowerks.com">
*     Metrowerks Codewarrior Pro 6 </a> (version 5 should work, too). Note that the Blitz++
*     may not compile in its current form with Codewarrior (untested), but Blitz++ is not
*     used by the library currently, so... I don't provide a project file yet.
*
*\section example Example
*     There is a short example in the ObjCryst/example directory. Under linux/unix,
*     (it works the same way under windows using the cygwin port of gcc)
*     change to the directory and simply type 'make all'. That should compile
*     all necessary libraries. Then run test-pbso4.exe, it should find the PbSO4 structure
*     (run it several times if necessary- it's fast).
*     from the experimental neutron powder spectrum used during the Rietveld Round Robin.
*     You can use the output crystal.pov file with <a href="http://www.povray.org/">
*     povray</a> to see the structure. The calculated powder spectrum is stored
*     in pbso4-xray/calc.out.
*     You can also try 'pbso4-neutron' (on neutron powder spectrum) and 
*     'make pbso4-xray2' (on extracted F's).
*     \subsection example_result
*     On a 500 Mhz linux box, pbso4-xray2 (using only the first 80 reflections and no
*     dynamical occupancy correction),
*     50000 trials are done in less than 30 s (for 5 independent atoms, each with 4
*     symetrics (not counting the center of symetry): that's 
*     80*50000*5*4/27.2 = 2.9 10^6 reflection.atom per second).
*     For pbso4-xray and pbso4-neutron, it takes 2 to 3 more time, using only the first
*     90 degrees of the powder spectrum (about 80 reflections again), because a full
*     powder spectrum is generated (applying profiles is slow), and dynamical occupancy
*     correction is used to take special positions into account (computing distances
*     takes a lot of time). On the final structure output, note that the population consists
*     of two factor, the first (1.0 here) is the 'real' occupancy of the site, and the
*     second is the dynamical correction due to the fact that several identical atoms
*     overlap (For PbSo4, all atoms should be at 0.5 dyn corr).
*/
