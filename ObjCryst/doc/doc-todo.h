/*! \page page_wishlist Wishlist/Todo page for ObjCryst++
*
* (stars indicate priority/complexity : (*****,**) is urgent,fairly easy)
*
*\section external Outside features
*\par Examples
*     Make a few examples for each feature...
*\par Compilation
*     Make the library and test programs compilable for other OS&configuration...
*     So far tests have been made on Linux (x86/ppc), MacOS and windows (cygwin)
*     Need to work on autoconf and automake (...)
*
*\section general Library Features
*\par Check Results (****,*) (in constant progress)
*     Compare Structure Factors with the results of other programs (Jana,...)
*     for a number of different spacegroups.
*\par Save & restore objects (******,***) DONE, using XML
*\par Move from atominfo+sglite to cctbx (******,***)
*  Change to the newest R. Grosse-Kunstleve library.
*\par CIF Import & Export (****,*****)
*     Ability to export (and import) Crystallographic Info Files. The importing will
*     be much harder, and is not a priority.
*\par Export to other LSQ refinement programs(****,**)
*     Ability to export for other refinement programs (fullprof, gsas, shellx,...).
*\par Anisotropic Thermic Factors (**,****)
*     Add support for anisotropic thermic factors. Spacegroup object should be able
*     to indicate the permutations needed for symetrical atoms. Also determine
*     the constraints between bij (eg beta12=beta13,etc...)
*     The Interface has been written, but no code.
*\par Powder Diffraction Background(**,**)
*     Use splines to interpolate background. Automagically determine background by filtering
*     the spectrum. So far only linearly-interpolated background is available.
*\par Intensity extraction from powder (*,***)
*     Extract F(hkl) from a powder spectrum. Not so interesting, there already are
*      numerous programs to do this. But the advantage to have it embeded is the possibility
*     to refine and refresh (recycle) these extracted intensities during an optimization.
*\par Multi-phase (****,**) - DONE
*     Multiple phase for powder diffraction.
*\par ZScatterer import(**,**)
*     Import Z-Matrix from file. Add the possibility to link two existing ZScatterer by
*     linking terminal atoms.
*\par Genetic Algorithm (****,***)
*      And compare the 3 different algorithms on several structures.
*\par Texture (***,***)
*      Include texture parameters in PowderPattern. Optionnaly, allow multiple
*      textured patterns (that would more difficult).
*
*\section internal Internal Design of the Library
*\par Vector & arrays (Blitz++ usage,...) (**,***)
*     Test again the use of the blitz++ library for vector computing.
*
*\section wild Wild ideas (interesting possibilities, but no priority given)
*\par Amorphous contribution
*     Add amorphous (disordered) contributions to powder spectrum. This could be easily
*     added from the PowderPatternComponent base class.
*\par Magnetic scattering
*     Add calculations for magnetic scattering (->scattering factors)
*/
