/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: vincefn@users.sourceforge.net
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*  libcryst-doc.h 
*  This file includes only documentation, to be formatted automatically by doxygen.
*/

/*! \mainpage ObjCryst++ Object-Oriented Crystallographic Library
*
*\section index Index
*  \li <a href="http://objcryst.sourceforge.net/">ObjCryst++ Homepage</a>
*  \li <a href="http://www.sourceforge.net/projects/objcryst/">Project page on SourceForge</a>
*  \li \ref page_install
*  \li \ref page_status
*  \li \ref page_wishlist
*  \li \ref page_legal
*  \li \ref page_oop
*  \li \ref page_design
*
*\section news News
*  A new version of the Library should be available in the second half of September,2001.
* It will include major changes for several classes (ScatteringPower, PowderPattern),
* and will include the graphical classes associated with all objects using
* <a href="http://www.wxWindows.org">wxWindows</a>, and saving all objects using a new,
* more expandable <a href="http://www.w3.org/XML/">XML</a>-based format.
*
*\section intro Introduction
*\par Overview of the project
*
* The aim of this project is to provide object-oriented tools for Crystallography. The
* focus of these tools is the analysis of diffraction data, with an emphasis
* on speed optimization for repetitive calculations (for global refinement methods such
* as simulated annealing, for example).
*\par 
* Even if we intend to use this library mainly for the development of a global optimization
* program from powder diffraction, this library is programmed in a general way so that
* other applications can make use of it.
* The library was designed sto be reusable, by adding new kind of experiments, new algorithms,
* new Scatterer type, new ScatteringPower,... See \ref page_design to learn more about
* the object-oriented design of this library and why it is good for its expandability.
* 
*\par Contributors
* Programmation: Vincent Favre-Nicolin (http://www.unige.ch/crystal/favrev/),
* <a href="mailto:vincefn@users.sourceforge.net">vincefn@users.sourceforge.net</a>
*\par 
* This project is being done in the laboratory of Crystallography of the University
* of Geneva (http://www.unige.ch/crystal/), and is part of the development of a global
* optimization program with Radovan Cerny (http://www.unige.ch/crystal/cerny/rcerny.htm).
* \par
* This project is supported by the <a href="http://www.snf.ch/">Swiss
* National Science Foundation</a>.
* 
*\par 
* This project also makes use of some other programs or libraries. Most notable in this 
* project are the spacegroup (SgLite) and atominfo packages from 
* <a href="mailto:rwgk@cci.lbl.gov">Ralf Grosse-Kunstleve</a>
* (see \ref page_legal).
*
*\par 
* This page & project is hosted on
* <A href="http://sourceforge.net"> <IMG src="http://sourceforge.net/sflogo.php?group_id=27546" * width="88" height="31" border="0" alt="SourceForge Logo"></A>
*
*/

/*! \page page_legal License Information
*\par Copyright
*  The ObjCryst++ library is copyright (2000-2001) 
* <a href="mailto:vincefn@users.sourceforge.net">Vincent FAVRE-NICOLIN</a>
*\par License
* This is <a href="http://www.gnu.org/philosophy/license-list.html">free software</a>, 
* but not public domain. By using it, you agree to the terms of the
* <a href="ObjCryst-license.txt">license</a>  (clarified artistic license),
* which you can find in the top directory of the ObjCryst++ package.
*\par Included libraries
* Software libraries distributed with this package (but not written by V. Favre-Nicolin)
* do not fall under the terms of the above copyright and license:
* these are the blitz++, atominfo, sglite and newmat package. You should
* refer to the documentation and legal notes included in their respective directories:
*\par 
* <a href="newmat.htm">The newmat matrix library</a> (Robert B Davies).
* This library is only used
* for matrix SVD decomposition and inversion in Least Squares methods (Least Squares
* method was written only for tests purposes... Don't expect anything from it...).
*\par 
* <a href="atominfo.html">The atominfo library</a>(Ralf Grosse-Kunstleve).
* This library is used to
* determine the atomic scattering factors for X-Ray, neutrons, as well as anomalous
* scattering factors and the atomic number.
*\par 
* <a href="sglite-license.txt">The SgLite library</a>(Ralf Grosse-Kunstleve).
* Note that this package
* is part of the <a href="http://pymol.sourceforge.net/">PyMOL Molecular 
* Graphics System</a> (used with permission),
* and is not free software. It is used to derive symetry operations from a given
* spacegroup symbol or number.
*
*\par 
* <a href="../../blitz/manual/index.html">The Blitz++ array library</a> (Todd Veldhuizen).
* This library
* is used for data storage and calculation in arrays and vectors. (Currently it
* is \e not used because of the huge memory requirements when compiling blitz++ expressions
* using the gcc compiler. Instead the CrystVector and CrystMatrix are used, emulating 
* the blitz interface, but without the smart handling of mathematical expressions.) Support
* for the Blitz++ should come back some time soon...
*/


/*! \page page_oop Why is Object Oriented Programming good for Crystallographic Computing ?
*\par Reusability 
* The \b encapsulation of data and \b inheritance 
* (and other features such as templates) allows a much better
* reusability of the code: adding new procedures does not require to understand how other
* function works, since these are only accessed through the object interface.
* And inheritance allows to create new object types (say new diffraction experiments) while
* re-using computing code from the base object: for example all diffraction analysis involves
* the computing of geometrical structure factors ( geometrical meaning before taking 
* scattering power into account), so that different diffraction ( X-Ray, neutron, electron,
* single crystal or powder) can re-use this code. See \ref page_design for examples
* on how to use inheritance to create new classes using inheritance.
* 
* \par Philosophy
* Re-usability in Object-oriented programming can be achieved by designing a base
* class interface which can be re-used for all derived classes using inheritance. It
* is much better than having isolated functions or small classes, since not only
* old code can be re-used, but the new code is compatible with the old because it has
* a common interface.. (of course the base interface must be well written and that
* is difficult to achieve... Hope I did not do a too bad job of it...).
*
*\par Performance 
* Although it has long been considered that the price of OOP was a slower
* execution speed, modern scientific & computing libraries have made that wrong, and c++
* is now widely used for large-scale computations. A short explanation is that c++
* is both a high-level language (object-oriented), and low-level since it uses
* pointers and thus programming can be very close to assembly. This latter point is a boon
* for numerically-intensive programs.
*\par 
* For good examples of numerically intensive calculations using c++, you take a look at:
*\par 
* -> The home page for Scientific Computing in Object-Oriented Languages:
* http://www.oonumerics.org/
*\par 
* -> \b POOMA (http://www.acl.lanl.gov/pooma/)Parallel Object-Oriented Methods and Applications
*\par 
* -> \b Blitz++ home page (http://www.oonumerics.org/blitz/):
* a 'smart' array computing library which yields performance on par
* with that of fortran.
*/


/*! \page page_design The Design of ObjCryst++
* If you have any comment or suggestion about this you can 
* <a href="mailto:vincefn@users.sourceforge.net">drop me an email</a>.
*
*  Also read \ref page_oop.
*
* \section overview Overview of the Library (Crystallographic classes)
* \par Scatterer
*   A Scatterer is the common denominator for any scattering object: all it includes
* is a function which gives a list of positions in fractional coordinates, with a
* scattering power associated to each position (Scatterer::GetScatteringComponentList()). 
* It also includes a few functions for the display (3D) of the object.
* \par 
*   All Scatterer can be derived from such a class: Atom, ZScatterer. The advantage of using
* inheritance is that all derived classes \b must re-use the functions declared in the base
* class, so that \e any function which knows what a generic 'Scatterer' object (but does not
* know what an Atom or ZScatterer is) can still use any derived class.
* \par 
*  \e Further \e development \e example: currently there is no 'rigid body' object: if any  
* developper wants to create such an object, he just needs to make sure he rewrites the function
* GetScatteringComponentList(). Thus without any modification, the Crystal and ScatteringData
* classes will automatically be able to use this new object... since this RigidBody object
* is derived from the Scatterer class (which Crystal and ScatteringData know).
*
* \par ScatteringPower
* This class can compute the scattering, resonant and thermic factor for any ScatteringData
* object (eg a list of reflections with some metric information). The three member function
* GetScatteringFactor(), GetTemperatureFactor(),GetResonantScattFactReal(), and 
* GetResonantScattFactImag() can be used to get the corresponding factors for a list
* of reflections (and wavelengths for resonant terms) in a ScatteringData object.
* \par 
* The base class is designed to handle anisotropic factors: for this the index of
* the symetric position in the Spacegroup must be given.
* \par 
*  \e Further \e development \e example: currently only the interface to handle anisotropy
* has been written, but no code or derived class. But no matter what kind of anisotropy
* is added, it will always work with the base class interface.
* \par 
* \e Note: why always use a ScatteringData object as input (to compute scattering factors,
* for example), rather than, say, a list of HKL or a list of sin(theta/lambda) ? The
* first reason is that from a ScatteringData you can extract all these
* hkl's and sin(theta)/lambda. The second reason is that with such an approach, no matter
* how complex the derived classes are, you can always use the same interface
* (for isotropic thermic factors as well as anharmonic !), so that any function
* written with only the knowledge of the base class can use any derived class.
*  
* \par Crystal
* a Crystal is a unit cell with an associated SpaceGroup with a list of Scatterer.
*
* \par ScatteringData
* The ScatteringData is a base class which is basically a list of reflections with the
* ability to compute structure factors. The DiffractionDataSingleCrystal and
* PowderPatternDiffraction classes are derived from it.
*
* \section optim Optimization design
* \par RefinableObj
* The RefinableObj is the base class for almost all objects in the library. The advantage
* of such a design (see the inheritance tree on the ObjCryst::RefinableObj page) is that
* when you design an algorithm, you do not need to know what kind of object is refined. All
* you need to know is (i) how many parameters there are (ii) how to move these parameters
* (iii) how to access one or several 'cost function' for the optimized object (to
* characterize 'how good' the current configuration is). Indeed, the global optimization
* class (for simulated annealing and parallel tempering) does not include any of the
* crystallographic headers, and yet it can refine the crystal strutures...
* \par 
* This design does not mean that only 'stupid' algorithms can be handled. Since the
* 'random moves' are handled by the refined objects, this 'random moves' can be very
* non-random (for example in the Crystal object, permutation of Scatterers is made from
* time to time...).
* \par 
* What if you are \e not interested in the RefinableObj functionnality ? You can simply
* ignore it, it will not do any harm. You can do other crystallographic work by 'forgetting'
* that a Crystal, PowderPattern, Atom is a RefinableObj. Even new derived objects do not
* have to declare their parameters as 'RefinablePar', if you want to save some memory.
* \par 
* what is currently lacking in RefinableObj is (i) a way to set constraints/restraints
* (currently there are only limits for individual parameters), the ability to have
* arrays of RefinablePar (to handle large structure without a significant memory penalty),
* and a design for analytical derivatives (well I've got a few ideas about this but
* this is not a priority yet...).
*/



