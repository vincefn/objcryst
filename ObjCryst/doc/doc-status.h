/*! \page page_status               Current status and development history of ObjCryst++
*
* If you have a question you can 
* <a href="mailto:vincefn@users.sourceforge.net"> drop me an email</a>
*
*/
/*\section status Current Status
* The Library is moving along beta stage. The interface has been changed to
* a better design, especially for the scattering/diffraction class. The RefinableObj
* is also being improved to be more 'generic' and easier to use with varied algorithms.
*\par Spacegroup
*  Thanks to the sglite package, all SpaceGroup can be recognized from their symbol. 
* An AsymmetricUnit is also generated for each spacegroup (actually it is a parallelepipedic
* approximation). Only 3D spacegroup are allowed.
* We will eventually (that is, when I find time... at the autumn 2001, move
* to the <a href="http://cctbx.sourceforge.net/"> cctbx library </a>
*\par Crystal
* For a Crystal, definition of spacegroup and adding of various types of scatterers 
* are working. Dynamical occupancy correction 
* for global optimizations has been implemeted  (to correct overlap of identical
* atoms due to sharing of corner atoms in polyhedra, or due to symmetry). An anti-bump
* cost function (simple) is also implemented, and also enables atoms to merge continuously.
* An output to povray file is provided, as well as an OpenGL display of the cystal structure
* (only tested using Mesa under linux, with glut). Most important lacking features are
* saving structure to CIF files.
*\par Scatterer
*  Various type of Scatterer are provided : simple Atom but also ZScatterer, using a
* Z-Matrix Description, and derived form it are ZPolyedron (tetrahedron...icosahedron..).
* other types of scatterer can easily be added. So far the scattering power and thermic
* factor must be isotropic. Scattering factors for X-Rays (Thomson using the interpolated
* values and resonant (anomalous using either Sasaki or Henke tables) and neutrons are
* provided.
*\par ScatteringData
*  PowderPattern is the most developped class and its support is sufficient
* for most Global Optimization work. DiffractionDataSingleCrystal is less being worked on,
* but it can be used provided that the data has been corrected and merged beforehand.
*\par Global optimization 
*  GlobalOptimObj provides an algorithm for the so-called 'ab initio' structure determination
* from diffraction data using either a simulated annealing or a parallel tempering
* algorithm. Genetic algorithm may be added later.
*\par Least Squares
*  LSQObjNum class provides a rough least squares support, but it is definitely not
* the purpose (as of yet, at least) of this library. The most important lack is
* that numerical derivatives are systematicllay used. I wrote this only for testing
* purposes and have not sed it for some time, so try it at your own risk !
*\par Graphical user interface
*  Most objects have a graphical counterpart, which has been created using the
* <a href="http://www.wxWindows.org"> wxWindows</a> (linux,unix,windows,mac,...)
* library. The result is the 
* <a href="http://objcryst.sourceforge.net"> Fox program</a>, but the GUI
* interface can be used for any other purpose. It is of course possible
* to use ObjCryst++ without the wxObjCryst part.
*/
/**
*
*\section history History
*\par 1.1 
* \li INCOMPATIBLE CHANGES:GlobalOptimObj has been replaced by MonteCarloObj,
* RefinableObj::GetClassName now returns a "const string&", and 
* RefinableObj::BeginOptimization has two parameters, to enable
* approximations and restraints
* \li now one can choose either float or double for precision. The default
* is float for global optimization.
* \li Changed RefinableObj::BeginOptimization to pass two flags to enable
* approximations (faster) and restraints before beginning an optimization.
* \li ZScatterer: worked on ZScatterer::GlobalOptimRandomMove(),
* to give the possibility of making 'smart' moves, especially for
* molecules.
* \li Scatterer: removed the ZScatterer::Update() function from the base
* Scatterer class, as previously scheduled.
* \li OptimizationObj: Forked the algorithms classes, with a base OptimizationObj
* class, derived (currently) to a MonteCarloObj class (which replaces the
* old GlobalOptimiObj class, for simulated annealing and parallel
* tempering). Also added a tentative GeneticAlgorithm class, still
* in very early development and not usable yet for the common mortal.
* \li restraints: added base Restraint class, not really used so far.
* \li cleaned up a bit wxCryst classes in the hope to remove any GUI
* bugs. Still some crashes under windows, unfortunately...
* 
*\par 1.02(2001-nov 12th)
* \li PowderPattern: Added the input format for Sietronics (.cpi) files.
* Now importing data will allow null points without crashing. Also
* when no phase (background, crystal) has been added to the pattern,
* the calculated pattern returns a constant value (1.).
* \li ZScatterer: corrected a bug in the fractionnal coordinates
* calculation for mono/triclinic unit cells (thanks Mark Edgar).
*
*\par 0.9.1(2001-09-20-)
* \li ZScatterer: corrected a nasty bug(thanks Yuri Andreev),
* which made the Z-Matrix interpretation completely wrong.
* \li wxRefPar: now all refinable parameters have a local menu
* which can be used to remove or change limits.
* \li wxZScatterer: now can globally change relative limits to
* bond lengths and angles.
*
*\par 0.9(2001-09-18-first FOX Release)
* \li wxCryst is now working, allowing the first compilation of Fox.
* \li The data format to save files has been changed to an xml format
* (see http://www.w3.org/XML/).
* \li the definition of the refinable parameter types (RefParType)
* are now made using pointers, in a tree-like fashion. This allows
* an easier modification of status for groups of parameters.
* \li Lots of other modifications...
*
*\par 0.5(july 2001)
* \li Design largely improved for reusability (based on inheritance).
* \li wxCryst (the GUI-Graphical User Interface) part is
* in fast-moving development, mostly working, but not included here as code is still
* ugly and its design is not stable. As of 10/july/2001, I can launch the program,
* create a crystal, add a PowderPattern object, a GlobalOptimization object and do
* the global optimization while looking at the 'live' evolution of the Crystal Structure
* (3D display with OpenGL) and of the PowderPattern !
* \li Saving structures, data to files has been implemented. See RefinableObj::XMLOutput().
* Refinement flags are not saved (yet).
* \li Changed all float to REAL precision. Performance hit=+15% on structure factor,
* and +25% on full powder spectrums.
*
*\par 0.6(march 2001)
* First public release.
*
*\par 0.5 (february, 2001): 
*\li ZPolyhedron
*  Added a ZPolyhedron class, based on the Zscatterer objects, to describe (almost) regular
* polyhedra.
*\li Interatomic distance calculations
* Improved speed computation of the internel distance table, using non-symetrical atoms
* and the Asymmetric unit cell.
*\li AsymmetricUnit Cell
* Finally implemented the AsymmetricUnit Cell class, to improve distance calculation.
* Currently, upon Spacegroup construction, the (or one of the) smallest parallelepiped
* (with sides parallel to the crystallographic axes) including an Asymmetric unit is computed,
* using steps of 1/24. (So technically it is not an asymmetric unit cell, if no asymmetric
* unit cell exists for the group with all vertices parallel to the unit cell faces)
*
*\par 0.4 (january 2001):
*\li ZScatterer
*  Added a Zscatterer class, were cluster of atoms or molecules can be defined using
* a Z-Matrix description (see http://chemistry.umeche.maine.edu/Modeling/GGZmat.html
* or http://www.arl.hpc.mil/PET/cta/ccm/training/tech-notes/model/node4.html for a description
* of what such a matrix is. The advantage of this description is that it will now be possible
* to describe molecules using relevant parameters (separating bond lengths, bond angles
* and dihedral angles), and that polyhedra (for inorganic crystals) can also be described by it.
*\li Special positions-Dynamical Correction
*     Added a function to compute a 'smooth' correction of population
*     as several atoms overlap (this is a dynamical correction during model search) 
*\li Compilation & Makefiles
*     Cleaned up directories & makefile to minimize platform dependencies. All 
* platform-dependant commands should be in rules.mak, in the top directory.
*
*\par 0.3 (2000, november 24) : 
*\li RefinablePar
*  Now Scatterer, DiffractionData and Crystal inherit are children of RefinableParList,
*which makes the 'refining' methods easier to write & use.
*\li Compilation on Win32-Cygwin
*  Compilation works using the win32 port of gcc by cygnus : http://www.cygwin.com.
*  
*
*\par 0.2 (2000, november 07) : 
*\li SpaceGroup
*     Moved to SgLite package, thanks to R. Grosse-Kunstleve.
*\li Powder Diffraction
*     Most of features needed for 'basic' powder profile spectrum generation/analysis
*     has been added (DiffractionDataPowder).
*\li Least-Squares
*     Added an LSQObjNum class, to perform Least-Squares Refinement. 
*     Num stands for 'numerical', since only numerical derivatives are used
*     for that object. Uses eigenvalue filtering to avoid correlations. 
*\li DiffractionData
*     is now an abstract base class, with children 
*     DiffractionDataSingleCrystal, DiffractionDataPowder.
*
*\par 0.1 :
*\li General
*     The basic classes (SpaceGroup, Crystal, Scatterer, Scatterer:Atom,
*     Scatterer:RegularPolyhedra, DiffractionData) are there. 
*\li Speed
*     Computing of structure factors works with a fair speed. 
*
*/
