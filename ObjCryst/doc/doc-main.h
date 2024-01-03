/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-     Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*  doc-main.h
*  This file includes only documentation, to be formatted automatically by doxygen.
*/

// Using thr namespace allows links to classes to be generated without the ObjCryst:: prefix
namespace ObjCryst
{

/*! \mainpage %ObjCryst++ Object-Oriented Crystallographic Library
**
*\section index Index
*  \li <a href="https://github.com/vincefn/objcryst">Project page on GitHub</a>
*  \li <a href="https://sourceforge.net/projects/objcryst/">Historical project on sourceforge (obsolete)</a>
*  \li \ref page_legal
*
*\section intro Introduction
*\par Overview of the project
*
* The aim of this project is to provide object-oriented tools for Crystallography. The
* focus of these tools is the analysis of diffraction data (mostly powder diffraction),
* with an emphasis on speed optimization for repetitive calculations (for global
* refinement methods such as simulated annealing, for example).
*
*\par
* To have an <b>overview of all available libraries</b>, go to the <a href="hierarchy.html">Class Hierarchy</a>.
* You are probably interested in the main Crystallographic classes, including:
* * Crystal: the base class describing a complete crystal structure
* * Scatterer: the base class for a scatterer, including:
*   * Atom
*   * Molecule
* * ScatteringData, including:
*   * DiffractionDataSingleCrystal
*   * PowderPatternDiffraction
* * PowderPattern: describing a complete powder pattern with multiple phases, background
* * PowderPatternComponent: can be derived to:
*   * PowderPatternBackground
*   * PowderPatternDiffraction
* * SpaceGroup: with all the spacegroup and symmetry information, interfacing with cctbx
* * RefinableObj: the base refinable object with a list of parameters - all
*   objects are derived from it.
* * RefinablePar: the base refinable paremeter, which is used in all RefinableObj
* * Algorithm classes:
*   * OptimizationObj, from which MonteCarloObj is derived
*   * LSQNumObj: least-squares minimisation
*   * CellExplorer: class for the volume dichotomy algorithm
*
*\par
* Even if we intend to use this library mainly for the development of a global optimization
* program from powder diffraction, this library is programmed in a general way so that
* other applications can make use of it.
* The library was designed to be reusable, by adding new kind of experiments, new algorithms,
* new Scatterer type, new ScatteringPower,...
*
*\par
* <b>There are two main usage of the ObjCryst++ library today, and you are encouraged
* to start using those, unless you already know you want to use the C++ library</b>:
* * \b Fox, 'Free Objects for Xtallography' is an application for Linux, macOS and windows
*   which can be <a href="https://github.com/vincefn/objcryst/releases">downloaded from GitHub</a>,
*   and can be used mainly to solve crystal structures from powder diffraction
* * \b pyobjcryst is the python module which can be used as an interface to a large
*   part of the ObjCryst++ library: see the <a href="https://pyobjcryst.readthedocs.io/">documentation</a>
*   and the <a href="https://github.com/diffpy/pyobjcryst">github project</a>), including
*   <a href="https://pyobjcryst.readthedocs.io/en/latest/examples/index.html"> example notebooks</a>.
*
*\par Contributors
* Design & implementation: Vincent Favre-Nicolin (http://v.favrenicolin.free.fr/),
* <a href="mailto:vincefn@users.sourceforge.net">vincefn@users.sourceforge.net</a>
*\par
* This project was initiated in the laboratory of Crystallography of the University
* of Geneva (http://www.unige.ch/), and is part of the development of a global
* optimization program with Radovan Cerny (https://www.unige.ch/societes/SSCr/w3cerny/rcerny.htm).
* \par
* This project has been supported by the <a href="http://www.snf.ch/">Swiss
* National Science Foundation (#21-53847.98)</a>.
*
*\par
* This project also makes use of some other programs or libraries. Most notable in this
* project is the cctbx library originally written by R. Grosse-Kunstleve
* (see \ref page_legal).
*
*/

/*! \page page_legal License Information
*\par Copyright
*  The ObjCryst++ library is copyright (2000-2002)
* <a href="mailto:vincefn@users.sourceforge.net">Vincent FAVRE-NICOLIN</a>
* and (2000-2001) The University of Geneva (Switzerland)
*\par License
* This is <a href="http://www.gnu.org/philosophy/license-list.html">free software</a>,
* but not public domain. By using it, you agree to the terms of the
* <a href="ObjCryst-license.txt">GPL license</a>  (General Public License),
* which you can find in the top directory of the ObjCryst++ package.
* Should you wish to use this library in a closed-source project,
* please contact the author for possible conditions.
*\par Included libraries
* Software libraries distributed with this package (but not written by V. Favre-Nicolin)
* do not fall under the terms of the above copyright and license:
* these are the cctbx and newmat package. You should
* refer to the documentation and legal notes included in their respective directories:
*\par
* <a href="newmat.htm">The newmat matrix library</a> (Robert B Davies).
* This library is only used
* for matrix SVD decomposition and inversion in Least Squares methods (Least Squares
* method was written only for tests purposes... Don't expect anything from it...).
*\par
* <a href="https://cctbx.github.io/">The cctbx library</a>(Ralf Grosse-Kunstleve).
* This library is used to determine the atomic scattering factors for X-Ray,
* neutrons, as well as anomalous scattering factors and the atomic number,
* and for all spacegroup interpretation.
*\par
*
*/

} //namespace
