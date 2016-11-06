=====================================
F.O.X.: Free Objects for Xtallography
=====================================

Introduction
============
Fox, 'Free Objects for Crystallography' is a free, open-source program for the ab initio structure determination from powder diffraction.

The FOX program (for Linux, MacOS X and windows) was made for the ab initio crystal structure solution from diffraction data (mostly powder diffraction data). 

Features
========
The most interesting features for ab initio structure determination are:

- a versatile description of the crystal contents: either isolated atoms , molecules described using a bond length, bond angles and dihedral angles, and polyhedra for inorganic compounds. You can describe your structure by using any combination of groups of atoms, using a chemist's or crystallographer knowledge about the connectivity in your sample to constrain possible solutions..
- an automatic correction for special positions and shared atoms between polyhedra, suitable for global optimization algorithms.
- the ability to find peaks in a powder pattern, index the unit cell, and explore possible spacegroups
- the ability to use simultaneously multiple powder patterns (X-rays, neutrons), as well as single crystal data (e.g. extracted from a powder pattern)
- smart global optimization algorithms which can get out of false minima.
- a graphical interface (see the screenshots) with a 3D crystal structure view, with live updates during the optimization process.
- the ability to view CIF files

So, if you:

- have an unknown compound but with (approximately) known composition.
- have a powder pattern (X-Ray or neutron or both)
- would like to solve the structure (i.e. find the atom -or group of atoms- positions, before refining them with another package).
 
Then Fox can help you.

Viewing crystal structure, diffraction data and CIF files
=========================================================
FOX can also be used to display CIF (Crystallographic Information Files).

It can thus be used also for educational purposes, to show a 3D display of Crystal structures, and the associated powder pattern(s) (see how adding atoms, changing the lattice, or changing the spacegroup affects the powder pattern and the 3D structure).

In the latest version (2016.1), it is possible to query directly the Crystallography Open Database (http://www.crystallography.net) from FOX, and open the associated crystal structures.

Development
===========
If you would also like to choose your own criterion and algorithm to solve the structure, then it will be even better : FOX is built on a very customizable and expandable library (ObjCryst++), which allows you to evaluate your Crystal structure following a combination of criteria from diffraction data and interatomic distances. Other criteria for optimization and other algorithms can easily be added.

The main API documentation can be found at: https://vincefn.net/ObjCryst/classes.html

Documentation
=============
The main documentation for FOX can be found at: https://fox.vincefn.net

Installation
============
Follow the instructions at: http://fox.vincefn.net/Install

Tutorials
=========
New to FOX ? See: http://fox.vincefn.net/Tutorials
