.. Fox documentation master file, created by
   sphinx-quickstart on Wed Jan  3 18:36:08 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

###############################
Welcome to Fox's documentation!
###############################
Fox, 'Free Objects for Crystallography' is a free, open-source program for the *ab initio*
**structure determination from powder diffraction.**


Beginning with Fox
==================

.. _introduction:

Introduction
------------

The FOX program (for Linux, MacOS X and windows) was made for the  *ab initio crystal structure solution from diffraction data* (mostly powder diffraction data). Its most interesting features for ab initio structure determination are:

* a **versatile description of the crystal contents**: either isolated `atoms` , `molecules` described using a bond length, bond angles and dihedral angles, and polyhedra for inorganic compounds. *You can describe your structure by using any combination of groups of atoms, using a chemist's or crystallographer knowledge about the connectivity in your sample to constrain possible solutions..*
* an **automatic correction for special positions and shared atoms between polyhedra**, suitable for global optimization algorithms.
* the ability to use simultaneously **multiple powder patterns** (X-rays, neutrons), as well as single crystal data (e.g. extracted from a powder pattern)
* **smart global optimisation algorithms** which can get out of false minima.
* a graphical interface (see the :ref:`screenshots <screenshots>`) with a 3D crystal structure view, with live updates during the optimization process.

So, if you:

* *have an unknown compound* but with (approximately) known composition.
* *have a powder pattern* (X-Ray or neutron or both), which you have already indexed (with the refined unit cell), and for which you have some idea of the possible spacegroup(s) (using systematic extinctions).
* *'would like to solve the structure* (i.e. find the atom -or group of atoms- positions, before refining them with another package).

Then Fox can help you.

This program can be used also for **educational purposes**, to show a 3D display of Crystal structures, and the associated powder pattern(s) (see how adding atoms, changing the lattice, or changing the spacegroup affects the powder spectrum and the 3D structure).

If you would also like to *choose your own criterion and algorithm to solve the structure*, then it will be even better: FOX is built on a **very customizable and expandable library** (`ObjCryst++ <https://objcryst.readthedocs.io/>`_), which allows you to evaluate your Crystal structure following a combination of criteria: currently available are the diffraction data Chi^2, an anti-bump and a bond-valence cost function, but it would be very easy to write a new criterion using (say) interatomic distances (energy of the configuration, analysis of the coordination). This, with a versatile way to describe the unit cell's contents, is what makes the algorithms used interesting. For example you can combine several datasets (X-Ray, Neutron...) for the same structure,...

.. toctree::
   :maxdepth: 1
   :caption: Get Fox

   get-fox

.. toctree::
   :maxdepth: 1
   :caption: Using Fox to solve Crystal structures

   tutorials
   manual-intro
   screenshots

.. toctree::
   :maxdepth: 1
   :caption: Bibliography

   biblio


Mailing List
============
To be kept informed about new releases, you are encouraged to **subscribe** to the (*very low*-volume)
FOX mailing list at: http://lists.sourceforge.net/mailman/listinfo/objcryst-foxx

The **archives** of the list can be found at https://sourceforge.net/p/objcryst/mailman/objcryst-foxx/
