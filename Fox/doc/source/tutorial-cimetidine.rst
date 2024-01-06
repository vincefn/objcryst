.. _tutorial_cimetidine:

Cimetidine tutorial: solve an organic structure using an x-ray powder pattern
-----------------------------------------------------------------------------

Prerequisites
^^^^^^^^^^^^^
* The refined unit cell parameters, and a possible spacegroup
* The molecule formula (ideally a z-matrix or pdb file)
* one or several diffraction data sets. They are in the ``Fox/examples/tutorial-cimetidine`` directory.

If you have all this, you can launch Fox !

First step: create your Crystal Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* From the top menu ``Objects``, **create a new Crystal Object**.
* Click on the ``Crystals`` tab to see your created crystal.
* You can already display you crystal structure: in the Crystal, use the ``Display menu-> 3D Display``.
  You can use the mouse to change the display:

  * drag with the left mouse button down to change the orientation o drag with the middle button to change
    the distance/aspect ratio. You can also use the '+' and '-' keys to change the distance, if you have a two-mouse button.
  * right-button click will display a popup menu to update the Crystal Display when you have changed
    some parameter (sometimes the program will try to do it by itself).
* You should **give your Crystal a name** ("Cimetidine" or "My beautiful Molecule", "Napoleon 1er"...)
  in the Field just after `Crystal` (where it is written `Change Me!`). This will be useful to choose
  your Crystal Structures thereafter (note that if you have several Crystal structures, they
  **must all have different names**).
* You can **change the spacegroup** by entering either a symbol or the spacegroup number in the
  `SpaceGroup` Field (in our case, `P121/a1` will do). [Note: If (after validation), the symbol reverts
  to the old entry, it means that it has not been understood.]
* You can set the **lattice parameters** (10.394, 18.819, 6.825, beta=106.44). Use the right-button menu to
  update the 3D Display of the structure.
* Generally before adding atoms you must first **create the atom types**, what is called `ScatteringPower`
  in the program (see the [:FoxTutPbSO4:PbSO4 tutorial] for that). As we will be importing the molecule
  structure from a z-matrix file, this will be done automatically.
* **Creating the Molecule**: The molecule is described using a list of atoms bound together by a list of
  bond lengths, bond angles and (more rarely) dihedral angle restraints. You can create one empty molecule,
  and add atoms and restraints manually, but **the simplest way to input a molecule is to import its structure
  from a Fenske-Hall Z-matrix**. A Z-matrix describes all atom positions from a first atom and bond distance,
  bond angles and dihedral angles with the other atoms (see http://chemistry.umeche.maine.edu/Modeling/GGZmat.html
  for a description). Use the Crystal ``Scatterer->Import Molecule from Fenske-Hall Z-matrix`` menu. FOX
  can (so far) only import Fenske-Hall Z-matrices, but it is possible to transform a wide range of molecule
  structures to this type of files using either **babel** (http://www.eyesopen.com/babel/) or **openbabel**
  (http://openbabel.sourceforge.net). To do this just use: ``babel -ipdb cime.pdb -ofhz cime.fhz -d``
  (the -d option will get rid of the hydrogens). You can see the cime.pdb and cime.fhz files in the
  example/tutorial-cimetidine directory.
* Check **Restraints** within the molecule: look at the list of restraints or, look at the displayed 3D
  structure if there are no extra or missing bonds (Fox tries to add bonds depending on interatomic distances
  and tabulated atomic radius, so it can be wrong). Normally it should be OK.
* Fox will automatically decide which bonds are free torsion angles for the Global Optimization.

You're done with the Crystal structure. You can (should) save using the top Fox menu ``File->Save``.
This will save everything as an xml file, using a specific format to ObjCryst++/Fox.

Second step: create the PowderPattern object (X-Ray)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Go to the next tab **Powder Diffraction**, and then use the top 'Object' menu to
  **create a PowderPattern object**. You can give this one a name, e.g. "Cimetidine X-ray"
* You can **import the X-Ray data** using the PowderPattern object *Data->Import Fullprof Pattern*,
  and select the cime.dat file in the Fox/example/tutorial-cimetidine directory.
* You can **display the powder pattern** using the `Pattern->Show Graph` menu. You can click & drag
  with the left button to zoom in (double-click to unzoom). You can right-click on the graph to update
  the graph if you have changed manually a parameter. The coordinates (2theta, intensity of the mouse
  pointer is displayed at the bottom).
* Then **add the background phase**, using the ``Component->Add Bayesian Background (automatic)`` menu.
  This will automatically estimate the best background.
* Now **add the Crystalline phase**, using the ``Component->Add Crystalline Phase``, which will prompt
  you to choose one crystal structure available (the one you have already defined) from its name.
  (*NOTE: for multi-phased powders, you can add several crystalline phases*)
* Now would be a good time to **set the correct wavelength**: just input 1.52904 in the wavelength field.
  This is synchrotron data, so you should put ~0.98 in the linear polarization rate field.
* OK, to obtain reasonable fit you need to **choose adequate profile parameters**. Approximate parameters will
  be sufficient for our needs, so using a W parameter of .001 with U=V=0, a pseudo-Voigt profile with
  Eta0=0.5 and Eta1 will do. (you can change the values and update the graph to choose the values).
  Of course, if you did profile fitting, you can use these values.
* Now if you want to **scale the calculated powder pattern to the observed one**, use the powder pattern
  menu ``Pattern -> Fit Scale Factor For R(w)``.
* Finally, **for a global optimization it is not necessary to use the entire pattern, so put 0.25 in the
  max sin(theta)/lambda**. (this discards a lot of information but as you'll see it was useless information).

Maybe you should save now ?

Third step: create the Global Optimization object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Go to the last tab of the Fox window 'Global Optimization', and use the main window's 'Objects'
  to **create a new Monte-Carlo Object**.
* Then we need to tell the algorithm what it is going to optimize (which objects), so we need to
  **declare the powder pattern and the Crystal Structure**. So use the
  *Optimized Objects->Add object to optimize* menu, and add both objects (one at a time).
  (This is where you see it is useful to set meaningful names for all objects).
  It will automatically use the powder pattern's Chi^2 as a criterion for convergence.
* OK, you're all set ! Normmaly, you should never change the choice of algorithm
  (**Parallel Tempering** is better than Simulated Annealing), nor the temperature or displacement
  amplitude parameters. They are supposed to work with any structure and combination of data.

  Now would be a good time to save using the top ``File->Save`` menu. You can activate the auto-save option
  (every hour is a fine setting if you think it's going to take a long time).

* Choose the **number of trials to use per run**: in this case 2 million trials will be ok.
* Now you can **Launch the Optimization**, using the ``Run->Multiple Runs`` menu of the Monte-Carlo object.
  If you have left the 3D Crystal structure window and the powder pattern graphs opened, they should be
  live-updated. Convergence is moderately slow in this case, on average 1.5 million trials, which can take
  5 minutes to 30 mn depending on your computer speed.
* You can follow the progression of the Chi^2^ statistics, Goodness Of Fit (Chi^2^/nb,,obs,,), Rwp and
  Rp,, in the powder pattern object. But there's nothing better than the eye to tell whether the fit is good or not.
* When satisfactory, use the menu to stop the optimization. You can compare to the optimized
  results in the fox examples.
* If you have done a "Multiple Run", you can **browse the solutions** using the Solutions menu.
  Click on any solution and it will automatically update the display of the structure and powder
  pattern(s). You can then choose whichever looks better.

Last step: export the solved crystal Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* You can **export the atomic coordinates** by going to the Crystal structure, and use the
  *File->Save as text* menu. this will save a file with all atom fractional coordinates and
  occupancies. You can compare the structure obtained with the one already refined in the fox example directory.
* You can also **export to a CIF file** from the same menu.
