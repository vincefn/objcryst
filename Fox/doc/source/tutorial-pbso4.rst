.. _tutorial_pbso4:

PbSO4 tutorial: solving an inorganic structure, using both X-ray and neutron powder patterns
--------------------------------------------------------------------------------------------

Prerequisites
^^^^^^^^^^^^^
* The refined unit cell parameters, and a possible spacegroup
* Some information about the unit cell contents (stoechiometry, possible building blocks)
* one or several diffraction data sets. **They are in the ``Fox/examples/tutorial-pbso4`` directory**.

If you have all this, you can launch Fox !

First step: create your Crystal Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* From the top menu ``Objects``, **create a new Crystal Object**.
* Click on the **Crystals** tab to see your created crystal.
* You can already **display your crystal structure**: in the Crystal, use the ``Display menu-> 3D Display``. You can use the mouse to change the display:

  * drag with the left mouse button down to change the orientation
  * drag with the middle button to change the distance/aspect ratio. You can also use the '+' and '-' keys to change the distance, if you have a two-mouse button.
  * right-button click will display a popup menu to update the Crystal Display when you have changed some parameter (sometimes the program will try to do it by itself).
* You can **name your Crystal** ("PbSO4" or "Lead Sulfate",...) in the Field just after 'Crystal' (where it is written 'Change Me!'). This will be useful to choose your Crystal Structures thereafter.
* You can **change the spacegroup** by entering either a symbol or the spacegroup number in the SpaceGroup Field (in our case, "Pnma" will do). [Note: If (after validation), the symbol reverts to the old entry, it means that it has not been understood.]
* You can **change the lattice parameters** (8.482 5.398 6.959). Use the right-button menu to update the 3D Display of the structure.
* To add atoms to your structure, you must first **create the atom types**, what is called *ScatteringPower* in the program. Go to the 'Scatterers' menu of the Crystal, and choose 'Add Atomic Scattering Power'. This will add a ScatteringPowerAtom a bit lower. You should change the name (eg 'Lead' or 'Pb'-this is free format, you can call it 'George' if you want) of this atom type, as well as the symbol corresponding to the atom - the symbol must correspond to a "standard" symbol or ion (e.g. "Pb" or "Pb2+"). Do this until you have the scattering powers for Pb, S and O.
* Then **add a Pb Atom** to the structure, using the ``Scatterers->Add Atom`` menu of the Crystal object. You will be prompted to choose one atom type from those you have entered (choose the lead atom type). The atom will be put by default at (0,0,0), which you can see by updating the 3D view (right-click...)
* Then **add a SO,,4,, tetrahedron**: to do this use the ``Scatterers->Add Tetrahedron`` menu of the Crystal Object. You will be prompted for the central atom type (the sulphur), the peripheral one (oxygen), and then the distance (roughly 1.5 Angstroem). a 'Molecule' will appear in the list of scatterers of the Crystal Structure. You can then change the name of the tetrahedron, as well as the names of the atoms, at your convenience. The geometry of the tetrahedron is set using bond lengths and bond angles restraints (it gives to the polyhedron limited flexibility).
* *NOTE*: you have **nothing to do for special positions**, as long as the **Use Dynamical Occupancy Correction** field is set to 'Yes' (the default).*

You're done with the Crystal structure. You can (should) save using the top Fox menu ``File->Save``. This will save everythoing as an xml file, using a specific format to ObjCryst++/Fox.

Second step: create the PowderPattern object (X-Ray)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Go to the next tab ``Powder Diffraction``, and then use the top ``Object`` menu to **create a PowderPattern object**. You can give this one a name, e.g. "PbSO4 X-ray"
* You can **import the X-Ray data** using the PowderPattern object ``Data->Import Fullprof Pattern``, and select the ``xray.dat`` file in the ``Fox/example/tutorial-pbso4`` directory.
* You can **display the powder pattern** using the ``Pattern->Show Graph`` menu. You can click &amp drag with the left button to zoom in (double-click to unzoom). You can right-click on the graph to update the graph if you have changed manually a parameter. The coordinates (2theta, intensity of the mouse pointer is displayed at the bottom).
* Then **add the background phase**, using the ``Component->Add Bayesian Background (automatic)`` menu. This will automatically estimate the best background.
* Now **add the Crystalline phase**, using the ``Component->Add Crystalline Phase``, which will prompt you to choose one crystal structure available (the one you have already defined) from its name. (NOTE: for multi-phased powders, you can add several crystalline phases)
* Now **set the correct wavelength**: use the PowderPattern ``Radiation->X-Ray Cu Tube Ka12`` to select the correct radiation. This will update the entry field, and normally also the graph.
* *Nota bene*: the scale factor is normally horribly wrong before launching the optimization. You can, however, fit manually the scale factor by using the ``Pattern`` menu of the PowderPattern object... the fit will of course remain wrong as long as the atom positions are far from the correct ones.
* OK, to obtain reasonable fit you need to **choose adequate profile parameters**. Approximate parameters will be sufficient for our needs, so using a W parameter of .01 with U=V=0, a pseudo-Voigt profile with Eta0=0.5 and Eta1=0 will do. (you can change the values and update the graph to choose the values). Of course, if you did profile fitting, you can use the obtained values.

* NOTE: generally it is a good idea to run a bit the optimization before choosing the correct profile parameters (this way there are more reflections at the correct scale).
* You can also see that a **2theta zero** (shift) of -0.02 is necessary (again, this is approximate but sufficient)
* Finally, for a global optimization it is not necessary to use the entire pattern, so **put 0.25 in the max sin(theta)/lambda**. (this discards a lot of information but as you'll see it was useless information).

Maybe you should save now ?

Second step (bis): create the PowderPattern object (neutron)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Just repeat the same steps (create a second PowderPattern object) (you can display both graphs at the same time), and choose neutron radiation with a wavelength of 1.909. A W parameter of 0.25, and a pseudo-Voigt witha Eta0=0.15 should be fine. Again, put a value of 0.25 in the max sin(theta)/lambda field.

Third step: create the Global Optimization object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Go to the last tab of the Fox window ``Global Optimization``, and use the main window's ``Objects`` to **create a new Monte-Carlo Object**.
* Then we need to tell the algorithm what it is going to optimize (which objects), so we need to declare both Powder Patterns and the Crystal Structure. So use the `Optimized Objects->Add object to optimize` menu, and add all three objects (one at a time). (This is where you see it is useful to set meaningful names for all objects !).
* As of Fox 1.5, it is not necessary any more to tell the algorithm what criterion it should use to test the configurations: it will automatically take the *-log(Likelihood)* of all added objects, in this case the sum of all Chi^2^.
* OK, you're all set ! Normaly, you should not change the choice of algorithm (**Parallel Tempering** is much more efficient than Simulated Annealing), nor the temperature or displacement amplitude parameters. They should work with any structure and combination of data.

Now would be a good time to save using the top ``File->Save`` menu.

* So **Launch the optimization**, using the ``Run->Multiple Runs`` menu of the Monte-Carlo object. If you have left the 3D Crystal structure window and the powder pattern graphs opened, they should be live-updated. Convergence is quick in this case, normally less than 50 000 trials. Normally, **you should first change the *number of trials per run** (to about 50 000 or 100 000 for this simple structure), and the algorithm would restart the optimization at the end of the run, providing one solution for each 'run'.

* You can also follow the progression of the Chi^2^ statistics, GoodnessOfFit (=Chi^2^/nb,,obs,,), R,,wp,, and R,,p,, in the powder pattern objects. But there's nothing better than the eye to tell whether the fit is good or not.
* **Stop the optimization** using the same menu when you have reached a satisfactory solution.
* Once the algorithm has converged, if you change the max sin(theta)/lambda fields in both powder patterns, and update the graphs (right-click...), you should see that the part of the pattern which was not used as a criterion is relatively well fitting (not great at high angle for the neutron, but this will improve when you do the Rietveld  refinement !).
* If you have done a "Multiple Run", you can **browse the solutions** using the Solutions menu. Click on any solution and it will automatically update the display of the structure and powder pattern(s). You can then choose whichever looks better.

Last step: export the crystal Structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* You can **export the atomic coordinates** by going to the Crystal structure, and use the ``File->Save as text``
  menu. this will save a file with all atom fractional coordinates and occupancies. You should see that all atoms
  have a 'dynamical occupancy' of 0.5, which is due to the fact that they are all (except two oxygens) on a special
  position. The other two oxygens are in fact equivalent, as the distance table included in the output shows.
* You can also **export to a CIF file** using the same menu. Normally, atoms which are fully overlapping are taken care of so that only a single atom remains.
