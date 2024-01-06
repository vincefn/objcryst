.. _tutorials:

*********
Tutorials
*********

Quickstart
==========
As an easy demonstration (even for non-crystallographers!), you can launch Fox, then use
the top menu File->Load and open the example file example/pbso4-joint.xml.
This is a file showing the global optimization of a very simple and well-known structure,
PbSO4 (lead sulfate). Then:

* Click on the **Crystal** tab, and in the window, choose the menu *Display->3D Display*. This opens
  a window with the crystal structure. You can click& drag with the mouse to rotate the structure
* Click on the **Powder Diffraction** tab, and for the two ``PowderPattern`` objects, choose the menu
  *Pattern->Show Graph*. This will show you the current observed (blue) and calculated (red) patterns.
* Click on the **Crystal** tab, and choose *Parameters->Randomize Configuration*
  (it would be too easy to start from the solution)
* Click on the **Global Optimization** tab. Then choose *Optimize->Run*. This will launch the global
  optimization running. You will see the structure evolve, as well as the Cost functions
  (the two Chi^2^ in that case).

.. toctree::
   :maxdepth: 1
   :caption: Tutorial structure determinations for organic and inorganic structures

   tutorial-cimetidine
   tutorial-pbso4
