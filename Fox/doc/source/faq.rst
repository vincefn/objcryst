.. _faq:

Frequently Asked Questions
==========================
General Questions
-------------------
* :ref:`My question is not answered here <gen1>`
* :ref:`How can I get news about the evolution (new releases etc..) of Fox ? <gen2>`
* :ref:`What is an "Object" ? <gen3>`

Global Optimization
-------------------

* :ref:`I know 'Simulated Annealing' (SA), why do you recommend to use "Parallel Tempering" (PT) ? <glob1>`
* :ref:`My optimization has reached (some cost). Has it converged enough ? <glob2>`
* :ref:`Fox does not optimize profile/unit cell parameters, even if I check the parameters to be optimized ! Why ? <glob3>`
* :ref:`What parameters can/should I optimize ? <glob4>`
* :ref:`What objects should I add in the Optimization object ? <glob5>`
* :ref:`How can I remove an object from an Optimization ? <glob6>`

Crystal Modeling
----------------

* :ref:`How can I take into account Special Positions ? What is the Dynamical Occupancy Correction ? <cryst1>`
* :ref:`What happens when two atoms are completely overlapping ? Does one disappear ? <cryst2>`
* :ref:`How can I ensure that the algorithm keeps the correct formula for my compound ? <cryst3>`
* :ref:`How can I use an Antibump criterion ? How can I see the Antibump parameters ? Does it prevent the "merging" of atoms using the Dynamical Occupancy Correction ? <cryst4>`
* :ref:`Should I restrict the positionnal parameters (fractionnal coordinates) to (0;1) using the parameter's limits ? Same for angles in (0,2pi)? <cryst5>`

Diffraction data
----------------

* :ref:`I copied the powder pattern 2theta zero (shift) from my profile fitting program, but it looks wrong ? <diff1>`
* :ref:`How can I change the powder pattern background points ? <diff2>`
* :ref:`Is it possible to exclude parts of a powder pattern ? Can I remove an exclude region ? <diff3>`
* :ref:`Can I use several data sets to refine a single crystal structure ? <diff4>`
* :ref:`You recommend to use only low-angle data, but why discard this information ? High-angle reflections *are* important for structure determination <diff5>`
* :ref:`What shall I put in the "linear polarization rate" field ? <diff6>`
* :ref:`Do you take into account the magnetic scattering of neutrons ? <diff7>`
* :ref:`The Rwp and graphs look fine, but the optimized preferred orientation parameters are wrong ! <diff8>`

General questions
-----------------

.. _gen1:

**My question is not answered here !**
  Don't hesitate to raise an issue on https://github.com/vincefn/objcryst/issues, or mail me.

.. _gen2:

**How can I get news about the evolution (new releases etc..) of Fox ?**
 To be kept informed about new releases, you are encouraged to subscribe to the (very) low-volume
 FOX mailing list at http://lists.sourceforge.net/mailman/listinfo/objcryst-foxx (the archives can be found at http://sourceforge.net/mailarchive/forum.php?forum_id=654).

 If you intend to use the !ObjCryst++ / Fox source code for your own projects, it is highly recommended to subscribe to the !ObjCryst++ mailing list at http://lists.sourceforge.net/mailman/listinfo/objcryst-devel, (the archives can be found at http://sourceforge.net/mailarchive/forum.php?forum_id=1689).

.. _gen3:

**What is an "Object" ?**
 Almost everything (Crystal, !PowderPattern, Scatterer, !ScatteringPower, !OptimizationObj) is considered an object. In Fox/!ObjCryst++ a generic approach is used, so that the algorithms can optimize them without actually knowing what they actually are: the idea is that the algorithm knows that these objects are "Refinable Objects", so that he can get from them (i) a list of parameters and (ii) a choice of Cost Functions. The advantage of such an approach is that if the modelization (parametrization) of some objects changes, the algorithms do not need any modification since they use a generic approach to the objects.

Global Optimization
-------------------

.. _glob1:

**I know 'Simulated Annealing' (SA), why do you recommend to use "Parallel Tempering" (PT) ?**
 In SA, the 'temperature' of the algorithm (think: high 'temperatures' allows for more improbable configurations, allowing to get of local minima) is decreased slowly following a predetermined law. If the decrease is too fast, you can easily be trapped in a local minimum.With PT, all temperatures are optimized in parallel, so that you can get of a local minimum at all times. The algorithm is invariant with time and therefore does not require from the user to choose the number of trial and the temperature schedule.

 In other words, always use Parallel Tempering, unless you are interested in algorithms and want to play with SA. **More general Idea: always keep the default algorithm and optimization schedule & parameters (temperature and displacement), they are tuned to work with any structure or data.**

.. _glob2:

**My optimization has reached [enter cost here]. Has it converged enough ?**
 Sorry, there is no general rule for just numbers... Nothing (for a rough idea at least) replaces the eye comparison of calculated and observed patterns (use the zoom by left-click & drag !). A R,,wp,, value of 0.15 can be sufficient (for the structure solution to be good) if the background is low, but it can also be crap if the background is high (with a background at 50 000 and the maximum intensity at 60 000, even Rwp=.05 can be very bad) !

.. _glob3:

**Fox does not optimize profile/unit cell parameters, even if I check the parameters to be optimized ! Why ?**
 "It's not a bug, it' a feature": Fox automatically fixes all profile and 2theta correction parameters, for several reasons:

 * It is a bad idea to use a Global Optimization algorithm derived from Monte-Carlo to search for this kind of parameters.
 * Profile parameters, 2theta shift (zero) can be determined before using Fox for the structure solution, using profile fitting (and therefore you should do it that way).
 * By using the "integrated R,,wp,," (iR,,wp,,) cost function (rather than R or R,,wp,,), the algorithm is almost insensitive to profile parameters. Personnally, for all structures I work on and for all examples given with Fox, I never use precise/refined profile parameters, but always adjust them "by eye" within Fox, which is sufficient thanks to the "i"Rwp approach. And generally I only use the 'W' term (constant width). You only need to refine these parameters to get a good full-profile agreement and please referees, but for the structure solution (finding approximate positions of atoms) part, only the integrated intensities matter.

.. _glob4:

**What parameters can/should I optimize ?**
 You can optimize positionnal parameters (translations, rotations and conformation), Biso's, occupancies, and texture parameters.

 However, you should begin to optimize positionnal parameters only (translations, rotations and conformation), after setting "expected" values for Biso , occupancies of 1, and without texture parameters.

 You may activate the optimization of occupancies after a first optimization if you believe it is necessary. Optimizing Biso is probably a bad idea, since at low angle (the only part used for global optimization) data is not sensitive to Biso.

 And do *not* activate texture parameters unless (i) you have a strong feeling it is necessary and (ii) you have not been able to obtain an un-textured data set. Having strong texture will significantly reduce the probability of finding the structure solution.

.. _glob5:

**What objects should I add in the Optimization object ?**
 If you are optimizing a Crystal structure with a single Powder Pattern, you should add these two objects. If you have several data sets, you must add all diffraction data objects.

.. _glob6:

**How can I remove an object from an Optimization ?**
 The only way to remove an object or cost function is to edit an xml file, and look at the very end: just remove the line corresponding to the object/cost function you want to remove, and voila ! NOTE: do not remove an object if you are still using on of its cost functions !!

Crystal Modeling
----------------

.. _cryst1:

**How can I take into account Special Positions ? What is the Dynamical Occupancy Correction ?**
 To take into account special positions, activate the "use dynamical occupancy correction" option in your Crystal structure interface. Do not change the occupancy of any atom. When atoms (of identical type) overlap, their occupancy will be automatically reduced (in the background, you won't see it) so that n atoms overlapping will all have an occupancy of 1/n. This allows to correct any type of overlapping, for an atom on a special position, or atoms being merged.

 It is strongly recommended to activate the Dynamical Occupancy Correction for all inorganic structures. Generally it will be useless for molecular structures, as special positions or shared atoms between molecules are extremenly rare.

.. _cryst2:

**What happens when two atoms are completely overlapping ? Does one disappear ?**
 When (say) two atoms fully overlap (assuming you have activated the Dynamical Occupancy Correction), both their occupancies are set to 0.5, so that the sum of the two atoms look just like one. However both the atoms are still present in the global optimization, and they will be able to move away from each other later if it gives a better solution.

.. _cryst3:

**How can I ensure that the algorithm keeps the correct formula for my compound ?**
 You have to put at least as many independent atoms as there are in the final result ; that is, you have to guess approximately how many atoms will be required. And thanks to the Dynamical Occupancy Correction, you can put too many atoms, and excess atoms will automatically be merged. Of course, do not exaggerate the number of excess atoms or the algorithm may converge very slowly.

.. _cryst4:

**How can I use an Antibump criterion ? How can I see the Antibump parameters ? Does it prevent the "merging" of atoms using the Dynamical Occupancy Correction ?**
 To use the antibump:

 * To add (or remove) an antibump distance, use the [[FoxRefGUICrystal#scattpow| Scattering Power window]]
 * the cost will automatically be added to the overall cost (as long as you added the Crystal structure to the list of objects to be optimized). If the cost is too small with regard to the Chi^2^ of your diffraction data, you can use the "scale" factor displayed beside the antibump cost to increase it. Be warned that a too large antibump factor will prevent the optimization to work correctly; therefore, do not mistake the antibump function for an energetic description of your crystal structure.

 For identical types of atoms, the antibump penalty vanishes to zero when the atoms get near each over (practically, the penalty goes from reaches its maximum when the distance decreases from the given d,,min,, to d,,min,,/2, and then decreases to 0 when the atoms fully merge). So the antibump cost function does not prevent merging.

.. _cryst5:

**Should I restrict the positionnal parameters (fractionnal coordinates) to (0;1) using the parameter's limits ? Same for angles in (0,2pi)?**
 No, these are "periodic" parameters, and are automatically corrected, e.g. changed to 0.01 if they get to 1.01 (using limits is then a bad idea, as it restricts movements).

Diffraction data
----------------

.. _diff1:

**I copied the powder pattern 2theta zero (shift) from my profile fitting program, but it looks wrong ?**
 Try to change the sign and/or to multiply it by a factor 2, and see if it is better (after each change, right-click in the powder graph window to update the display).

.. _diff2:

**How can I change the powder pattern background points ?**
 Simply change the 2-column (2theta - intensity) file you used to import the points, and re-import it. It does not matter if the number of points used has changed.

.. _diff3:

**Is it possible to exclude parts of a powder pattern ? Can I remove an exclude region ?**
 Yes, in the !PowderPattern Object, use the Pattern->"Add 2theta exclude region", and enter the 2theta limits. You can add several ones.

 NOTE: the powder patterns in these regions are still computed, so if you want to avoid computing the higher parts of the powder pattern, rather use the "Max sin(theta)/lambda" entry field.

 To remove it, you must save in an xml file, and then remove the line declaring the region, which should be something like : <!Exclude2Theta>0 5</!Exclude2Theta>

.. _diff4:

**Can I use several data sets to refine a single crystal structure ?**
 Yes, generally speaking, you can add any number of data sets. You just have to add them all to the !GlobalOptimization object, and declare a cost function (the iRwp preferably) for each data set.

.. _diff5:

**You recommend to use only low-angle data, but why discard this information ? High-angle reflections *are* important for structure determination !**
 First of all, a clear distinction should be made between (i) structure solution and (ii) structure refinement, which are both part of a complete structure determination:

 1. "Structure solution" consists in starting from relatively little information about the structure (one data set, and chemical composition), and one way or another to determine the approximate (0.1-0.5 Angstroem) position of all atoms.
 2. "Structure refinement" consist in starting from known approximate atoms positions, and to determine these positions, as well as temperature factors, occupancies, with the best precision possible.

 Fox only targets (1), and therefore only requires enough data to find these "approximate" positions, which obviously is low-angle data (for a low-resolution structure determination). You should then use sin(theta)/lambda limits of between 0.25 (2 Angstroem) and 0.5 (1. Angstroem). I suggest you try the examples to convince yourself. After the optimization, you can change the sin(theta)/lambda limit to check that it is indeed the correct solution.

 Using too high-angle data will slow the structure solution by two factors, (a) because the number of reflections is proportionnal to the cubic power of sin(theta)/lambda (and therefore much slower), and (b) because high-angle data will only contribute to confuse the algorithm with useless information (to make a simple analogy, our goal consists to drive by car from Paris, to find a street in Geneva ; although all we want is to find the correct street, having the "high resolution data" is equivalent to know the color of the phone booth in that street, which is perfectly useless).

 Of course, the high angle data will be useful... but for the refinement stage (ii), i.e. in another program.

.. _diff6:

**What shall I put in the "linear polarization rate" field ?**
 The linear polarization rate is around 1.0 for synchrotron data (typically 0.95-0.98), and 0 for laboratory diffractometer (unless you use a monochromator). A precise value is not critical for global optimization, as we mostly use low-angle data, and as errors can be swallowed by the overall B-factor.

 The calculation of the polarization factor assumes that the polarization plane is horizontal, and the diffracting plane is vertical (i.e. standard equatorial geometry to maximize intensity on a synchrotron).

 This is ignored for neutron diffraction.

.. _diff7:

**Do you take into account the magnetic scattering of neutrons ?**
 No.

.. _diff8:

**The Rwp and graphs look fine, but the optimized preferred orientation parameters are wrong !**
 With preferred orientation, it is easy to find several close configuration, by inverting the March coefficient and using a perpendicular vector of preferred orientation. This cannot happen if you use the limits on the March coefficient, either to be above (needles) or below (plate) 1.0.
