.. _manual_tips:

Fox tips
========
* you can optimize powder patterns with multiple phases. Just add several crystal phases to the
  PowderPattern object. There can be several unknown phases (if you solve several unknown phases at the same
  time, be sure to tell us !).

* you can use several powder patterns to solve one crystal structure. Just add all the powder
  patterns (X-ray, neutron, neutron TOF) to the Monte-Carlo object.

* the *Dynamical Occupancy Correction* can automatically correct atoms overlapping due to special
  positions or shared atoms between polyhedra. It is activated by default as an option in the Crystal object.
  However, if no overlap or special position is expected, do not activate it as it slows down a lot the computation.

* keep in mind that Fox is **open-source**: you can access and modify the source code if you want !
