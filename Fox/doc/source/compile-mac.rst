.. _compile_macos:

Compiling Fox under macOS
=========================

Requirements
------------
* XCode

Other requirements (wxWidgets, fftw, cctbx, newmat) will be automatically downloaded and compiled.

Compile Fox
-----------

* First download the linux package for Fox ```Fox-xxxx.tar.bz2``` from
  `GitHub <https://github.com/vincefn/objcryst/releases>`_, and uncompress it.
* Go into the Fox/Fox subdirectory, then compile Fox using ```make -f macosx.mak dist```
* You will get a disk image (\*.dmg) including the Fox application
