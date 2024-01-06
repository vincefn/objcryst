.. _compile_linux:

Compiling Fox under Linux
=========================

Requirements
------------
To install F.O.X. under Linux, you will need:
* standard development packages (gcc,..)
* OpenGL development libraries, e.g. `libopengl-dev`

You will also need more specific dependencies:
* wxGTK (>=3.0.2)
* newmat
* freeglut
* fftw3 (development library)

Not that the latter can be automatically donwloaded and linked statically as part of the Fox installation,
but you can use the system libraries when they are available.

Debian/Ubuntu
^^^^^^^^^^^^^
Install the following packages (22.04): ``gcc g++ libwxgtk3.0-gtk3 ibwxgtk3.0-gtk3-dev freeglut3 freeglut3-dev``

Download Fox
------------
Fox release from github
^^^^^^^^^^^^^^^^^^^^^^^
Get the Fox.tar.bz2 from `GitHub <https://github.com/vincefn/objcryst/releases>`_,
and uncompress:

.. code-block:: shell

  tar -xjf Fox-VERSION.tar.bz2

Get the latest Fox from the git repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: shell

  git clone https://github.com/vincefn/objcryst.git
  mv objcryst Fox
  cd Fox/ObjCryst
  ln -sf rules-gnu.mak rules.mak
  cd ../Fox
  ln -sf gnu.mak Makefile

Compile & install Fox
---------------------

To compile and install Fox, go to the Fox/Fox subdirectory (or Fox-VERSION/Fox/), and then type:

.. code-block:: shell

  make shared=1  # See note below about shared libraries options
  sudo make install #(for this last step you must be root).

**Note 1**: there are several options for the ```make``` command:

If fftw, wxgtk (>=3.0.2), freeglut or newmat are already available as shared libraries on your computer,
you can use them without recompiling them. You can just append ```shared-wxgtk=1```, ```shared-glut=1```,
``shared-fftw=1``` and ```shared-newmat=1```.

If you have all five shared libraries (as possible for Debian and Ubuntu), just use ```shared=1```. Examples:

.. code-block:: shell

  make shared-wxgtk=1 shared-glut=1 shared-fftw=1

  make shared=1


Fox will be installed in /usr/local/bin/Fox. Otherwise the compiled file is in the `src/Fox` subdirectory.

**Note 2**: if during the early stages of the compilation (when freeglut is compiled), you get
an error message about ```libGL.la```, this can mean that ```libGL.la``` is not in the same
directory as other x11 libraries. On my computer, where ```libGL.la``` is provided by the nVidia
installer, I needed to do "```ln -s /usr/lib/libGL.la /usr/X11R6/lib/```".

Compile & install Fox without GUI (Fox-nogui)
---------------------------------------------

It is possible to compile a version of Fox that does not use a GUI and can only be used
from the command-line. This is useful if you want to use Fox on a cluster.
To do that just you just need to compile the ```Fox-nogui``` target:

.. code-block:: shell

  tar -xjf Fox-VERSION.tar.bz2
  cd Fox-VERSION
  cd Fox
  make Fox-nogui shared=1

The resulting application is ```src/Fox-nogui```

Note that when switching from building ```Fox``` to building ```Fox-nogui``` (and vice-versa)
you must do a ``` make tidy``` first in the build directory.

Compiling an optimized version of Fox
-------------------------------------
*Note: the instructions below are most likely obsolete as they were written > 15 years ago, they are
left just in case they can be usefu. Feel free to experiment and propose an update !*
Compile options
^^^^^^^^^^^^^^^
The first set of optimizations can be activated by using processor-specific optimizations.
If you edit the ```Fox/ObjCryst/rules.mak``` file, and search for the part where the ```CPPFLAGS``` are defined:

.. code-block:: shell

  CPPFLAGS = -O3 -w -ffast-math -fstrict-aliasing -pipe -fomit-frame-pointer -funroll-loops

*Auto-vectorization*: starting with gcc 4.0.0, it is possible to automatically **vectorize** some loops:
append ```-ftree-vectorize``` to the ```CPPFLAGS``` options to do that
(see the examples commented out in the ```rules.mak``` file).

Profile driven optimizations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Fox can be further optimized by making "test runs" which are used to give hints to the compiler
on how to best optimize the code. To do this you (i) compile Fox by enabling the "recording"
of the optimization, then (ii) you run a few optimizations, then (iii) you recompile using the
recorded profile. To do that from the Fox subdirectory, do:

.. code-block:: shell

  make clean
  make Fox profile=2
  src/Fox --nogui example/pbso4-joint.xml --randomize -n 50000 -o /dev/null
  src/Fox --nogui example/Cimetidine-powder.xml --randomize -n 50000 -o /dev/null
  src/Fox --speedtest
  make clean
  make Fox profile=3
  make install

This yields about 10% faster code.

If you also launch the Fox GUI and do a profile-fitting (between the 'profile=2' and 'profile=3' compilations),
this will also **accelerate profile fitting  and least squares operations**.
