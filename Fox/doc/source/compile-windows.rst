.. _compile_windows:

Compiling Fox on windows (windows 10 and visual studio 2019 instructions)
=========================================================================
* Get **Visual Studio 2019** from https://visualstudio.microsoft.com/fr/downloads/.
  Install it with C++ desktop tools, include CLI (command=line interface), C++ clang
* Download **objcryst** from git (https://github.com/vincefn/objcryst), rename the root 'objcryst' directory to 'Fox'
* Download **cctbx.tar.bz2** and **newmat.tar.bz2** (from https://github.com/vincefn/objcryst/releases/tag/v2021-3rdPartyLibs) 
  in the Fox directory and uncompress them (e.g. with 7-zip)
* Delete the Fox/cctbx/include/boost directory
* Download **boost** 1.74 (7z format) from https://sourceforge.net/projects/boost/files/boost/1.74.0/.
  Unzip in Fox and rename the directory to 'boost'
* Open the visual studio command prompt. In the Fox/boost directory, run ```bootstrap``` and then
  ```b2 --build_dir=..\boost_build toolset=msvc release threading=multi --build_type=complete```.
* Next download the **wxwidgets 3.1.4** installer and run it with default installation,
  select the installation directory as Fox/wxWidgets
* Open the visual studio command prompt. In Fox\wxWidgets\build\msw, run
  ```nmake /f makefile.vc BUILD=release RUNTIME_LIBS=static``` as per https://github.com/wxWidgets/wxWidgets/blob/v3.1.2/docs/msw/install.md
* Download **fftw-3.3.5-dll32.zip** from http://www.fftw.org/install/windows.html, uncompress it as Fox/fftw
* In a visual studio command prompt, in the Fox/fftw directory, execute ```lib /def:libfftw3f-3.def```
* Open **Fox_vc12.sln** from Fox/Fox/, and then you can build Fox in debug or release version
