#------------------------------------------------------------------------------
# INFO
#------------------------------------------------------------------------------
Summary: F.O.X., Free Objects for Crystallography, a program for crystal structure analysis, and ab initio structure determination from powder diffraction
Name: ObjCryst-Fox
Version: 1.7.0
Release: 20070211
Copyright: GPL
Packager: Vincent Favre-Nicolin <vincefn@users.sourceforge.net>
Group: Sciences/Chemistry
Source: Fox.tar.bz2
#Patch: ObjCryst-Fox.patch
BuildRoot: %{_tmppath}/%{name}-buildroot
URL: http://objcryst.sourceforge.net/Fox/
#Provides: (what does it provide - binaries and/or libraries)
Prefix: /usr/local
DOCDIR: /usr/local/doc

#------------------------------------------------------------------------------
# REQUIRES (DISTRIB-DEPENDANT)
#------------------------------------------------------------------------------
#MANDRIVA 2007
Requires: libgtk+2.0_0 libmesaglu1 libmesagl1
BuildRequires: gcc gcc-c++ gcc-cpp libxorg-x11-devel libgtk+2.0_0-devel libmesaglu1-devel libmesagl1-devel
#FEDORA Core 3,4
#Requires: xorg-x11-Mesa-libGL xorg-x11-Mesa-libGLU
#BuildRequires:gcc gcc-c++ xorg-x11-devel

#------------------------------------------------------------------------------
# DESCRIPTION
#------------------------------------------------------------------------------
%description
FOX is an open-source program for crystallographers who wish to solve crystal 
structures from powder diffraction data (X-ray, neutron, neutron TOF).

Additionally, it can be used to create crystal structures and display them 
in 3D, simulate powder patterns and calculate structure factors. It can also
import crystal structures from CIF (Crystallographic Information File) files.

FOX is built on top of an object-oriented crystallographic library (ObjCryst++),
which can be used for other purposes than structure determination.

#------------------------------------------------------------------------------
# BUILD
#------------------------------------------------------------------------------
%prep
%setup -q -n Fox
#%patch -p1 -b .buildroot

%build
make RPM_OPT_FLAGS="$RPM_OPT_FLAGS"

%install
[ "$RPM_BUILD_ROOT" != "/" ] && rm -rf ${RPM_BUILD_ROOT}%{prefix}
mkdir -p ${RPM_BUILD_ROOT}%{prefix}/bin
#mkdir -p ${RPM_BUILD_ROOT}%{prefix}/man/man1

install -s -m 755 src/Fox ${RPM_BUILD_ROOT}%{prefix}/bin/Fox

%clean
[ "${RPM_BUILD_ROOT}%{prefix}" != "/" ] && rm -rf ${RPM_BUILD_ROOT}%{prefix}

%files
%defattr(-,root,root)
%doc README LICENSE ChangeLog
/usr/local/bin/Fox

%changelog
* Sun Feb 11 2007 Vincent Favre-Nicolin <vincefn@users.sf.net>
- Release 1.7

* Thu Jul 06 2005 Vincent Favre-Nicolin <vincefn@users.sf.net>
- First rpm package !

