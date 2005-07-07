#------------------------------------------------------------------------------
# INFO
#------------------------------------------------------------------------------
Summary: F.O.X., Free Objects for Crystallography, a program for the ab initio structure determination from powder diffraction
Name: ObjCryst-Fox
Version: 1.7.0
Release: CVS20050706
Copyright: GPL
Packager: Vincent Favre-Nicolin <vincefn@users.sourceforge.net
Group: Sciences/Physics
Source: Fox.tar.bz2
#Patch: ObjCryst-Fox.patch
BuildRoot: /var/tmp/%{name}-buildroot
URL: http://objcryst.sourceforge.net/cgi-bin/moin.cgi
#Provides: (what does it provide - binaries and/or libraries)
Prefix: /usr/local
DOCDIR: /usr/local/doc

#------------------------------------------------------------------------------
# REQUIRES (DISTRIB-DEPENDANT)
#------------------------------------------------------------------------------
#MANDRIVA 10.2
Requires: libwxgtk2.6 libwxgtkgl2.6 wxGTK2.6 libwxgtkgl2.6 libMesaGLU1 libMesaglut3 
BuildRequires: gcc gcc-c++ gcc-cpp gcc-g77 libwxgtk2.6-devel libMesaGLU1-devel libMesaglut3-devel libxorg-x11-devel
#FEDORA Core 3,4
#Requires: wxGTK2>=2.4.0 wxGTK2-gl>=2.4.0 wxGTK2-stc>=2.4.0 wxGTK2-xrc>=2.4.0 xorg-x11-Mesa-libGL xorg-x11-Mesa-libGLU freeglut
#BuildRequires:gcc gcc-c++ gcc-g77 wxGTK2-devel>=2.4.0 freeglut-devel xorg-x11-devel

#------------------------------------------------------------------------------
# DESCRIPTION
#------------------------------------------------------------------------------
%description
FOX is an open-source program for Crystallographers who wish to solve crystal structures
from powder diffraction data (X-ray, neutron, neutron TOF).

Additionally, it can be used to create crystal structures and display them in 3D,
simulate powder patterns and calculate structure factors.

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
%doc README LICENSE ChangeLog html
/usr/local/bin/Fox

%changelog
* Thu Jul 06 2005 Vincent Favre-Nicolin <vincefn@users.sf.net>
- First rpm package !

