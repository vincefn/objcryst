/* 
* ObjCryst++ : a Crystallographic computing library in C++
*			http://objcryst.sourceforge.net
*			http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/
*
*  (c) 2000-2001 Vincent FAVRE-NICOLIN vincefn@users.sourceforge.net
*
*/
/*   General.h
*  header file for various global declaration. Most of these declarations
* should be moved to other headers, and this file could eventually disappear.
*
*/

#ifndef _VFN_OBJCRYST_H_
#define _VFN_OBJCRYST_H_


//#include <stdlib.h>
#include <string>
//#include <iomanip>
#include <cmath>
//#include <typeinfo>
//#include <fstream>
//#include <ctime>

//profiling
#ifdef __MWERKS__
#include <Profiler.h>
#endif

#ifdef PROFILING_ON
#include "tau/include/Profile/Profiler.h"
#else

#define TAU_PROFILE(name, type, group) 
#define TAU_PROFILE_START(var)
#define TAU_PROFILE_TIMER(var, name, type, group)
#define TAU_PROFILE_STOP(var)
#define TAU_PROFILE_INIT(argc, argv)
#define TAU_PROFILE_SET_NODE(node)
#define TAU_PROFILE_SET_CONTEXT(context)
#define TAU_EVENT(event, data)
#define TAU_REPORT_STATISTICS()
#define TAU_REPORT_THREAD_STATISTICS()

#endif

#include "Quirks/VFNDebug.h"

using namespace std;
/** The namespace which includes all objects (crystallographic and
* algorithmic) in ObjCryst++.
*
* \note It may be a good idea to separate in 3 namespaces for ObjCryst, WXObjCryst,
* and RefinableObj
*/
namespace ObjCryst
{

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define DEG2RAD (M_PI/180.)
#define RAD2DEG (180./M_PI)

// forward declarations... These should be removed...
   class AsymmetricUnit;
   //class Atom;
   class Crystal ;
   class IOCrystTag;
   class PowderPattern;
   class RefinableObj;
   class RefinableObjClock;
   template<class T> class ObjRegistry;
   class RefinablePar;
   class Scatterer ;
   class ScatteringComponent;
   class ScatteringComponentList;
   class ScatteringData ;
   class ScatteringPower;
   class ScatteringPowerAtom;
   class SpaceGroup;
   class ZPolyhedron ;
   class ZScatterer ;
   #ifdef __WX__CRYST__
   class WXCrystObj;
   class WXRefinableObj;
   class WXCrystRegistry;
   #endif
/*
*/

//######################################################################

/// Type of radiation used. Only neutrons and X-Rays are used so far,
/// electrons would require a very different treatment.
enum RadiationType { RAD_NEUTRON, RAD_XRAY, RAD_ELECTRON};
/// Sample type (not used yet)
enum SampleType { SAMPLE_SINGLE_CRYSTAL, SAMPLE_POWDER};
/// Incident beam characteristics : monochromatic, X-Ray tube with Alpha1 and alpha2,
/// MAD (a few wavelengths-UNUSED YET), DAFS (continuous wavelength range-UNUSED YET)
/// LAUE (UNUSED YET)
enum WavelengthType { WAVELENGTH_MONOCHROMATIC, WAVELENGTH_ALPHA12, WAVELENGTH_MAD,
                      WAVELENGTH_DAFS , WAVELENGTH_LAUE};
/// Profile type for powder (could it be used fopr single crystals on 2D detectors ?)
enum ReflectionProfileType { PROFILE_GAUSSIAN, PROFILE_LORENTZIAN, PROFILE_PSEUDO_VOIGT,
                             PROFILE_PSEUDO_VOIGT_FINGER_COX_JEPHCOAT,
                             PROFILE_PEARSON_VII };

enum PowderBackgroundInterpType{ POWDER_BACKGROUND_LINEAR,
                                 POWDER_BACKGROUND_CUBIC_SPLINE};

#define XRAY_WAVELENGTH_TO_ENERGY 12398.4

//######################################################################
//  Exception.
/** \brief Exception class for ObjCryst++ library
*
*/
//######################################################################
//:TODO: This should go into VFNDebug.h
class ObjCrystException
{
   public:
      ObjCrystException();
      ObjCrystException(const string & message);
      ~ObjCrystException();
   protected:
   private:
};

}//Namespace

#endif //_VFN_OBJCRYST_H_
