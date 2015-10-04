/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2002 Vincent Favre-Nicolin vincefn@users.sourceforge.net
        2000-2001 University of Geneva (Switzerland)

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
/*   IO.h
*
*/
#ifndef _OBJCRYST_IOCRYST_H_
#define _OBJCRYST_IOCRYST_H_
#endif //_OBJCRYST_IOCRYST_H_

#include "ObjCryst/RefinableObj/IO.h"
#include "ObjCryst/RefinableObj/RefinableObj.h"


namespace ObjCryst
{
/** \brief Save all Objcryst++ objects.
*
* This saves all Crystal, PowderPattern, DiffDataSingleCrystal and GlobalOptimObj objects,
* using the global registries for these classes. All other objects (Scatterer,
* ScatteringPower, PowderPatternComponent are saved as well since they are sub-objects
* of Crystal or PowderPattern objects).
*
* Saving is done in well-formed xml format.
*/
void XMLCrystFileSaveGlobal(const string & filename);
/** \brief Save all Objcryst++ objects.
*
* This saves all Crystal, PowderPattern, DiffDataSingleCrystal and GlobalOptimObj objects,
* using the global registries for these classes. All other objects (Scatterer,
* ScatteringPower, PowderPatternComponent are saved as well since they are sub-objects
* of Crystal or PowderPattern objects).
*
* Saving is done in well-formed xml format.
*/
void XMLCrystFileSaveGlobal(std::ostream &out);
/** \brief Get the list (tags) of ObjCryst objects in a file
*
* This will recognize only certain tags in the file (Crystal,PowderPattern,
* DiffDataSingleCrystal, GlobalOptimObj). Eventually it should include also
* the ZScatterer objects.
*
* \note It will be the duty of the caller to destroy all the tags which have been
* allocated.
*
* NOT TESTED YET !
*/
ObjRegistry<XMLCrystTag> XMLCrystFileLoadObjectList(const string & filename);

/** \brief Load an object from a file, identifying it from its tag
*
* \param file: the filename from which the object will be loaded.
* \param tagName: the name of the tag
* \param name: the name of the object to be found (in a 'Name' attribute)
* \param obj: the pointer to the object to be loaded. The allocation will be done
* by the function, and the pointer changed accordingly.
*
* NOT TESTED YET !
*/
template<class T> void XMLCrystFileLoadObject(const string & file,
                                              const string &tagName,
                                              const string &name, T*obj);

/** \brief Load all 'top' objects from a file (Crystal, PowderPattern, DiffDataSingleCrystal
*  and GlobalOptimObj objects). All objects are directly allocated, and can be accessed through
* their respective global registry (eg gCrystalRegistry fro a Crysta, etc...)
*
* \param file: the filename from which the objects will be loaded.
*/
void XMLCrystFileLoadAllObject(const string & file);
/** \brief Load all 'top' objects from a file (Crystal, PowderPattern, DiffDataSingleCrystal
*  and GlobalOptimObj objects). All objects are directly allocated, and can be accessed through
* their respective global registry (eg gCrystalRegistry fro a Crysta, etc...)
*
* \param file: the filename from which the objects will be loaded.
*/
void XMLCrystFileLoadAllObject(std::istream &is);

}
