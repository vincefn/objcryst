#ifndef _OBJCRYST_IOCRYST_H_
#define _OBJCRYST_IOCRYST_H_
#endif //_OBJCRYST_IOCRYST_H_

#include "RefinableObj/IO.h"
#include "RefinableObj/RefinableObj.h"
namespace ObjCryst
{
/** \brief Save the complete Objcryst environment.
*
* This saves all Crystal, PowderPattern, DiffDataSingleCrystal and GlobalOptimObj objects,
* using the global registries for these classes
*/
void XMLCrystFileSaveGlobal(const string & filename);
/** \brief Get the list (tags) of ObjCryst objects in a file
*
* This will recognize only certain tags in the file (Crystal,PowderPattern,
* DiffDataSingleCrystal, GlobalOptimObj). Eventually it should include also
* the ZScatterer objects.
*
* \note It will be the duty of the caller to destroy all the tags which have been
* allocated.
*/
ObjRegistry<XMLCrystTag> XMLCrystFileLoadObjectList(const string & filename);

/** \brief Load an object from a file, identifying it from its tag
*
* \param file: the filename from which the object will be loaded.
* \param tagName: the name of the tag
* \param name: the name of the object to be found (in a 'Name' attribute)
* \param obj: the pointer to the object to be loaded. The allocation will be done
* by the function, and the pointer changed accordingly.
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



void IOCrystFileSaveGlobal(const string & filename);
ObjRegistry<IOCrystTag> IOCrystFileLoadObjectList(const string & filename);
template<class T> void IOCrystFileLoadObject(const string & file,const IOCrystTag &tag, T*obj);
void IOCrystFileLoadAllObject(const string & file);

}
