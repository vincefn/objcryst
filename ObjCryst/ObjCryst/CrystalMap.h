/* 
* ObjCryst++ : a Crystallographic computing library in C++
*			http://objcryst.sourceforge.net
*			http://www.ccp14.ac.uk/ccp/web-mirrors/objcryst/
*
*  (c) 2000-2001 Vincent FAVRE-NICOLIN vincefn@users.sourceforge.net
*
*/
/*   CrystalMap.h header file for the CrystalMap object
*
*/

#include "ObjCryst/Crystal.h"

namespace ObjCryst
{

/** Base class for Crystal maps (electronic density, Patterson,...)
* This is an abstract base class.
*
*
*/
class CrystalMap
{
	public:
		/// Constructor, specifying the associated Crystal Structure
		CrystalMap(const Crystal &cryst);
		/// Constructor, specifying the associated Crystal Structure
		/// and the number of grid points along all directions
		CrystalMap(const Crystal &cryst,
					  const unsigned int nx,
					  const unsigned int ny,
					  const unsigned int nz);
		/// Destructor
		virtual ~CrystalMap();
		/// Set the number of grid points along all directions
		void SetNbGridPoint(const unsigned int nx,
					   		  const unsigned int ny,
					   		  const unsigned int nz);
		/// Get the map as a 3d array
		const CrystArray3D_REAL& GetMap3D()const;
		/// Get a section of the map 
		const CrystArray3D_REAL& GetMap2Dxy(const REAL z)const;
		/// Get a section of the map 
		const CrystArray3D_REAL& GetMap2Dxz(const REAL y)const;
		/// Get a section of the map 
		const CrystArray3D_REAL& GetMap2Dyz(const REAL x)const;
		
		/// Get a value at a given point (with grid coordinates)
		/// Note that the map is periodic,
		/// so any coordinate in 3d space can be inputted.
		REAL GetValueGrid(const int  x, const int  y, const int  z)const;
		/// Get a value at a given point (with cartesian coordinates in Angstroems)
		/// Note that the map is periodic,
		/// so any coordinate in 3d space can be inputted.
		REAL GetValueCart(const REAL x, const REAL y, const REAL z)const;
		/// Get a value at a given point (with fractionnal coordinates)
		/// Note that the map is periodic,
		/// so any coordinate in 3d space can be inputted.
		REAL GetValueFrac(const REAL x, const REAL y, const REAL z)const;
		/// Init the OpenGL display list (displays the map as a mesh)
		void GLInitDisplayList(const REAL xMin=-.1,const REAL xMax=1.1,
                             const REAL yMin=-.1,const REAL yMax=1.1,
                             const REAL zMin=-.1,const REAL zMax=1.1);
	
	private:
		/// Get the list of polygons to display
		void MarchingTetrahedra()const;
		/// Actual computation of the map.
		virtual void CalcMap()=0;
		/// The Crystal associated with this map
		const Crystal* mpCrystal;
		/// The number of points along each direction
		unsigned int mNbPointX,mNbPointY,mNbPointZ;
}
