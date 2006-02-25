/////////////////////////////////////////////////////////////////////////////////////////////
//	FileName:	MarchingCubes.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//	website	:	www.angelfire.com/linux/myp
//	date	:	July 2002
//	
//	Description:	'Straight' and Recursive Marching Cubes Algorithms
//				Normal vectors are defined for each vertex as a gradients
//				For function definitions see MarchingCubes.h
//				For a tutorial on Marching Cubes please visit www.angelfire.com/myp/linux
//
//	Please email me with any suggestion/bugs.
/////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include "wxCryst/mpVector.h"

//struct for storing triangle information - 3 vertices and 3 normal vectors for each vertex
typedef struct {
	mpVector p[3];
	mpVector norm[3];
} TRIANGLE;

//does Linear Interpolation between points p1 and p2 (they already contain their computed values)
mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value);

////////////////////////////////////////////////////////////////////////////////////////
//POINTERS TO FUNCTIONS
//pointer to function which computes the value at point p
typedef float (*FORMULA)(mpVector);
////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////
///// MARCHING CUBES ALGORITHM /////
/////////////////////////////////////////////////////////////////////////////////////////

// 'STRAIGHT' Marching Cubes Algorithm //////////////////////////////////////////////////
//takes number of cells (ncellsX, ncellsY, ncellsZ) to subdivide on each axis
// minValue used to pass into LinearInterp
// gradFactor for each axis (multiplies each component of gradient vector by 1/(2*gradFactor) ).
//		Should be set to the length of a side (or close to it)
// array of length (ncellsX+1)(ncellsY+1)(ncellsZ+1) of mp4Vector points containing coordinates and values
//returns pointer to triangle array and the number of triangles in numTriangles
//note: array of points is first taken on z axis, then y and then x. So for example, if you iterate through it in a
//       for loop, have indexes i, j, k for x, y, z respectively, then to get the point you will have to make the
//		 following index: i*(ncellsY+1)*(ncellsZ+1) + j*(ncellsZ+1) + k .
//		Also, the array starts at the minimum on all axes.
TRIANGLE* MC(int ncellsX, int ncellsY, int ncellsZ, float gradFactorX, float gradFactorY, float gradFactorZ,
										float minValue, mp4Vector * points, int &numTriangles);
/////////////////////////////////////////////////////////////////////////////////////////


// RECURSIVE Marching Cubes	/////////////////////////////////////////////////////////////
//DrawBacks: must know how many objects are drawn

//This function starts the Recursive Marching Cubes
// numCubes: number of intersecting cubes passed as starting points
// ii, jj, kk: arrays of size numCubes. Contain respective indexes of cubes that are intersected
TRIANGLE* MarchingCubesRec(int ncellsX, int ncellsY, int ncellsZ, 
							float gradFactorX, float gradFactorY, float gradFactorZ,
							int numCubes, int *ii, int *jj, int *kk, 
							float minValue, mp4Vector * points, int &numTriangles);

//Next 6 functions are called by the corresponding face of each cube (e.g. face 0 calls MCFace0 etc...)
// Each function accepts the following information as arguments:
//	First 3 arguments: number of cells to subdivide on each axis
//	gradFactor: factor to scale the gradient vectors (multiplies by 1/(2*gradFactor) )
//	ind: index of this cube in the points array (this is so it doesnt have to be computed again)
//	i,j,k: indexes of the this cube on x,y,z axis respectively
//	minValue: minimum used in LinearInterp
//	points: array of type mp4Vector and size (ncellsX+1)*(ncellsY+1)*(ncellsZ+1) that contains
//		coordinates and values at each point
//	triangle: array of triangles which is being built and returned at the end of recursion
//	numTriangles: number of triangles formed
//	prevVerts: adjacent 4 vertices of the previous cube, passed in the special order which the correspoding
//		MCFace function recognizes. For specificatino on which indexes passed from the last cube go to 
//		which vertexes of the current cube visit my website at www.angelfire.com/linux/myp
//	prevIntVerts: array of 4 linearly interpolated vertices on 4 adjacent edges
//	edgeIndex: integer, bits of which correspond to the edges of the current cube which have been computed
//		from the previous cube
//	prevGradVerts: array of 4 gradient vectors at 4 adjacent vertices
//	prevGrads: linearly interpolated gradient vectors on 4 adjacent edges
//	gradIndex: integer bits of which correspond to already computed vertices of the current cube
//	marchedCubes: bool array of size ncellsX*ncellsY*ncellsZ which stores which cubes already have been marched
//		through. Initialized to FALSE. This is not the most effective way in terms of memory management, but
//		faster than using STLs vector<bool> which stores booleans as bits.
//* NOTE *: Each of these 6 functions run marching cubes on the 'next'cube adjacent to the surface number
//		that is specified in their names (for numbering see my webpage). Each assigns the previous computed
//		values to the face 'opposite' of that specified. Then each one runs Marching Cubes for the current cube.
//		It returns doing nothing if the cube is not intersected, however, if it is then other recursive
//		functions are called using macros MC_FACE defined for each face. See MarchingCubes.cpp
//		For example MCFace0 will initialize surface 2 using the previous values that are passed to it. Then if
//		the current cube is intersected it will run MC_FACE0, MC_FACE1, MC_FACE3, MC_FACE4, and MC_FACE5. Each
//		returnes the pointer to the array of formed triangles.

//FACE 0 Marching Cubes
TRIANGLE* MCFace0(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 1 Marching Cubes
TRIANGLE* MCFace1(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 2 Marching Cubes
TRIANGLE* MCFace2(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 3 Marching Cubes
TRIANGLE* MCFace3(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 4 Marching Cubes
TRIANGLE* MCFace4(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);
//FACE 5 Marching Cubes
TRIANGLE* MCFace5(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 									
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex,
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes);



//Does Marching Cubes on cube (i, j, k) in points
// Requirements: needs index ind which specifies which cube it is (its 0th point in points array)
//	If cube is intersected the triangles are added to triangles array passed to it, and numTriangles is
//	incremented appropriately. verts, array of vertices of this cube should be initialized.
//	intVerts, interpolated vertices on edges, gradVerts, gradients on vertices, and grads, interpolated
//	gradient vectors on edges dont have to be initialized, but if they are corresponding bits in edgeIndex
//	and indGrad index should show it. (for numbering see my web page)
//	Global variables YtimeZ should be initialized to (ncellsZ+1)*(ncellsY+1) and pointsZ should be ncellsZ+1
TRIANGLE* MarchOneCube(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector verts[8], mpVector intVerts[12], int &edgeIndex, 
						mp4Vector gradVerts[8], mpVector grads[12], int &indGrad);


//Find the first intersecting cube (if it exists).
//Starts at (i,j,k)=(0,0,0) at the minimum x,y,z and then goes linearly searching for an intersecting cube.
//Returns array of indices i,j,k of the first found cube or -1,-1,-1 if cube was not found
float* MCFind(int ncellsX, int ncellsY, int ncellsZ, float minValue, mp4Vector * points);

//Calls the MCFind and if found calls MarchingCubesRec with the returned indices
//If not found NULL is returned and numTriangles is set to ZERO.
TRIANGLE* MCRecFind(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						float minValue, mp4Vector * points, int &numTriangles);

#endif
