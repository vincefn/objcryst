/*  This program is free software; you can redistribute it and/or modify
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
/////////////////////////////////////////////////////////////////////////////////////////////
//	FileName:	MarchingCubes.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu  or  mikepolyakov@hotmail.com
//	website	:	www.angelfire.com/linux/myp
//	date	:	July 2002
//	
//	Description:	'Straight' and Recursive Marching Cubes Algorithms
//					Recursive method is faster than the 'straight' one, especially when intersection does not 
//						have to be searched for every time.
//				Normal vectors are defined for each vertex as a gradients
//				For function definitions see MarchingCubes.h
//				For a tutorial on Marching Cubes please visit www.angelfire.com/myp/linux
//
//	Please email me with any suggestion/bugs.
/////////////////////////////////////////////////////////////////////////////////////////////

#include "MC.h"
#include "MCTable.h"
#include <math.h>

//Linear Interpolation between two points
mpVector LinearInterp(mp4Vector p1, mp4Vector p2, float value)
{
	mpVector p;
	if(fabs(p1.val - p2.val) > 0.00001)
		p = (mpVector)p1 + ((mpVector)p2 - (mpVector)p1)/(p2.val - p1.val)*(value - p1.val);
	else 
		p = (mpVector)p1;
	return p;
}


//Macros used to compute gradient vector on each vertex of a cube
//argument should be the name of array of vertices
//can be verts or *verts if done by reference
#define CALC_GRAD_VERT_0(verts) mp4Vector(points[ind-YtimeZ].val-(verts[1]).val,points[ind-pointsZ].val-(verts[4]).val,points[ind-1].val-(verts[3]).val, (verts[0]).val);
#define CALC_GRAD_VERT_1(verts) mp4Vector((verts[0]).val-points[ind+2*YtimeZ].val,points[ind+YtimeZ-pointsZ].val-(verts[5]).val,points[ind+YtimeZ-1].val-(verts[2]).val, (verts[1]).val);
#define CALC_GRAD_VERT_2(verts) mp4Vector((verts[3]).val-points[ind+2*YtimeZ+1].val,points[ind+YtimeZ-ncellsZ].val-(verts[6]).val,(verts[1]).val-points[ind+YtimeZ+2].val, (verts[2]).val);
#define CALC_GRAD_VERT_3(verts) mp4Vector(points[ind-YtimeZ+1].val-(verts[2]).val,points[ind-ncellsZ].val-(verts[7]).val,(verts[0]).val-points[ind+2].val, (verts[3]).val);
#define CALC_GRAD_VERT_4(verts) mp4Vector(points[ind-YtimeZ+ncellsZ+1].val-(verts[5]).val,(verts[0]).val-points[ind+2*pointsZ].val,points[ind+ncellsZ].val-(verts[7]).val, (verts[4]).val);
#define CALC_GRAD_VERT_5(verts) mp4Vector((verts[4]).val-points[ind+2*YtimeZ+ncellsZ+1].val,(verts[1]).val-points[ind+YtimeZ+2*pointsZ].val,points[ind+YtimeZ+ncellsZ].val-(verts[6]).val, (verts[5]).val);
#define CALC_GRAD_VERT_6(verts) mp4Vector((verts[7]).val-points[ind+2*YtimeZ+ncellsZ+2].val,(verts[2]).val-points[ind+YtimeZ+2*ncellsZ+3].val,(verts[5]).val-points[ind+YtimeZ+ncellsZ+3].val, (verts[6]).val);
#define CALC_GRAD_VERT_7(verts) mp4Vector(points[ind-YtimeZ+ncellsZ+2].val-(verts[6]).val,(verts[3]).val-points[ind+2*ncellsZ+3].val,(verts[4]).val-points[ind+ncellsZ+3].val, (verts[7]).val);

///////////////////////////////////////////////////////////////////////////////////////////////////////
// GLOBAL //
//Global variables - so they dont have to be passed into functions
int pointsZ;	//number of points on Z zxis (equal to ncellsZ+1)
int YtimeZ;		//'plane' of cubes on YZ (equal to (ncellsY+1)*pointsZ )
///////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
//	'STRAIGHT' MARCHING CUBES	ALGORITHM  ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
//for gradients at the edges values 1.0, 1.0, 1.0, 1.0  are given
TRIANGLE* MC(int ncellsX, int ncellsY, int ncellsZ, 
						float gradFactorX, float gradFactorY, float gradFactorZ,
						float minValue, mp4Vector * points, int &numTriangles)
{
	//this should be enough space, if not change 3 to 4
	TRIANGLE * triangles = new TRIANGLE[3*ncellsX*ncellsY*ncellsZ];
	numTriangles = int(0);
	
	pointsZ = ncellsZ+1;			//initialize global variable (for extra speed) 
	YtimeZ = (ncellsY+1)*pointsZ;	
	int lastX = ncellsX;			//left from older version
	int lastY = ncellsY;
	int lastZ = ncellsZ;

	mp4Vector *verts[8];			//vertices of a cube (array of pointers for extra speed)
	mpVector intVerts[12];			//linearly interpolated vertices on each edge
	int cubeIndex;					//shows which vertices are outside/inside
	int edgeIndex;					//index returned by edgeTable[cubeIndex]
	mp4Vector gradVerts[8];			//gradients at each vertex of a cube		
	mpVector grads[12];				//linearly interpolated gradients on each edge
	int indGrad;					//shows which gradients already have been computed
	int ind, ni, nj;				//ind: index of vertex 0
	//factor by which corresponding coordinates of gradient vectors are scaled
	mpVector factor(1.0/(2.0*gradFactorX), 1.0/(2.0*gradFactorY), 1.0/(2.0*gradFactorZ));

	//MAIN LOOP: goes through all the points
	for(int i=0; i < lastX; i++) {			//x axis
		ni = i*YtimeZ;
		for(int j=0; j < lastY; j++) {		//y axis
			nj = j*pointsZ;
			for(int k=0; k < lastZ; k++, ind++)	//z axis
			{
				//initialize vertices
				ind = ni + nj + k;
				verts[0] = &points[ind];
				verts[1] = &points[ind + YtimeZ];
				verts[4] = &points[ind + pointsZ];
				verts[5] = &points[ind + YtimeZ + pointsZ];
				verts[2] = &points[ind + YtimeZ + 1];
				verts[3] = &points[ind + 1];
				verts[6] = &points[ind + YtimeZ + pointsZ + 1];
				verts[7] = &points[ind + pointsZ + 1];

				//get the index
				cubeIndex = int(0);
				for(int n=0; n < 8; n++)
					if(verts[n]->val <= minValue) cubeIndex |= (1 << n);
				
				//check if its completely inside or outside
				if(!edgeTable[cubeIndex]) continue;
			
				indGrad = int(0);
				edgeIndex = edgeTable[cubeIndex];
				
				if(edgeIndex & 1) {
					intVerts[0] = LinearInterp(*verts[0], *verts[1], minValue);
					if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
					else gradVerts[0] = mp4Vector(1.0, 1.0, 1.0, 1.0);	//for now do not wrap around
					if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
					else gradVerts[1] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					indGrad |= 3;
					grads[0] = LinearInterp(gradVerts[0], gradVerts[1], minValue);
					grads[0].x *= factor.x; grads[0].y *= factor.y; grads[0].z *= factor.z;
				}
				if(edgeIndex & 2) {
					intVerts[1] = LinearInterp(*verts[1], *verts[2], minValue);
					if(! (indGrad & 2)) {
						if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
						else gradVerts[1] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 2;
					}
					if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
					else gradVerts[2] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					indGrad |= 4;
					grads[1] = LinearInterp(gradVerts[1], gradVerts[2], minValue);
					grads[1].x *= factor.x; grads[1].y *= factor.y; grads[1].z *= factor.z;
				}
				if(edgeIndex & 4) {
					intVerts[2] = LinearInterp(*verts[2], *verts[3], minValue);
					if(! (indGrad & 4)) {
						if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
						else gradVerts[2] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 4;
					}
					if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
					else gradVerts[3] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					indGrad |= 8;
					grads[2] = LinearInterp(gradVerts[2], gradVerts[3], minValue);
					grads[2].x *= factor.x; grads[2].y *= factor.y; grads[2].z *= factor.z;
				}
				if(edgeIndex & 8) {
					intVerts[3] = LinearInterp(*verts[3], *verts[0], minValue);
					if(! (indGrad & 8)) {
						if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
						else gradVerts[3] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 8;
					}
					if(! (indGrad & 1)) {
						if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
						else gradVerts[0] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 1;
					}
					grads[3] = LinearInterp(gradVerts[3], gradVerts[0], minValue);
					grads[3].x *= factor.x; grads[3].y *= factor.y; grads[3].z *= factor.z;
				}
				if(edgeIndex & 16) {
					intVerts[4] = LinearInterp(*verts[4], *verts[5], minValue);
					
					if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
					else gradVerts[4] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					
					if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
					else gradVerts[5] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					
					indGrad |= 48;
					grads[4] = LinearInterp(gradVerts[4], gradVerts[5], minValue);
					grads[4].x *= factor.x; grads[4].y *= factor.y; grads[4].z *= factor.z;
				}
				if(edgeIndex & 32) {
					intVerts[5] = LinearInterp(*verts[5], *verts[6], minValue);
					if(! (indGrad & 32)) {
						if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
						else gradVerts[5] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 32;
					}
					
					if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
					else gradVerts[6] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					indGrad |= 64;
					grads[5] = LinearInterp(gradVerts[5], gradVerts[6], minValue);
					grads[5].x *= factor.x; grads[5].y *= factor.y; grads[5].z *= factor.z;
				}
				if(edgeIndex & 64) {
					intVerts[6] = LinearInterp(*verts[6], *verts[7], minValue);
					if(! (indGrad & 64)) {
						if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
						else gradVerts[6] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 64;
					}
					
					if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
					else gradVerts[7] = mp4Vector(1.0, 1.0, 1.0, 1.0);
					indGrad |= 128;
					grads[6] = LinearInterp(gradVerts[6], gradVerts[7], minValue);
					grads[6].x *= factor.x; grads[6].y *= factor.y; grads[6].z *= factor.z;
				}
				if(edgeIndex & 128) {
					intVerts[7] = LinearInterp(*verts[7], *verts[4], minValue);
					if(! (indGrad & 128)) {
						if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
						else gradVerts[7] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 128;
					}
					if(! (indGrad & 16)) {
						if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
						else gradVerts[4] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 16;
					}
					grads[7] = LinearInterp(gradVerts[7], gradVerts[4], minValue);
					grads[7].x *= factor.x; grads[7].y *= factor.y; grads[7].z *= factor.z;
				}
				if(edgeIndex & 256) {
					intVerts[8] = LinearInterp(*verts[0], *verts[4], minValue);
					if(! (indGrad & 1)) {
						if(i != 0 && j != 0 && k != 0) gradVerts[0] = CALC_GRAD_VERT_0(*verts)
						else gradVerts[0] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 1;
					}
					if(! (indGrad & 16)) {
						if(i != 0 && j != lastY-1 && k != 0) gradVerts[4] = CALC_GRAD_VERT_4(*verts)
						else gradVerts[4] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 16;
					}
					grads[8] = LinearInterp(gradVerts[0], gradVerts[4], minValue);
					grads[8].x *= factor.x; grads[8].y *= factor.y; grads[8].z *= factor.z;
				}
				if(edgeIndex & 512) {
					intVerts[9] = LinearInterp(*verts[1], *verts[5], minValue);
					if(! (indGrad & 2)) {
						if(i != lastX-1 && j != 0 && k != 0) gradVerts[1] = CALC_GRAD_VERT_1(*verts)
						else gradVerts[1] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 2;
					}
					if(! (indGrad & 32)) {
						if(i != lastX-1 && j != lastY-1 && k != 0) gradVerts[5] = CALC_GRAD_VERT_5(*verts)
						else gradVerts[5] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 32;
					}
					grads[9] = LinearInterp(gradVerts[1], gradVerts[5], minValue);
					grads[9].x *= factor.x; grads[9].y *= factor.y; grads[9].z *= factor.z;
				}
				if(edgeIndex & 1024) {
					intVerts[10] = LinearInterp(*verts[2], *verts[6], minValue);
					if(! (indGrad & 4)) {
						if(i != lastX-1 && j != 0 && k != 0) gradVerts[2] = CALC_GRAD_VERT_2(*verts)
						else gradVerts[5] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 4;
					}
					if(! (indGrad & 64)) {
						if(i != lastX-1 && j != lastY-1 && k != lastZ-1) gradVerts[6] = CALC_GRAD_VERT_6(*verts)
						else gradVerts[6] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 64;
					}
					grads[10] = LinearInterp(gradVerts[2], gradVerts[6], minValue);
					grads[10].x *= factor.x; grads[10].y *= factor.y; grads[10].z *= factor.z;
				}
				if(edgeIndex & 2048) {
					intVerts[11] = LinearInterp(*verts[3], *verts[7], minValue);
					if(! (indGrad & 8)) {
						if(i != 0 && j != 0 && k != lastZ-1) gradVerts[3] = CALC_GRAD_VERT_3(*verts)
						else gradVerts[3] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 8;
					}
					if(! (indGrad & 128)) {
						if(i != 0 && j != lastY-1 && k != lastZ-1) gradVerts[7] = CALC_GRAD_VERT_7(*verts)
						else gradVerts[7] = mp4Vector(1.0, 1.0, 1.0, 1.0);
						indGrad |= 128;
					}
					grads[11] = LinearInterp(gradVerts[3], gradVerts[7], minValue);
					grads[11].x *= factor.x; grads[11].y *= factor.y; grads[11].z *= factor.z;
				}

				//now build the triangles using triTable
				for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
					int index[3] = {triTable[cubeIndex][n+2], triTable[cubeIndex][n+1], triTable[cubeIndex][n]};
					for(int h=0; h < 3; h++) {	//copy vertices and normals into triangles array
						triangles[numTriangles].p[h] = intVerts[index[h]];
						triangles[numTriangles].norm[h] = grads[index[h]];
					}
					numTriangles++;	//one more triangle has been added
				}			
			}	//END OF FOR LOOP ON Z AXIS
		}	//END OF FOR LOOP ON Y AXIS
	}	//END OF FOR LOOP ON X AXIS

	//free all wasted space
	TRIANGLE * retTriangles = new TRIANGLE[numTriangles];
	for(int i=0; i < numTriangles; i++)
		retTriangles[i] = triangles[i];
	delete [] triangles;

	return retTriangles;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	RECURSIVE MARCHING CUBES ALGORITHM ////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  MACROS  ///////////////////////////////////////////////////////////////////////////////////////////////////
//these macros initialize data and then run marching cubes on the cube with the surface having the specified 
//	number as 'recieving' data (but number of that surface for the current cube is going to be 'opposite').
//	Each runs the corresponding recursive function
//	For numbering, to see which indices of prevVerts,... correspong to indices of the current cube, see
//	my webpage at www.angelfire.com/linux/myp

#define MC_FACE0																							\
{																											\
  if(! marchedCubes[ind - 1]) {																				\
	passGradIndex = 0;																						\
	if(gradIndex & 1) passGradIndex |= 8;																	\
	if(gradIndex & 2) passGradIndex |= 4;																	\
	if(gradIndex & 32) passGradIndex |= 64;																	\
	if(gradIndex & 16) passGradIndex |= 128;																\
	passEdgeIndex = 0;																						\
	if(edgeIndex & 1) passEdgeIndex |= 4;																	\
	if(edgeIndex & 512) passGradIndex |= 1024;																\
	if(edgeIndex & 16) passEdgeIndex |= 64;																	\
	if(edgeIndex & 256) passGradIndex |= 2048;																\
	prevVerts[0] = verts[0]; prevVerts[1] = verts[1]; prevVerts[2] = verts[5]; prevVerts[3] = verts[4];		\
	prevIntVerts[0] = intVerts[0]; prevIntVerts[1] = intVerts[9];											\
	prevIntVerts[2] = intVerts[4]; prevIntVerts[3] = intVerts[8];											\
	prevGradVerts[0] = gradVerts[0]; prevGradVerts[1] = gradVerts[1];										\
	prevGradVerts[2] = gradVerts[5]; prevGradVerts[3] = gradVerts[4];										\
	prevGrads[0] = grads[0]; prevGrads[1] = grads[9]; prevGrads[2] = grads[4]; prevGrads[3] = grads[8];		\
	triangles = MCFace0(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind-1, i, j, k-1, minValue, points, triangles, numTriangles,					\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex, marchedCubes);							\
  }																											\
}

#define MC_FACE1																							\
{																											\
  if(! marchedCubes[ind + YtimeZ]) {																		\
	passGradIndex = 0;																						\
	if(gradIndex & 4) passGradIndex |= 8;																	\
	if(gradIndex & 2) passGradIndex |= 1;																	\
	if(gradIndex & 32) passGradIndex |= 16;																	\
	if(gradIndex & 64) passGradIndex |= 128;																\
	passEdgeIndex = 0;																						\
	if(edgeIndex & 2) passEdgeIndex |= 8;																	\
	if(edgeIndex & 512) passEdgeIndex |= 256;																\
	if(edgeIndex & 32) passEdgeIndex |= 128;																\
	if(edgeIndex & 1024) passEdgeIndex |= 2048;																\
	prevVerts[0] = verts[2]; prevVerts[1] = verts[1]; prevVerts[2] = verts[5]; prevVerts[3] = verts[6];		\
	prevIntVerts[0] = intVerts[1]; prevIntVerts[1] = intVerts[9];											\
	prevIntVerts[2] = intVerts[5]; prevIntVerts[3] = intVerts[10];											\
	prevGradVerts[0] = gradVerts[2]; prevGradVerts[1] = gradVerts[1];										\
	prevGradVerts[2] = gradVerts[5]; prevGradVerts[3] = gradVerts[6];										\
	prevGrads[0] = grads[1]; prevGrads[1] = grads[9]; prevGrads[2] = grads[5]; prevGrads[3] = grads[10];	\
	triangles = MCFace1(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind+YtimeZ, i+1, j, k, minValue, points, triangles, numTriangles,				\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex,  marchedCubes);						\
  }																											\
}

#define MC_FACE2																							\
{																											\
  if(! marchedCubes[ind + 1]) {																				\
	passGradIndex = 0;																						\
	if(gradIndex & 8) passGradIndex |= 1;																	\
	if(gradIndex & 4) passGradIndex |= 2;																	\
	if(gradIndex & 64) passGradIndex |= 32;																	\
	if(gradIndex & 128) passGradIndex |= 16;																\
	passEdgeIndex = 0;																						\
	if(edgeIndex & 4) passEdgeIndex |= 1;																	\
	if(edgeIndex & 1024) passEdgeIndex |= 512;																\
	if(edgeIndex & 64) passEdgeIndex |= 16;																	\
	if(edgeIndex & 2048) passEdgeIndex |= 256;																\
	prevVerts[0] = verts[3]; prevVerts[1] = verts[2]; prevVerts[2] = verts[6]; prevVerts[3] = verts[7];		\
	prevIntVerts[0] = intVerts[2]; prevIntVerts[1] = intVerts[10];											\
	prevIntVerts[2] = intVerts[6]; prevIntVerts[3] = intVerts[11];											\
	prevGradVerts[0] = gradVerts[3]; prevGradVerts[1] = gradVerts[2];										\
	prevGradVerts[2] = gradVerts[6]; prevGradVerts[3] = gradVerts[7];										\
	prevGrads[0] = grads[2]; prevGrads[1] = grads[10]; prevGrads[2] = grads[6]; prevGrads[3] = grads[11];	\
	triangles = MCFace2(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind+1, i, j, k+1, minValue, points, triangles, numTriangles,					\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex, marchedCubes);							\
  }																											\
}

#define MC_FACE3																							\
{																											\
  if(! marchedCubes[ind - YtimeZ]) {																		\
	passGradIndex = 0;																						\
	if(gradIndex & 8) passGradIndex |= 4;																	\
	if(gradIndex & 1) passGradIndex |= 2;																	\
	if(gradIndex & 128) passGradIndex |= 64;																\
	if(gradIndex & 16) passGradIndex |= 32;																	\
	passEdgeIndex = 0;																						\
	if(edgeIndex & 8) passEdgeIndex |= 2;																	\
	if(edgeIndex & 256) passEdgeIndex |= 512;																\
	if(edgeIndex & 128) passEdgeIndex |= 32;																\
	if(edgeIndex & 2048) passEdgeIndex |= 1024;																\
	prevVerts[0] = verts[3]; prevVerts[1] = verts[0]; prevVerts[2] = verts[4]; prevVerts[3] = verts[7];		\
	prevIntVerts[0] = intVerts[3]; prevIntVerts[1] = intVerts[8];											\
	prevIntVerts[2] = intVerts[7]; prevIntVerts[3] = intVerts[11];											\
	prevGradVerts[0] = gradVerts[3]; prevGradVerts[1] = gradVerts[0];										\
	prevGradVerts[2] = gradVerts[4]; prevGradVerts[3] = gradVerts[7];										\
	prevGrads[0] = grads[3]; prevGrads[1] = grads[8]; prevGrads[2] = grads[7]; prevGrads[3] = grads[11];	\
	triangles = MCFace3(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind-YtimeZ, i-1, j, k, minValue, points, triangles, numTriangles,				\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex, marchedCubes);							\
  }																											\
}

//done
#define MC_FACE4																							\
{																											\
  if(! marchedCubes[ind + pointsZ]) {																		\
	passGradIndex = 0;																						\
	if(gradIndex & 128) passGradIndex |= 8;																	\
	if(gradIndex & 64) passGradIndex |= 4;																	\
	if(gradIndex & 32) passGradIndex |= 2;																	\
	if(gradIndex & 16) passGradIndex |= 1;																	\
	passEdgeIndex = 0;																						\
	if(edgeIndex & 128) passEdgeIndex |= 8;																	\
	if(edgeIndex & 64) passEdgeIndex |= 4;																	\
	if(edgeIndex & 32) passEdgeIndex |= 2;																	\
	if(edgeIndex & 16) passEdgeIndex |= 1;																	\
	prevVerts[0] = verts[7]; prevVerts[1] = verts[6]; prevVerts[2] = verts[5]; prevVerts[3] = verts[4];		\
	prevIntVerts[0] = intVerts[6]; prevIntVerts[1] = intVerts[5];											\
	prevIntVerts[2] = intVerts[4]; prevIntVerts[3] = intVerts[7];											\
	prevGradVerts[0] = gradVerts[7]; prevGradVerts[1] = gradVerts[6];										\
	prevGradVerts[2] = gradVerts[5]; prevGradVerts[3] = gradVerts[4];										\
	prevGrads[0] = grads[6]; prevGrads[1] = grads[5]; prevGrads[2] = grads[4]; prevGrads[3] = grads[7];		\
	triangles = MCFace4(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind+pointsZ, i, j+1, k, minValue, points, triangles, numTriangles,				\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex, marchedCubes);							\
  }																											\
}

#define MC_FACE5																							\
{																											\
  if(! marchedCubes[ind - ncellsZ - 1]) {																	\
	passGradIndex = (gradIndex << 4) & 240;																	\
	passEdgeIndex = (edgeIndex << 4) & 240;																	\
	prevVerts[0] = verts[3]; prevVerts[1] = verts[2]; prevVerts[2] = verts[1]; prevVerts[3] = verts[0];		\
	prevIntVerts[0] = intVerts[2]; prevIntVerts[1] = intVerts[1];											\
	prevIntVerts[2] = intVerts[0]; prevIntVerts[3] = intVerts[3];											\
	prevGradVerts[0] = gradVerts[3]; prevGradVerts[1] = gradVerts[2];										\
	prevGradVerts[2] = gradVerts[1]; prevGradVerts[3] = gradVerts[0];										\
	prevGrads[0] = grads[2]; prevGrads[1] = grads[1]; prevGrads[2] = grads[0]; prevGrads[3] = grads[3];		\
	triangles = MCFace5(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,					\
							ind-ncellsZ-1, i, j-1, k, minValue, points, triangles, numTriangles,			\
							prevVerts, prevIntVerts, passEdgeIndex,											\
							prevGradVerts, prevGrads, passGradIndex, marchedCubes);							\
  }																											\
}
/// END FACE MACROS /////////////////////////////////////////////////////////////////////////////////////////





// RECURSIVE Marching Cubes Function - cubes at indexes ii, jj, kk intersect the surface
//	Number of intersecting cubes = numCubes
TRIANGLE* MarchingCubesRec(int ncellsX, int ncellsY, int ncellsZ, 
							float gradFactorX, float gradFactorY, float gradFactorZ,
							int numCubes, int *ii, int *jj, int *kk, 
							float minValue, mp4Vector * points, int &numTriangles)
{
	TRIANGLE * triangles = new TRIANGLE[3*ncellsX*ncellsY*ncellsZ];
	numTriangles = int(0);
	//check if none of the starting points are on the outside perimeter
	for(int n=0; n < numCubes; n++) {
		if(ii[n] <= 0 || ii[n] >= ncellsX-1) return NULL;
		if(jj[n] <= 0 || jj[n] >= ncellsY-1) return NULL;
		if(kk[n] <= 0 || kk[n] >= ncellsZ-1) return NULL;
	}
	//array stores which cubes have been marched through - initialized to FALSE
	int all = ncellsX*ncellsY*ncellsZ;
	bool* marchedCubes = new bool[all];
	for(int i=0; i < all; i++) marchedCubes[i] = FALSE;		//initialize
		
	mp4Vector verts[8];				//vertices of a starting cube
	mpVector intVerts[12];			//linearly interpolated vertices on each edge
	int edgeIndex;					//shows which edges are intersected
	mp4Vector gradVerts[8];			//gradients at each vertex of a cube		
	mpVector grads[12];				//linearly interpolated gradients on each edge
	int gradIndex;					//show on which vertices gradients have been computed
	//initialize global variables - for speed - these would be used by all the recursive functions
	pointsZ = ncellsZ+1;
	YtimeZ = (ncellsY+1)*pointsZ;
	//arrays used to pass already computed values from this initial cube to the next ones
	mp4Vector prevVerts[4];			//passes vertices
	mpVector prevIntVerts[4];		//passes interpolated vertices on edges
	mp4Vector prevGradVerts[4];		//passes gradients at vertices
	mpVector prevGrads[4];			//passes interpolated gradients on edges
	
	//two new indexes formed for each face
	int passEdgeIndex, passGradIndex;	//used to tell which vertices and which edges have been initialized
	//indices
	int i, j, k;
	//initialize first cubes and 'launch' the recursion for each
	for(int n=0; n < numCubes; n++) {
		//init vertices
		i = ii[n]; j = jj[n]; k = kk[n];
		int ind = i*YtimeZ + j*pointsZ + k;
		verts[0] = points[ind];
		verts[1] = points[ind + YtimeZ];
		verts[2] = points[ind + YtimeZ + 1];
		verts[3] = points[ind + 1];
		verts[4] = points[ind + pointsZ];
		verts[5] = points[ind + YtimeZ + pointsZ];
		verts[6] = points[ind + YtimeZ + pointsZ + 1];
		verts[7] = points[ind + pointsZ + 1];
		
		//first check if this cube wasnt marched in recursive calls of the previous cube
		if(! marchedCubes[ind]) {
			//run marching cubes on the initial cube
			gradIndex = edgeIndex = 0;
			triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,
							ind, i, j, k, minValue, points, triangles, numTriangles, 
							verts, intVerts, edgeIndex, gradVerts, grads, gradIndex);
			//this cube has been done:
			marchedCubes[ind] = TRUE;
			//run M.C. on all 6 faces
			MC_FACE0
			MC_FACE1
			MC_FACE2
			MC_FACE3
			MC_FACE4
			MC_FACE5
		}
	}
	//free wasted space
	TRIANGLE * retTriangles = new TRIANGLE[numTriangles];
	for(int i=0; i < numTriangles; i++) retTriangles[i] = triangles[i];
	delete [] triangles;
	delete [] marchedCubes;
	return retTriangles;	//done
}

//SURFACE 0 - Cube ran on surface 0 of previous cube. Recieving side: 2.
TRIANGLE* MCFace0(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = points[ind];
	verts[1] = points[ind + YtimeZ];
	verts[2] = prevVerts[1];
	verts[3] = prevVerts[0];
	verts[4] = points[ind + ncellsZ + 1];
	verts[5] = points[ind + YtimeZ + ncellsZ + 1];
	verts[6] = prevVerts[2];
	verts[7] = prevVerts[3];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[2] = prevIntVerts[0]; intVerts[10] = prevIntVerts[1];
	intVerts[6] = prevIntVerts[2]; intVerts[11] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[3] = prevGradVerts[0]; gradVerts[2] = prevGradVerts[1];
	gradVerts[6] = prevGradVerts[2]; gradVerts[7] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[2] = prevGrads[0]; grads[10] = prevGrads[1]; grads[6] = prevGrads[2]; grads[11] = prevGrads[3];
	//for test if this cube is intersected:
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face to be passed into other recursive functions
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE0
	MC_FACE1
	MC_FACE3
	MC_FACE4
	MC_FACE5
	return triangles;
}

//SURFACE 1 - Cube ran on surface 1 of previous cube. Recieving side: 3.
TRIANGLE* MCFace1(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = prevVerts[1];
	verts[1] = points[ind + YtimeZ];
	verts[2] = points[ind + YtimeZ + 1];
	verts[3] = prevVerts[0];
	verts[4] = prevVerts[2];
	verts[5] = points[ind + YtimeZ + ncellsZ + 1];
	verts[6] = points[ind + YtimeZ + ncellsZ + 2];
	verts[7] = prevVerts[3];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[3] = prevIntVerts[0]; intVerts[8] = prevIntVerts[1];
	intVerts[7] = prevIntVerts[2]; intVerts[11] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[3] = prevGradVerts[0]; gradVerts[0] = prevGradVerts[1];
	gradVerts[4] = prevGradVerts[2]; gradVerts[7] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[3] = prevGrads[0]; grads[8] = prevGrads[1]; grads[7] = prevGrads[2]; grads[11] = prevGrads[3];
	//for test if this cube is intersected:
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face to be passed into other recursive functions
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE0
	MC_FACE1
	MC_FACE2
	MC_FACE4
	MC_FACE5
	return triangles;
}

//SURFACE 2 - Cube ran on surface 2 of previous cube. Recieving side: 0.
TRIANGLE* MCFace2(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = prevVerts[0];
	verts[1] = prevVerts[1];
	verts[2] = points[ind + YtimeZ + 1];
	verts[3] = points[ind + 1];
	verts[4] = prevVerts[3];
	verts[5] = prevVerts[2];
	verts[6] = points[ind + YtimeZ + ncellsZ + 2];
	verts[7] = points[ind + ncellsZ + 2];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[0] = prevIntVerts[0]; intVerts[9] = prevIntVerts[1];
	intVerts[4] = prevIntVerts[2]; intVerts[8] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[0] = prevGradVerts[0]; gradVerts[1] = prevGradVerts[1];
	gradVerts[5] = prevGradVerts[2]; gradVerts[4] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[0] = prevGrads[0]; grads[9] = prevGrads[1]; grads[4] = prevGrads[2]; grads[8] = prevGrads[3];
	//for test if this cube is intersected
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE1
	MC_FACE2
	MC_FACE3
	MC_FACE4
	MC_FACE5
	return triangles;
}

//SURFACE 3 - Cube ran on surface 3 of previous cube. Recieving side: 1.
TRIANGLE* MCFace3(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = points[ind];
	verts[1] = prevVerts[1];	
	verts[2] = prevVerts[0];	
	verts[3] = points[ind + 1];
	verts[4] = points[ind + ncellsZ + 1];
	verts[5] = prevVerts[2];
	verts[6] = prevVerts[3];
	verts[7] = points[ind + ncellsZ + 2];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[1] = prevIntVerts[0]; intVerts[9] = prevIntVerts[1];
	intVerts[5] = prevIntVerts[2]; intVerts[10] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[2] = prevGradVerts[0]; gradVerts[1] = prevGradVerts[1];
	gradVerts[5] = prevGradVerts[2]; gradVerts[6] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[1] = prevGrads[0]; grads[9] = prevGrads[1]; grads[5] = prevGrads[2]; grads[10] = prevGrads[3];
	//for test if this cube is intersected
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE0
	MC_FACE2
	MC_FACE3
	MC_FACE4
	MC_FACE5
	return triangles;
}

//SURFACE 4 - Cube ran on surface 4 of previous cube. Recieving side: 5.
TRIANGLE* MCFace4(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = prevVerts[3];
	verts[1] = prevVerts[2];	
	verts[2] = prevVerts[1];	
	verts[3] = prevVerts[0];
	verts[4] = points[ind + ncellsZ + 1];
	verts[5] = points[ind + YtimeZ + ncellsZ + 1];
	verts[6] = points[ind + YtimeZ + ncellsZ + 2];
	verts[7] = points[ind + ncellsZ + 2];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[2] = prevIntVerts[0]; intVerts[1] = prevIntVerts[1];
	intVerts[0] = prevIntVerts[2]; intVerts[3] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[3] = prevGradVerts[0]; gradVerts[2] = prevGradVerts[1];
	gradVerts[1] = prevGradVerts[2]; gradVerts[0] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[2] = prevGrads[0]; grads[1] = prevGrads[1]; grads[0] = prevGrads[2]; grads[3] = prevGrads[3];
	//for test if this cube is intersected
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE0
	MC_FACE1
	MC_FACE2
	MC_FACE3
	MC_FACE4
	return triangles;
}

//SURFACE 5 - Cube ran on surface 5 of previous cube. Recieving side: 4.
TRIANGLE* MCFace5(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector prevVerts[4], mpVector prevIntVerts[4], int edgeIndex, 
						mp4Vector prevGradVerts[4], mpVector prevGrads[4], int gradIndex, bool* marchedCubes)
{
	//first check if not outside the region
	if(i <= 0 || i >= ncellsX-1 || j <= 0 || j >= ncellsY-1 || k <= 0 || k >= ncellsZ-1)
		return triangles;
	//make sure we save that this cube was marched through
	marchedCubes[ind] = TRUE;
	//initialize vertices
	mp4Vector verts[8];
	verts[0] = points[ind];
	verts[1] = points[ind + YtimeZ];	
	verts[2] = points[ind + YtimeZ + 1];	
	verts[3] = points[ind + 1];
	verts[4] = prevVerts[3];
	verts[5] = prevVerts[2];
	verts[6] = prevVerts[1];
	verts[7] = prevVerts[0];
	//initialize edges from the last cube
	mpVector intVerts[12];
	intVerts[6] = prevIntVerts[0]; intVerts[5] = prevIntVerts[1];
	intVerts[4] = prevIntVerts[2]; intVerts[7] = prevIntVerts[3];
	//initialize gradients on vertices
	mp4Vector gradVerts[8];
	gradVerts[7] = prevGradVerts[0]; gradVerts[6] = prevGradVerts[1];
	gradVerts[5] = prevGradVerts[2]; gradVerts[4] = prevGradVerts[3];
	//initialize gradients on edges from the last cube
	mpVector grads[12];
	grads[6] = prevGrads[0]; grads[5] = prevGrads[1]; grads[4] = prevGrads[2]; grads[7] = prevGrads[3];
	//for test if this cube is intersected
	int oldNumTriangles = numTriangles;				
	//run marching cubes on this cube
	triangles = MarchOneCube(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ, ind, i, j, k,
								minValue, points, triangles, numTriangles, verts, intVerts, edgeIndex,
								gradVerts, grads, gradIndex);
	//check if this cube is intersected
	if(numTriangles == oldNumTriangles)
		return triangles;
	//two new indexes formed for each face
	int passEdgeIndex, passGradIndex;
	//run recursive functions on each surface
	MC_FACE0
	MC_FACE1
	MC_FACE2
	MC_FACE3
	MC_FACE5
	return triangles;
}

//Marching Cubes on a single cube (i, j, k)
//	Verts should be initialized before hand
//  Global variables YtimeZ and pointsZ should be defined and initialized
TRIANGLE* MarchOneCube(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						int ind, int i, int j, int k, 
						float minValue, mp4Vector * points, TRIANGLE *triangles, int &numTriangles,
						mp4Vector verts[8], mpVector intVerts[12], int &edgeIndex, 
						mp4Vector gradVerts[8], mpVector grads[12], int &indGrad)
{
	//factor by which gradients are scaled
	mpVector factor(1.0/(2.0*gradFactorX), 1.0/(2.0*gradFactorY), 1.0/(2.0*gradFactorZ));
	//get the index
	int cubeIndex = int(0);
	for(int n=0; n < 8; n++) if(verts[n].val <= minValue) cubeIndex |= (1 << n);
	//check if its completely inside or outside
	if(!edgeTable[cubeIndex]) return triangles;

	int prevEdgeIndex = edgeIndex;
	edgeIndex = edgeTable[cubeIndex];

	if((edgeIndex & 1) && !(prevEdgeIndex & 1)) {
		intVerts[0] = LinearInterp(verts[0], verts[1], minValue);
		if(! (indGrad & 1)) {
			gradVerts[0] = CALC_GRAD_VERT_0(verts)
			indGrad |= 1;
		}
		if(! (indGrad & 2)) {
			gradVerts[1] = CALC_GRAD_VERT_1(verts)
			indGrad |= 2;
		}
		grads[0] = LinearInterp(gradVerts[0], gradVerts[1], minValue);
		grads[0].x *= factor.x; grads[0].y *= factor.y; grads[0].z *= factor.z;
	}
	if((edgeIndex & 2) && !(prevEdgeIndex & 2)) {
		intVerts[1] = LinearInterp(verts[1], verts[2], minValue);
		if(! (indGrad & 2)) {
			gradVerts[1] = CALC_GRAD_VERT_1(verts)
			indGrad |= 2;
		}
		if(! (indGrad & 4)) {
			gradVerts[2] = CALC_GRAD_VERT_2(verts)
			indGrad |= 4;
		}
		grads[1] = LinearInterp(gradVerts[1], gradVerts[2], minValue);
		grads[1].x *= factor.x; grads[1].y *= factor.y; grads[1].z *= factor.z;
	}
	if(edgeIndex & 4) {
		intVerts[2] = LinearInterp(verts[2], verts[3], minValue);
		if(! (indGrad & 4)) {
			gradVerts[2] = CALC_GRAD_VERT_2(verts)
			indGrad |= 4;
		}
		if(! (indGrad & 8)) {
			gradVerts[3] = CALC_GRAD_VERT_3(verts)
			indGrad |= 8;
		}
		grads[2] = LinearInterp(gradVerts[2], gradVerts[3], minValue);
		grads[2].x *= factor.x; grads[2].y *= factor.y; grads[2].z *= factor.z;
	}
	if((edgeIndex & 8) && !(prevEdgeIndex & 8)) {
		intVerts[3] = LinearInterp(verts[3], verts[0], minValue);
		if(! (indGrad & 8)) {
			gradVerts[3] = CALC_GRAD_VERT_3(verts)
			indGrad |= 8;
		}
		if(! (indGrad & 1)) {
			gradVerts[0] = CALC_GRAD_VERT_0(verts)
			indGrad |= 1;
		}
		grads[3] = LinearInterp(gradVerts[3], gradVerts[0], minValue);
		grads[3].x *= factor.x; grads[3].y *= factor.y; grads[3].z *= factor.z;
	}
	if((edgeIndex & 16) && !(prevEdgeIndex & 16)) {
		intVerts[4] = LinearInterp(verts[4], verts[5], minValue);
		if(! (indGrad & 16)) {
			gradVerts[4] = CALC_GRAD_VERT_4(verts)
			indGrad |= 16;
		}
		if(! (indGrad & 32)) {
			gradVerts[5] = CALC_GRAD_VERT_5(verts)
			indGrad |= 32;
		}
		grads[4] = LinearInterp(gradVerts[4], gradVerts[5], minValue);
		grads[4].x *= factor.x; grads[4].y *= factor.y; grads[4].z *= factor.z;
	}
	if((edgeIndex & 32) && !(prevEdgeIndex & 32)) {
		intVerts[5] = LinearInterp(verts[5], verts[6], minValue);
		if(! (indGrad & 32)) {
			gradVerts[5] = CALC_GRAD_VERT_5(verts)
			indGrad |= 32;
		}
		if(! (indGrad & 64)) {
			gradVerts[6] = CALC_GRAD_VERT_6(verts)
			indGrad |= 64;
		}
		grads[5] = LinearInterp(gradVerts[5], gradVerts[6], minValue);
		grads[5].x *= factor.x; grads[5].y *= factor.y; grads[5].z *= factor.z;
	}
	if((edgeIndex & 64) && !(prevEdgeIndex & 64)) {
		intVerts[6] = LinearInterp(verts[6], verts[7], minValue);
		if(! (indGrad & 64)) {
			gradVerts[6] = CALC_GRAD_VERT_6(verts)
			indGrad |= 64;
		}
		if(! (indGrad & 128)) {
			gradVerts[7] = CALC_GRAD_VERT_7(verts)
			indGrad |= 128;
		}
		grads[6] = LinearInterp(gradVerts[6], gradVerts[7], minValue);
		grads[6].x *= factor.x; grads[6].y *= factor.y; grads[6].z *= factor.z;
	}
	if((edgeIndex & 128) && !(prevEdgeIndex & 128)) {
		intVerts[7] = LinearInterp(verts[7], verts[4], minValue);
		if(! (indGrad & 128)) {
			gradVerts[7] = CALC_GRAD_VERT_7(verts)
			indGrad |= 128;
		}
		if(! (indGrad & 16)) {
			gradVerts[4] = CALC_GRAD_VERT_4(verts)
			indGrad |= 16;
		}
		grads[7] = LinearInterp(gradVerts[7], gradVerts[4], minValue);
		grads[7].x *= factor.x; grads[7].y *= factor.y; grads[7].z *= factor.z;
	}
	if((edgeIndex & 256) && !(prevEdgeIndex & 256)) {
		intVerts[8] = LinearInterp(verts[0], verts[4], minValue);
		if(! (indGrad & 1)) {
			gradVerts[0] = CALC_GRAD_VERT_0(verts)
			indGrad |= 1;
		}
		if(! (indGrad & 16)) {
			gradVerts[4] = CALC_GRAD_VERT_4(verts)
			indGrad |= 16;
		}
		grads[8] = LinearInterp(gradVerts[0], gradVerts[4], minValue);
		grads[8].x *= factor.x; grads[8].y *= factor.y; grads[8].z *= factor.z;
	}
	if((edgeIndex & 512) && !(prevEdgeIndex & 512)) {
		intVerts[9] = LinearInterp(verts[1], verts[5], minValue);
		if(! (indGrad & 2)) {
			gradVerts[1] = CALC_GRAD_VERT_1(verts)
			indGrad |= 2;
		}
		if(! (indGrad & 32)) {
			gradVerts[5] = CALC_GRAD_VERT_5(verts)
			indGrad |= 32;
		}
		grads[9] = LinearInterp(gradVerts[1], gradVerts[5], minValue);
		grads[9].x *= factor.x; grads[9].y *= factor.y; grads[9].z *= factor.z;
	}
	if((edgeIndex & 1024) && !(prevEdgeIndex & 1024)) {
		intVerts[10] = LinearInterp(verts[2], verts[6], minValue);
		if(! (indGrad & 4)) {
			gradVerts[2] = CALC_GRAD_VERT_2(verts)
			indGrad |= 4;
		}
		if(! (indGrad & 64)) {
			gradVerts[6] = CALC_GRAD_VERT_6(verts)
			indGrad |= 64;
		}
		grads[10] = LinearInterp(gradVerts[2], gradVerts[6], minValue);
		grads[10].x *= factor.x; grads[10].y *= factor.y; grads[10].z *= factor.z;
	}
	if((edgeIndex & 2048) && !(prevEdgeIndex & 2048)) {
		intVerts[11] = LinearInterp(verts[3], verts[7], minValue);
		if(! (indGrad & 8)) {
			gradVerts[3] = CALC_GRAD_VERT_3(verts)
			indGrad |= 8;
		}
		if(! (indGrad & 128)) {
			gradVerts[7] = CALC_GRAD_VERT_7(verts)
			indGrad |= 128;
		}
		grads[11] = LinearInterp(gradVerts[3], gradVerts[7], minValue);
		grads[11].x *= factor.x; grads[11].y *= factor.y; grads[11].z *= factor.z;
	}
	//now build the triangles using triTable and add them to the triangle array
	for (int n=0; triTable[cubeIndex][n] != -1; n+=3) {
		int index[3] = {triTable[cubeIndex][n+2], triTable[cubeIndex][n+1], triTable[cubeIndex][n]};
		for(int h=0; h < 3; h++) {
			triangles[numTriangles].p[h] = intVerts[index[h]];
			triangles[numTriangles].norm[h] = grads[index[h]];
		}
		numTriangles++;	//one more triangle...
	}
	return triangles;
}



//Finds first intersecting cube and returns its indices (or -1 if nothing is found)
float* MCFind(int ncellsX, int ncellsY, int ncellsZ, float minValue, mp4Vector * points)
{
	pointsZ = ncellsZ+1;			//initializes global variables
	YtimeZ = (ncellsY+1)*pointsZ;	
	int lastX = ncellsX - 1, lastY = ncellsY - 1, lastZ = ncellsZ - 1;	//for a little extra speed
	mp4Vector *verts[8];		//store address of each point rather than copying x,y,z, and val
	int cubeIndex, ind, ni;
	float *found = new float[3];//returned indices: initialized to -1
	found[0] = -1; found[0] = -1; found[0] = -1; 		

	for(int i=1; i < lastX; i++) {
		ni = i*YtimeZ;
		for(int j=1; j < lastY; j++) {
			//get the first 4 vertices (0, 1, 4, 5) and place them into verts at indices 3, 2, 7, 6
			ind = ni + j*pointsZ + 1;
			verts[3] = &points[ind];
			verts[2] = &points[ind + YtimeZ];
			verts[7] = &points[ind + pointsZ];
			verts[6] = &points[ind + YtimeZ + pointsZ];
			for(int k=1; k < lastZ; k++, ind++) {
				//initialize vertices
				verts[0] = verts[3];	//these are saved from the last iteration
				verts[1] = verts[2];
				verts[4] = verts[7];
				verts[5] = verts[6];
				verts[4] = &points[ind + pointsZ];
				verts[5] = &points[ind + YtimeZ + pointsZ];
				verts[6] = &points[ind + YtimeZ + ncellsZ + 2];
				verts[7] = &points[ind + ncellsZ + 2];	
				//build index - shows which vertices are outside/inside
				cubeIndex = int(0);
				for(int n=0; n < 8; n++) if(verts[n]->val <= minValue) cubeIndex |= (1 << n);
				if(edgeTable[cubeIndex]) {	//if found then return them
					found[0] = i; found[1] = j; found[2] = k;
					return found;
				}
			}
		}
	}
	return found;
}
		
//First finds intersecting cube and runs recursive Marching Cubes.
TRIANGLE* MCRecFind(int ncellsX, int ncellsY, int ncellsZ,
						float gradFactorX, float gradFactorY, float gradFactorZ,
						float minValue, mp4Vector * points, int &numTriangles)
{
	float *found = MCFind(ncellsX, ncellsY, ncellsZ, minValue, points);
	if(found[0] == -1) {	//if nothing is found return NULL
		numTriangles = 0;
		return NULL;
	}
	int i[1] = {(int)found[0]}, j[1] = {(int)found[1]}, k[1] = {(int)found[2]};
	delete [] found;
	return MarchingCubesRec(ncellsX, ncellsY, ncellsZ, gradFactorX, gradFactorY, gradFactorZ,
							1, i, j, k, minValue, points, numTriangles);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// THE END ///////////////////////////////////////////////////////////////
