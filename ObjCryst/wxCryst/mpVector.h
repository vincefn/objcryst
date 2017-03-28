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
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	FileName:	mpVector.h
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu	or  mikepolyakov@hotmail.com
//	Website	:	www.angelfire.com/linux/myp
//	Date	:	7/16/2002
//	
//		Provides basic vector handling. 
//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef MPVECTOR_H
#define MPVECTOR_H

#include <math.h>
#include <stdio.h>
#include <wx/string.h>

class mp4Vector;

class mpVector
{
friend class mp4Vector;
public:
	float x, y, z;
	
	mpVector();
	mpVector(float xx, float yy, float zz);
	mpVector(const mpVector& other);

	mpVector& Normalize();
	float Magnitude();
	mpVector Cross(const mpVector& other);
	mpVector operator - (mpVector v);
	mpVector operator + (mpVector v);
	float operator * (mpVector v);
	mpVector operator - (float c);
	mpVector operator + (float c);
	mpVector operator / (float c);
	mpVector operator * (float c);
	void operator = (const mpVector& other);
	operator mp4Vector() const;

	operator char*()  const;
	
};

class mp4Vector
{
public:
	float x, y, z, val;

	mp4Vector();
	mp4Vector(float aa, float bb, float cc, float dd);
	mp4Vector(const mp4Vector& other);
	mp4Vector(const mpVector& v, const float value);
	
	void operator = (const mp4Vector& v);
	void operator = (const mpVector& v);

	bool operator == (const mp4Vector& v) const;

	operator mpVector() const;
};
#endif
