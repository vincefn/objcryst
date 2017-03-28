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
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
//	FileName:	mpVector.cpp
//	Author	:	Michael Y. Polyakov
//	email	:	myp@andrew.cmu.edu	or  mikepolyakov@hotmail.com
//	Website	:	www.angelfire.com/linux/myp
//	Date	:	7/16/2002
//
//		Provides basic vector handling.
//////////////////////////////////////////////////////////////////////////////////////////////////////

#include "mpVector.h"

mpVector::mpVector(float xx, float yy, float zz) :
	x(xx), y(yy), z(zz)
	{ }

mpVector::mpVector() : x(0), y(0), z(0)
	{ }
mpVector::mpVector(const mpVector& other) : x(other.x), y(other.y), z(other.z)
	{ }

mpVector& mpVector::Normalize()
{
	float length = sqrt(x*x + y*y + z*z);
	if(!length) return *this;
	x /= length;
	y /= length;
	z /= length;
	return *this;
}

float mpVector::Magnitude()
{
	return sqrt(x*x + y*y + z*z);
}

mpVector mpVector::Cross(const mpVector& other)
{
	return mpVector(y*other.z-z*other.y, z*other.x-x*other.z, x*other.y-y*other.x);
}

mpVector mpVector::operator - (mpVector v)
{
	return mpVector(x - v.x, y - v.y, z - v.z);
}

mpVector mpVector::operator + (mpVector v)
{
	return mpVector(x + v.x, y + v.y, z + v.z);
}

float mpVector ::operator * (mpVector v)
{
	return x*v.x + y*v.y + z*v.z;
}
mpVector mpVector::operator - (float c)
{
	return mpVector(x-c, y-c, z-c);
}

mpVector mpVector::operator + (float c)
{
	return mpVector(x+c, y+c, z+c);
}

mpVector mpVector::operator / (float c)
{
	return mpVector(x/c, y/c, z/c);
}

mpVector mpVector::operator * (float c)
{
	return mpVector(x*c, y*c, z*c);
}

void mpVector::operator = (const mpVector& other)
{
	x = other.x;
	y = other.y;
	z = other.z;
}

mpVector::operator mp4Vector() const
{
	return mp4Vector(*this);
}

mpVector::operator char*()  const
{
	return (char*)wxString::Format(_T("(%f %f %f)"), x, y, z).char_str();
}





mp4Vector::mp4Vector() : x(0), y(0), z(0), val(0)
{ }

mp4Vector::mp4Vector(float aa, float bb, float cc, float dd) :
	x(aa), y(bb), z(cc), val(dd)
{ }

mp4Vector::mp4Vector(const mp4Vector& other) :
	x(other.x), y(other.y), z(other.z), val(other.val)
{ }

mp4Vector::mp4Vector(const mpVector& v, const float value) :
	x(v.x), y(v.x), z(v.z), val(value)
{ }

void mp4Vector::operator = (const mp4Vector& v)
{
	x = v.x; y = v.y; z = v.z; val = v.val;
}

void mp4Vector::operator = (const mpVector& v)
{
	x = v.x; y = v.y; z = v.z;
}

mp4Vector::operator mpVector() const
{
	return mpVector(x, y, z);
}

bool mp4Vector::operator == (const mp4Vector& v) const
{
	return x == v.x && y == v.y && z == v.z && val == v.val;
}
