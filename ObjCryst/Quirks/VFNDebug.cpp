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
/*
*  source file for LibCryst++  LibCryst Debug messages functions
*
*/

#ifdef __DEBUG__

#include "Quirks/VFNDebug.h"
//#include <stdlib.h>
#include <iostream.h>
int gDebugMessageGlobalLevel=10;
int gDebugMessageLevel=gDebugMessageLevel;
unsigned int gVFNDebugMessageIndent=0;
/*
void LibCrystDebugMessage(const string &message, const int level=0)
{  
   if(level >= sDebugMessageLevel) cout << "DEBUG MSG:" << message <<endl;
}
*/
void LibCrystDebugGlobalLevel(const int level)
{
   gDebugMessageGlobalLevel=level;
   LibCrystDebugLocalLevel(-1);
};

//Use this for a local modification of debug level messages. Call this at the
//beginning of the function, and call it without argument at the ned of the function
//to revert to the default global debug level.
void LibCrystDebugLocalLevel(const int level)
{
   if(level != -1) gDebugMessageLevel=level; 
   else gDebugMessageLevel=gDebugMessageGlobalLevel;
   cout << "DEBUG MSG: Setting debug level to " << gDebugMessageLevel <<endl;
}

int gVFNDebugMessageGlobalLevel=10;//The default : few messages.
int gVFNDebugMessageLevel=gVFNDebugMessageGlobalLevel;

#endif //__DEBUG__

