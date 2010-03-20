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
*  source file for LibCryst++ ObjCrystException class
*
*/

#include "ObjCryst/ObjCryst/General.h"
#include <iostream>
#include <ctime>
#include "ObjCryst/ObjCryst/IO.h"

namespace ObjCryst
{

bool ObjCrystException::verbose = true;

ObjCrystException::ObjCrystException() : message()
{
   if (ObjCrystException::verbose)
   {
      cout << "LibCryst ++ exception thrown!!" << endl;
   }
}

ObjCrystException::ObjCrystException(const string & _message)
{

   if (!ObjCrystException::verbose)
   {
       message = _message;
       return;
   }

   static bool inException;
   cout << "LibCryst ++ exception thrown!!" << endl;
   cout << "  Message: " + message <<endl;
   if(false==inException)
   {
      inException=true;
      string saveFileName="ObjCryst";
      time_t date=time(0);
      char strDate[40];
      strftime(strDate,sizeof(strDate),"%Y-%m-%d_%H-%M-%S",gmtime(&date));//%Y-%m-%dT%H:%M:%S%Z
      saveFileName=saveFileName+strDate+".xml";
      cout << "Attempting to save ObjCryst++ environment to file:"<<saveFileName<<endl;
      try
      {
         XMLCrystFileSaveGlobal(saveFileName);
      }
      catch(...)
      {cout<<"Sorry, failed to save ObjCryst++ environment"<<endl;}
      inException=false;
   }
}

ObjCrystException::~ObjCrystException(){}

//######################################################################
void (*fpObjCrystInformUser)(const string &)=ObjCrystInformUserStdOut;

void ObjCrystInformUserStdOut(const string &str)
{
   cout <<str<<endl;
}

}//namespace


