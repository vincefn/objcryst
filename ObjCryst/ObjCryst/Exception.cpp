/*
* ObjCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*
*  source file for LibCryst++ ObjCrystException class
*
*/

#include "ObjCryst/General.h"
#include <iostream>
#include <ctime>
#include <ObjCryst/IO.h>

namespace ObjCryst
{

ObjCrystException::ObjCrystException()
{
   cout << "LibCryst ++ exception thrown!!" << endl;
}

ObjCrystException::ObjCrystException(const string & message)
{
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


