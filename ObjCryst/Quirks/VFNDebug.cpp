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

