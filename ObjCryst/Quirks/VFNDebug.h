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
/*   VFNDebug.h
*  header file for debugging objects & functions
*
*/
#ifndef __VFN_DEBUG__H__
#define __VFN_DEBUG__H__

#ifdef __DEBUG__

/// Set the Global debug level for messages
void LibCrystDebugGlobalLevel(const int level);

/// Use this for a local modification of debug level messages. Call this at the
/// beginning of the function, and call it without argument at the end of the function
/// to revert to the default global debug level.
void LibCrystDebugLocalLevel(const int level);

extern int gVFNDebugMessageGlobalLevel;
extern int gVFNDebugMessageLevel;
extern unsigned int gVFNDebugMessageIndent;

//Debug messages are printed only if __DEBUG__ is on, and if their level
//is greater or equal to debugMessageGlobalLevel
//  0 : messages from low-level routines, 
//  5
// 10 : messages from top LibCryst++ routines 

#define VFN_DEBUG_MESSAGE(message,level) \
   if(level >= gVFNDebugMessageLevel) \
   {\
      for(unsigned int iii=0;iii<gVFNDebugMessageIndent;iii++) cout <<"  ";\
      cout << "%DEBUG:"<< level << "  "\
      << message << " (at " << __FILE__ << "," << __LINE__ << ")" <<endl;\
   }
      
#define VFN_DEBUG_MESSAGE_SHORT(message,level) \
   if(level >= gVFNDebugMessageLevel) cout << message;
   
#define VFN_DEBUG_ENTRY(message,level) \
   if(level >= gVFNDebugMessageLevel) \
   {\
      for(unsigned int iii=0;iii<gVFNDebugMessageIndent;iii++) cout <<"  ";\
      cout << "%DEBUG:"<< level << " <"\
      << message << " (at " << __FILE__ << "," << __LINE__ << ")" <<endl;\
      gVFNDebugMessageIndent++;\
   }

#define VFN_DEBUG_EXIT(message,level) \
    if(level >= gVFNDebugMessageLevel) \
   {\
      if(gVFNDebugMessageIndent>0) gVFNDebugMessageIndent--;\
      for(unsigned int iii=0;iii<gVFNDebugMessageIndent;iii++) cout <<"  ";\
      cout << "%DEBUG:"<< level << " \\"\
      << message << "> (at " << __FILE__ << "," << __LINE__ << ")" <<endl;\
   }

    
#define VFN_DEBUG_GLOBAL_LEVEL(level) gVFNDebugMessageGlobalLevel=level;\
   gVFNDebugMessageLevel=gVFNDebugMessageGlobalLevel;
   
#define VFN_DEBUG_LOCAL_LEVEL(level) if(level != -1) gVFNDebugMessageLevel=level; else gVFNDebugMessageLevel=gVFNDebugMessageGlobalLevel;

#else //__DEBUG__

#define VFN_DEBUG_MESSAGE(message,level)
#define VFN_DEBUG_MESSAGE_SHORT(message,level)
#define VFN_DEBUG_ENTRY(message,level)
#define VFN_DEBUG_EXIT(message,level)
#define VFN_DEBUG_GLOBAL_LEVEL(level)
#define VFN_DEBUG_LOCAL_LEVEL(level)


#endif //__DEBUG__

#endif // __VFN_DEBUG__H__
