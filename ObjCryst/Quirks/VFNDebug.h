/*
* Part of LibCryst++ : a Crystallographic computing library in C++
*
*  (c) 2000 Vincent FAVRE-NICOLIN
*           Laboratoire de Cristallographie
*           24, quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
*  Contact: Vincent.Favre-Nicolin@cryst.unige.ch
*           Radovan.Cerny@cryst.unige.ch
*
*/
/*   VFNDebug.h
*  header file for debugging objects & functions
*
*/
#ifndef __VFN_DEBUG__H__
#define __VFN_DEBUG__H__

#ifdef __DEBUG__

///Set the Global debug level for messages
void LibCrystDebugGlobalLevel(const int level);

///Use this for a local modification of debug level messages. Call this at the
///beginning of the function, and call it without argument at the ned of the function
///to revert to the default global debug level.
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

#define VFN_DEBUG_MESSAGE(level,message)
#define VFN_DEBUG_MESSAGE_SHORT(message,level)
#define VFN_DEBUG_ENTRY(level,message)
#define VFN_DEBUG_EXIT(level,message)
#define VFN_DEBUG_GLOBAL_LEVEL(level)
#define VFN_DEBUG_LOCAL_LEVEL(level)


#endif //__DEBUG__

#endif // __VFN_DEBUG__H__
