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
#ifndef __VFN_CHRONOMETER__
#define __VFN_CHRONOMETER__


#include <stdlib.h>
#include <time.h>
#include <iostream.h>

class Chronometer
{
   public:
      Chronometer(){this->start();};
      ~Chronometer(){};
      void start() {mPaused=false;mTime0=clock();mTimeSec0=time(0);}
      void pause() {mTime1=clock();mTimeSec1=time(0);mPaused=true;}
      void resume()
      {
         mTime0=clock()-(mTime1-mTime0);
         mTimeSec0=time(0)-(mTimeSec1-mTimeSec0);
         mPaused=false;
      }
      void print() 
      {
         if(mPaused == false) mTime1=clock();
         cout.setf(ios::fixed);
         int tmp=cout.precision(2);
         cout << "Elapsed time : " << this->seconds() << " s."<<endl ;
         cout.precision(tmp);
         cout.unsetf(ios::fixed);
      }
      float seconds()
      {
         if(mPaused ==false)
         {
            mTime1=clock();
            mTimeSec1=time(0);
         }
         if((mTimeSec1-mTimeSec0)>100) return (mTimeSec1-mTimeSec0);
         return (mTime1-mTime0)/(float)CLOCKS_PER_SEC;
      }
   private:
      clock_t mTime0;
      clock_t mTime1;
      time_t mTimeSec0;
      time_t mTimeSec1;
      bool mPaused;
};

#endif
