/*  ObjCryst++ Object-Oriented Crystallographic Library
    (c) 2000-2007 Vincent Favre-Nicolin vincefn@users.sourceforge.net
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
#include <iostream>
#include "boost/date_time/posix_time/posix_time_types.hpp"
using namespace std;

/** Simple chronometer class, with microsecond precision
*
* Reported time correspond to \e real time, i.e. not the time the program
* has been running independently from other programs.
*/
class Chronometer
{
   public:
      Chronometer(){this->start();};
      ~Chronometer(){};
      void start() {mPaused=false;mTime0=boost::posix_time::microsec_clock::local_time();}
      void pause() {mTime1=boost::posix_time::microsec_clock::local_time();mPaused=true;}
      void resume()
      {
         mTime0=boost::posix_time::microsec_clock::local_time()-(mTime1-mTime0);
         mPaused=false;
      }
      void print() 
      {
         if(mPaused == false) mTime1=boost::posix_time::microsec_clock::local_time();
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
            mTime1=boost::posix_time::microsec_clock::local_time();
         }
         return (mTime1-mTime0).total_microseconds()/1.0e6;
      }
   private:
      bool mPaused;
      boost::posix_time::ptime mTime0;
      boost::posix_time::ptime mTime1;
};

#endif
