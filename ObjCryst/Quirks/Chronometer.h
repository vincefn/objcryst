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
		REAL seconds()
      {
         if(mPaused ==false)
			{
				mTime1=clock();
				mTimeSec1=time(0);
			}
			if((mTimeSec1-mTimeSec0)>100) return (mTimeSec1-mTimeSec0);
         return (mTime1-mTime0)/(REAL)CLOCKS_PER_SEC;
      }
	private:
		clock_t mTime0;
		clock_t mTime1;
		time_t mTimeSec0;
		time_t mTimeSec1;
      bool mPaused;
};

#endif
