#ifndef __VFN_CHRONOMETER__
#define __VFN_CHRONOMETER__


#include <stdlib.h>
#include <ctime>
#include <iostream.h>

class Chronometer;
class Chronometer
{
	public:
		Chronometer(){this->start();};
		~Chronometer(){};
		void start() {mPaused=false;time0=clock();};
		void pause() {time1=clock();mPaused=true;};
		void resume(){time0=clock()-(time1-time0);mPaused=false;};
		void print() 
		{
         if(mPaused == false) time1=clock();
			cout.setf(ios::fixed);
			int tmp=cout.precision(2);
			cout << "Elapsed time : " << (time1-time0)/(double)CLOCKS_PER_SEC << " s."<<endl ;
			cout.precision(tmp);
			cout.unsetf(ios::fixed);
		};
		double seconds()
      {
         if(mPaused ==false) time1=clock();
         return (time1-time0)/(double)CLOCKS_PER_SEC;
      };
	private:
		clock_t time0;
		clock_t time1;
      bool mPaused;
};

#endif
