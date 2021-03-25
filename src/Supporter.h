 #ifndef  _SUPPORTER_H_
#define  _SUPPORTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include<iostream> 
#include <string>

#ifndef _MSC_VER
#include <chrono>
#else
#include <time.h>
#endif

#define INFTY 1.e30

using namespace std;
class Timer{
public:
	static double CPU_time(){
		#ifndef _MSC_VER
		  using namespace std::chrono;
		  auto now = time_point_cast<milliseconds>(system_clock::now());
		  using sys_milliseconds = decltype(now);
		  return now.time_since_epoch().count()/1000.0;
		#else
			return (double)clock()/(double)CLOCKS_PER_SEC;
		#endif
	}
};

class Format{
public:
	static void pretty_print(string str){
		cout << " ------>> " << str << endl;
	};
};




#endif  // _SUPPORTER_H_