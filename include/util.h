#ifndef __IBM_MERLIN_UTIL__
#define __IBM_MERLIN_UTIL__

#include<cmath>
#include<ctime>
#include<sys/time.h>
#include<cassert>
#include<cstdlib>
#include<stdint.h>
#include<limits>

#include<string>
#include<sstream>
#include<vector>

namespace merlin {


/*
#ifdef WINDOWS
    #include <windows.h>
//    #include <boost/math/special_functions/atanh.hpp>  // for atanh
//    #include <boost/math/special_functions/log1p.hpp>  // for log1p
    #include <float.h>  // for _isnan
#else
    // Assume POSIX compliant system. We need the following for querying the system time
    #include <sys/time.h>
#endif
*/


/* // from libDAI
#ifdef WINDOWS
double atanh( double x ) {
    return boost::math::atanh( x );
}
double log1p( double x ) {
    return boost::math::log1p( x );
}
#endif
*/


//static inline bool isfinite(value v) { return std::abs(v)!=std::numeric_limits<value>::infinity(); }
inline bool isfinite(double v) { return (v <= std::numeric_limits<double>::max() && v >= -std::numeric_limits<double>::max()); }
//inline bool isnan(double v)    { return (v!=v); }
inline double infty()          { return std::numeric_limits<double>::infinity(); }


// Returns system (wall clock) time in seconds
inline double timeSystem() {
#ifdef WINDOWS
    SYSTEMTIME tbuf;
    GetSystemTime(&tbuf);
    return( (double)(tbuf.wSecond + (double)tbuf.wMilliseconds / 1000.0) );
#else
    struct timeval tv;
    gettimeofday( &tv, NULL );
    return( (double)(tv.tv_sec + (double)tv.tv_usec / 1000000.0) );
#endif
}

inline std::string timestamp() {
	struct tm * dt;
	char buffer[30];
	time_t t = time(NULL);
	dt = localtime(&t);
	strftime(buffer, sizeof(buffer), "[%F %X]", dt);
	return std::string(buffer);
}

inline double timeProcess() {
#ifdef WINDOWS
    throw std::runtime_error("No process time implemented for Windows");
#else
    clock_t tv( clock() );
    return( (double)tv );
#endif
}


//
// String splitting functions
//
inline std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}




//
// Random number classes
//
inline void rand_seed() { srand(time(0)); }
inline void rand_seed(size_t s) { srand(s); }
inline double randu()  {
	return rand() / double(RAND_MAX); 
}
// randi returns a random integer in 0..imax-1
inline int randi(int imax) { 
 assert(imax>0);
 imax--;
 if (imax==0) return 0;
 int guard = (int) (randu()*imax)+1;
 return (guard > imax)? imax : guard;
}
inline double randn() {  // Marsaglia polar method
  double u,v,s; 
	u=2*randu()-1; v=2*randu()-1; s=u*u+v*v;
	return u*std::sqrt(-2*std::log(s)/s);
} 


} // namespace

#endif // re-include

