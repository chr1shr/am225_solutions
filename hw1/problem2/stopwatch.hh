#ifndef STOPWATCH_HH
#define STOPWATCH_HH

#include "omp.h"

struct stopwatch {

	double t0;
	double el;

	stopwatch() : t0(0),el(0) {}

	inline void tic() {t0 = omp_get_wtime();}
	inline double toc() {
		double amt = omp_get_wtime()-t0; 
		el += amt;
		return amt;
	}
	inline double elapsed() {return el;}
	inline void reset() {el=0;}

};

#endif
