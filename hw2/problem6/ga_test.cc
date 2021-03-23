#include "galaxy.hh"
#include <cstdio>

int main() {


	// part b

    galaxy_sun ga(1.25,1,0.75,1,1,0.25);
    ga.solve_fixed(2000,40000,true);
	/*printf("\n\n");

    galaxy_rk4 ga_r(1.25,1,0.75,1,1,0.25);
    ga_r.solve_fixed(2000,40000,true);
	

	// part c
	galaxy_sun ga(1.25,1,0.75,1,1,0.25);
	ga.solve_fixed(1e5,20e5,true);*/
}
