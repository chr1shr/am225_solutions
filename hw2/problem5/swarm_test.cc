#include "swarm.hh"
#include <cstdio>
#include "omp.h"

int main() {
	
	double J[] = {0.5,0.3,1.};
	double K[] = {0.5,-0.2,-0.2};
	char odir[256];

	for (int i=0;i<3;i++) {
		printf("%d\n",i);
		double t0=omp_get_wtime();
		sprintf(odir,"sw.%d",i);
		swarm_fsal sw(1250,J[i],K[i],odir);
		sw.solve(200,1e-6,1e-6,false,600);
		printf("T %g\n",omp_get_wtime()-t0);
	}
}
