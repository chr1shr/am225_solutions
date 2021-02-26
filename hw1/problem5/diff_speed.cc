#include "diffuse.hh"

#include <cstdio>

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif

int main() {
    double t0=wtime(),tbase;

    // Loop over a range of grid sizes
    for(int i=64;i<=8192;i<<=1) {
        diffuse di(i,1.25,-1.,1.);
        for(int j=1;j<=8;j++) {

            // Set threads and set up initial condition for diffusion solver
            di.num_t=j;
            di.init();

            // Measure the time for 1000 timesteps 
            t0=wtime();
            di.integrate(1000,false);
            t0=wtime()-t0;

            // Print timing info
            if(j==1) tbase=t0;
            printf("%d %d %g %g\n",i,j,t0,tbase/(t0*j));
        }
        puts("\n");
    }
}
