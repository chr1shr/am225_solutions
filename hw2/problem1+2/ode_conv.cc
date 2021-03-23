#include "sol_fsal.hh"
#include "brusselator.hh"

#include <cstdio>
#include <cmath>

#include "omp.h"

int main() {

    // Compute reference solution
//    brus_ckr br;
	brus_fsal br;
    br.solve(20.,3e-15,3e-15);
    double ref0=br.q[0],ref1=br.q[1];

    // Allocate space for storing accuracy
	const int N=105;
    int fcount[N+1];
    double err[N+1],t0=omp_get_wtime();

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for schedule(dynamic)
    for(int i=0;i<=N;i++) {

        // Compute the number of steps according to a power law. Minimum
        // steps of 100, and maximum steps of 100000.
//        brus_ckr *sb=new brus_ckr();
        brus_fsal *sb=new brus_fsal();

        // Perform integration and compute the difference to the reference
        // solution
        double tol=0.001*pow(0.1,0.1*i);
        sb->solve(20.,tol,tol);
        double dy0=ref0-sb->q[0],
               dy1=ref1-sb->q[1];
        err[i]=sqrt(dy0*dy0+dy1*dy1);
        fcount[i]=sb->fcount;
        delete sb;
    }

    // Output the results in a Gnuplot-readable format
    printf("# Time taken : %g ms\n",1e3*(omp_get_wtime()-t0));
    for(int i=0;i<=N;i++) printf("%d %g\n",fcount[i],err[i]);
}
