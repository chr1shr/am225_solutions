#include <cstdio>
#include <cmath>

#include "omp.h"

#include "squircle.hh"

int main() {

#pragma omp parallel for ordered schedule(dynamic)
    for(int i=0;i<=20;i++) {

        // Create the finite-element problem and set the Neumann boundary
        // condition to match the manufactured solution
        int j=int(4*pow(100,(1/20.)*i)+0.5);
        double t0=omp_get_wtime(),t1,t2,nor;
        squircle sq(j,8);

        // Initialize the source term for the manufactured solution, solve, and
        // calculate the L2 error
        sq.set_source(1);
        t1=omp_get_wtime();
        sq.solve();
        t2=omp_get_wtime();
        nor=sq.l2_norm_mms(1);

        // Print the grid size, L2 error. Print the time taken for the problem
        // setup and the conjugate gradient solve.
#pragma omp ordered
        printf("%d %g %g %g %g\n",j,sq.h,nor,t1-t0,t2-t1);
    }
}
