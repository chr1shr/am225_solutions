#include <cmath>
#include <cstring>

#include "lax_wend.hh"

int main() {

    // Period of the solution
    const double T=3*M_PI/sqrt(5);

    // Loop over a variety of different grid sizes
#pragma omp parallel for schedule(dynamic) ordered
    for(int i=0;i<=15;i++) {

        // Create the finite-element problem and set the Neumann boundary
        // condition to match the manufactured solution
        int k,j=int(1024*pow(32,(1./15.)*i)+0.5);
        lax_wend lw(j);
        double *z=new double[j],l2=0.,d;

        // Set up the initial condition and store a copy of it
        lw.init_triangle();
        memcpy(z,lw.a,sizeof(double)*j);

        // Integrate the PDE for one complete period and compute L2 norm
        // between the initial condition and solution
        lw.solve(T,1/3.);
        for(k=0;k<j;k++) {d=z[k]-lw.a[k];l2+=d*d;}
        l2=sqrt(l2*lw.dx);

        // Print the L2 norm (ordered according to grid size and free the
        // dynamically allocated memory
#pragma omp ordered
        printf("%d %g %g\n",j,lw.dx,l2);
        delete [] z;
    }
}
