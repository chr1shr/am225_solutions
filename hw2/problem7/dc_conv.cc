#include "discont.hh"

#include <cstdio>
#include <cmath>

#include "omp.h"

int main() {
    const double min_steps=1e3,max_steps=1e7,
                 ratio_steps=max_steps/min_steps;

    const double min_atol=1e-2,
                 max_atol=1e-12,
                 ratio_atol=max_atol/min_atol;

    const double t_end=48+exp(-1);
    const bool step=true;

    // Allocate space for storing accuracy
    int fcount[303];
    double err[303],t0=omp_get_wtime();

	char fnames[3][256] = {"rk4.conv","fsal.conv","ckr.conv"};

    // Loop over a variety of grid sizes (given by i index) and integration
    // schemes (given by j index)
#pragma omp parallel for collapse(2) schedule(dynamic)
    for(int i=0;i<=100;i++) {
        for(int j=0;j<3;j++) {
            int fc;double xc,yc,fr=0.01*i;

            // Dynamically allocate a solver
            if(j==0) {
                discont_rk4 drk4(step);
                drk4.solve_fixed(t_end,int(min_steps*pow(ratio_steps,fr)),false);
                xc=drk4.q[0];yc=drk4.q[1];fc=drk4.fcount;
            } else {
                sol_adapt* sa=j==1?(sol_adapt*) new discont_fsal(step):
                                   (sol_adapt*) new discont_ckr(step);
                sa->solve(t_end,min_atol*pow(ratio_atol,fr),0.,false);
                xc=sa->q[0];yc=sa->q[1];fc=sa->fcount;
                delete sa;
            }

            // Perform integration and compute the difference to the reference
            double dx=xc-1,dy=yc-exp(-1);
            err[i+101*j]=sqrt(dx*dx+dy*dy);
            fcount[i+101*j]=fc;
        }
    }

    // Output the results in a Gnuplot-readable format
    printf("# Time taken : %g ms\n",1e3*(omp_get_wtime()-t0));
    for(int j=0;j<3;j++) {
		FILE *fh = fopen(fnames[j],"w");
        for(int i=0;i<=100;i++) fprintf(fh,"%d %g\n",fcount[i+101*j],err[i+101*j]);
		fclose(fh);
    }
}
