#include "extra.hh"
#include "standard.hh"
#include "strassen.hh"

#include "omp.h"

// Amount of time (in seconds) to run each test for
const double duration=1.;

int main() {

    // Open output file
    FILE *f=fopen("mm_time.out","w");
    if(f==NULL) {
        fputs("Can't open output file\n",stderr);
        return 1;
    }

    // Loop over a variety of different array sizes
    double t0,t1,t2,t3;
    int n_strassen,n_standard,n_blas;
    for(int n=4;n<=8192;n<<=1) {
        n_strassen=n_standard=n_blas=0;

        // Allocate arrays
        int nn=n*n;
        double *a=random_matrix(nn),*b=random_matrix(nn),
               *c=new double[nn];

        // Time Strassen's algorithm
        t0=omp_get_wtime();
        do {
            strassen(n,a,b,c);
            t1=omp_get_wtime();n_strassen++;
        } while(t1<t0+duration);

        // Time standard O(n^3) algorithm
        do {
            standard(n,a,b,c);
            t2=omp_get_wtime();n_standard++;
        } while(t2<t1+duration);

        // Time BLAS routine
        do {
            mat_mul_blas(n,a,b,c);
            t3=omp_get_wtime();n_blas++;
        } while(t3<t2+duration);

        // Print timing info, dividing test durations by the number of trials
        t0=(t1-t0)/n_strassen;
        t1=(t2-t1)/n_standard;
        t2=(t3-t2)/n_blas;
        printf("Dim %d: Strassen %g s, standard %g s, BLAS %g s\n",n,t0,t1,t2);
        fprintf(f,"%d %g %g %g %d %d %d\n",n,t0,t1,t2,n_strassen,n_standard,n_blas);

        // Free dynamically allocated memory
        delete [] a;
        delete [] b;
        delete [] c;
    }
    fclose(f);
}
