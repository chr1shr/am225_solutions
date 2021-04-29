#ifndef EXTRA_HH
#define EXTRA_HH

#include <cstdio>
#include <cstdlib>

inline double urand() {
    return -1+2./RAND_MAX*static_cast<double>(rand());
}

double* random_matrix(int nn);
void print_matrix(int n,double *a);
void mat_mul_blas(int n,double *a,double *b,double *c);

#endif
