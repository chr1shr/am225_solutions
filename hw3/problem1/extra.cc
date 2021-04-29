#include "extra.hh"

extern "C" {
    void dgemm_(const char *transa, const char *transb,
                const int *m, const int *n, const int *k,
                const double *alpha, const double *a,
                const int *lda, const double *b, const int *ldb,
                const double *beta, double *c, const int *ldc);
}

double* random_matrix(int nn) {
    double* a=new double[nn];
    for(double *ap=a;ap<a+nn;ap++) *ap=urand();
    return a;
}

void print_matrix(int n,double *a) {
    for(double *ae=a+n;a<ae;a++) {
        printf("%7.4f",*a);
        for(double *aa=a+n;aa<a+n*n;aa+=n) printf(" %7.4f",*aa);
        fputc('\n',stdout);
    }
}

void mat_mul_blas(int n,double *a,double *b,double *c) {
    char trans='n';
    double alpha=1,beta=0;
    dgemm_(&trans,&trans,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}
