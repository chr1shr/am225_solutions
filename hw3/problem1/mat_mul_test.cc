#include <cstdio>
#include <cstdlib>

#include "extra.hh"
#include "standard.hh"
#include "strassen.hh"

int main() {
    const int n=8,nn=n*n;
    double *a=random_matrix(nn),*b=random_matrix(nn),
           *c=new double[nn],*d=new double[nn];

    standard(n,a,b,c);
    strassen(n,a,b,d);
    print_matrix(n,a);puts("");
    print_matrix(n,b);puts("");
    print_matrix(n,c);puts("");
    print_matrix(n,d);puts("");

    delete [] a;
    delete [] b;
    delete [] c;
    delete [] d;
}
