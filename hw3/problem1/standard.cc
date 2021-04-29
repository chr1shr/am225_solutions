#include "standard.hh"

void standard(int n,double *a,double *b,double *c) {
    for(double *ce=c+n*n;c<ce;b+=n) for(double *aa=a,*cr=c+n;c<cr;c++,aa++) {
        *c=*aa*(*b);
        for(double *bb=b+1,*aaa=aa+n;aaa<a+n*n;aaa+=n,bb++) *c+=*aaa*(*bb);
    }
}
