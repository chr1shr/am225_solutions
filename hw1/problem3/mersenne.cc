#include "mersenne.hh"
#include "mersenne.hh"

#include <cstdio>
#include <cstdlib>

void mersenne::initialize(int pow) {
    unsigned int* dp=d;

    if(pow>n*32) {
        fputs("Memory allocation not large enough\n",stderr);
        exit(1);
    }

    while(pow>=32) {
        *(dp++)=~0;pow-=32;
    }

    if(pow>0) *(dp++)=(1<<pow)-1;
    while(dp<d+n) *(dp++)=0;
}

unsigned int mersenne::remainder(unsigned int x) {

    unsigned long r=0;
    for(int k=n-1;k>=0;k--) {
        r<<=32;r+=static_cast<unsigned long>(d[k]);
        d[k]=static_cast<unsigned int>(r/x);
        r%=x;
    }
    return static_cast<unsigned int>(r);
}
