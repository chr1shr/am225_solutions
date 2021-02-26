#include "sieve.hh"

#include <cstdio>
#include <cmath>

void sieve::find_primes() {
    bool *a=new bool[n];

    for(bool *ap=a+2;ap<a+n;) *(ap++)=true;
    int sn=static_cast<int>(sqrt(n));

#pragma omp parallel for schedule(dynamic,8)
    for(int i=2;i<=sn;i++) {
        for(bool *ap=a+i*i;ap<a+n;ap+=i) {
#pragma omp atomic write
            *ap=false;
        }
    }

    tot=0;
    for(bool *ap=a+2;ap<a+n;ap++) if(*ap) tot++;

    printf("%d\n",tot);
    if(prime!=NULL) delete [] prime;
    prime=new int[tot];

    int *ptr=prime;
    for(int j=2;j<n;j++) if(a[j]) *(ptr++)=j;

    delete [] a;
}
