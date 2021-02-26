#include <cstdio>

#include "sieve.hh"
#include "mersenne.hh"
#include "omp.h"

int main() {

//    sieve s(1000000);
	sieve s(200000);
    s.find_primes();
    double t0=omp_get_wtime();

	puts("M:");
	int count=0;
#pragma omp parallel reduction(+:count)
    {
        mersenne m(82589935/32+1);

#pragma omp for ordered schedule(dynamic)
        for(int i=0;i<s.tot;i++) {
            m.initialize(82589933);
            int j=s.prime[i];
            unsigned int r=m.remainder(j);
#pragma omp ordered
            if(r==0) {
				printf("%d\n",j);
				count++;
			}
        }
    }

    printf("# Time taken: %g s. Primes found: %d\n",omp_get_wtime()-t0,count);

	puts("N:");
	count=0;
#pragma omp parallel reduction(+:count)
    {
        mersenne m(82589935/32+1);

#pragma omp for ordered schedule(dynamic)
        for(int i=0;i<s.tot;i++) {
            m.initialize(82589932);
            int j=s.prime[i];
            unsigned int r=m.remainder(j);
#pragma omp ordered
            if(r==0) {
				printf("%d\n",j);
				count++;
			}
        }
    }

    printf("# Time taken: %g s. Primes found: %d\n",omp_get_wtime()-t0,count);
}
