#include <cstdio>

#include "custom_rng.hh"

#ifdef _OPENMP
#include "omp.h"
inline int max_t() {return omp_get_max_threads();}
inline int t_num() {return omp_get_thread_num();}
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline int max_t() {return 1;}
inline int t_num() {return 0;}
inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif

int main() {
    double t0=wtime();
    long d=0;
#pragma omp parallel
    {

        // Create random number generator within each thread. Each
        // is seeded differently.
        custom_rng cr(t_num());

#pragma omp for reduction(+:d)
        for(long i=0;i<1000000000L;i++) {

            // Since minimum number of tries is two, compute those first
            // without doing any checks
            int k=2;
            double q=cr.doub()+cr.doub();

            // Keep adding numbers until the sum exceeds 1
            while(q<1.) {
                q+=cr.doub();k++;
            }
            d+=k;
        }
    }
    printf("Time taken: %g s\n",wtime()-t0);
    printf("Expected winnings: $%.10g\n",double(d)*1e-7);

}
