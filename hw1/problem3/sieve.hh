#ifndef SIEVE_HH
#define SIEVE_HH

#include <cstdlib>

class sieve {
    public:
        /** The maximum number to search for. */
        const int n;
        /** The list of primes. */
        int* prime;
        /** The total number of primes. */
        int tot;
        sieve(int n_) : n(n_), prime(NULL) {}
        ~sieve() {
            if(prime!=NULL) delete [] prime;
        }
        void find_primes();
};

#endif
