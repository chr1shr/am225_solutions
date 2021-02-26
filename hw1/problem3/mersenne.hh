#ifndef MERSENNE_HH
#define MERSENNE_HH

class mersenne {
    public:
        int n;
        unsigned int* d;
        mersenne(int n_) : n(n_), d(new unsigned int[n]) {}
        ~mersenne() {
            delete [] d;
        }
        void initialize(int power);
        unsigned int remainder(unsigned int x);
};

#endif
