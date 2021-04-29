#ifndef RBF_HILBERT_HH
#define RBF_HILBERT_HH

#include "rbf.hh"

struct hil_coord {
    int ind;
    unsigned int q;
    friend bool operator<(const hil_coord& a,const hil_coord& b) {
        return a.q<b.q;
    }
};

class rbf_hilbert : public rbf {
	public:
        rbf_hilbert(int n_,int type_) : rbf(n_,type_) {}
        void init_hilbert();
        void matrix_analysis(int bls,int &tot,int &pre);
        unsigned int pos(double x,double y);
};

#endif
