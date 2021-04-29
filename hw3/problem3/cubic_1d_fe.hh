#ifndef CUBIC_1D_FE_HH
#define CUBIC_1D_FE_HH

#include "conj_grad.hh"

#include <cstdio>
#include <cmath>

/** Class for solving an elliptic PDE problem over the domain [1,2] using
 * piecewise cubic basis functions. */
class cubic_1d_fe : public conj_grad {
    public:
        /** The number of intervals to divide the domain into. */
        const int n;
        /** The source function. */
        double* const f;
        /** The grid spacing. */
        const double h;
        /** The Neumann condition to apply at x=2. */
        double g;
        cubic_1d_fe(int n_) : conj_grad(2*n_+1),
           n(n_), f(new double[2*n+2]), h(1./n) {}
        virtual ~cubic_1d_fe() {delete [] f;}
        void init_const();
        void init_slope();
        void init_mms();
        double l2_norm_mms();
        void print_matrix();
        void print(FILE *fp);
        void print(const char* filename);
        inline void print() {print(stdout);}
        void print_fe_function(FILE *fp,int I);
        void print_fe_function(const char* filename,int I);
        inline void print_fe_function(int I) {print_fe_function(stdout,I);}
        virtual void mul_A(double *in,double *out);
    private:
        double sol(double xx);
        inline double mms_dsq(double xx,double s) {
            double del=exp(1-xx)*sin(5*M_PI*xx)-s;
            return del*del;
        }
        void assemble_b();
};

#endif
