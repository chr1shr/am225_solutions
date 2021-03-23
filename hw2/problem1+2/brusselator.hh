#ifndef BRUSSELATOR_HH
#define BRUSSELATOR_HH

#include "sol_fsal.hh"
#include "sol_ckr.hh"

/** This class has functions to specify the test Brusselator problem. */
class brusselator {
    public:
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void brus_ff(double t_,double *in,double *out) {
            double &y1=*in,&y2=in[1];
            *out=1+y1*(y1*y2-4);
            out[1]=y1*(3-y1*y2);
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void brus_init(double *q) {
            *q=1.5;
            q[1]=3.;
        }
};

/** Class to solve the Brusselator problem with the Euler method. */
class brus_fsal : public sol_fsal, public brusselator {
    public:
        brus_fsal() : sol_fsal(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

/** Class to solve the Brusselator problem with the Ralston method. */
class brus_ckr : public sol_ckr, public brusselator {
    public:
        brus_ckr() : sol_ckr(2) {}
        virtual void ff(double t_,double *in,double *out) {
            brus_ff(t_,in,out);
        }
        virtual void init() {brus_init(q);}
};

#endif
