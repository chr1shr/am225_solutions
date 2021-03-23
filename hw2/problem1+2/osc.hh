#ifndef OSC_HH
#define OSC_HH

#include "sol_fsal.hh"
#include "sol_ckr.hh"

#include <cstdio>
#include <cmath>

/** This class has functions to specify the test oscillator problem. */
class osc {
    public:
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void osc_ff(double t_,double *in,double *out) {
            *out=-in[1]*t_;
            out[1]=*in*t_;
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void osc_init(double *q) {
            *q=1;
            q[1]=0;
        }
        inline void osc_print_step(double t,double *q) {
            printf("%g %g %g %g %g\n",t,*q,q[1],*q-sol0(t),q[1]-sol1(t));
        }
        inline void osc_print_dense(double t_,double *in) {
            printf("%g %g %g %g %g\n",t_,*in,in[1],*in-sol0(t_),in[1]-sol1(t_));
        }
        inline double sol0(double t) {return cos(0.5*t*t);}
        inline double sol1(double t) {return sin(0.5*t*t);}
};

/** Class to solve the oscillator problem with the Euler method. */
class osc_fsal : public sol_fsal, public osc {
    public:
        osc_fsal() : sol_fsal(2) {}
        virtual void ff(double t_,double *in,double *out) {
            osc_ff(t_,in,out);
        }
        virtual void init() {osc_init(q);}
        virtual void print_step() {osc_print_step(t,q);}
        virtual void print_dense(double t_,double *in) {osc_print_dense(t_,in);}
};

/** Class to solve the oscillator problem with the Ralston method. */
class osc_ckr : public sol_ckr, public osc {
    public:
        osc_ckr() : sol_ckr(2) {}
        virtual void ff(double t_,double *in,double *out) {
            osc_ff(t_,in,out);
        }
        virtual void init() {osc_init(q);}
        virtual void print_step() {osc_print_step(t,q);}
        virtual void print_dense(double t_,double *in) {osc_print_dense(t_,in);}
};



#endif
