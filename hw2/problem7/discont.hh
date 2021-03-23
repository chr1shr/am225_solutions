#ifndef DISCONT_HH
#define DISCONT_HH

#include <cstdio>
#include <cmath>

#include "sol_rk4.hh"
#include "sol_fsal.hh"
#include "sol_ckr.hh"

class discont {
    public:
        bool step;
        discont(bool step_) : step(step_) {}
        inline void discont_ff(double t_,double *in,double *out) {
            if(fabs(*in)>fabs(in[1])) {
                *out=0;
                out[1]=f(*in);
            } else {
                *out=-f(in[1]);
                out[1]=0;
            }
        }
        inline double f(double x) {
            return step?(x>0?1:-1):x;
        }
        inline void discont_init(double *q) {*q=1;q[1]=0;}
};

class discont_fsal : public sol_fsal, public discont {
    public:
        discont_fsal(bool step_) : sol_fsal(2), discont(step_) {}
        virtual void ff(double t_,double *in,double *out) {discont_ff(t_,in,out);}
        virtual void init() {discont_init(q);}
        virtual void print_dense(double t_,double *in) {}
};

class discont_ckr : public sol_ckr, public discont {
    public:
        discont_ckr(bool step_) : sol_ckr(2), discont(step_) {}
        virtual void ff(double t_,double *in,double *out) {discont_ff(t_,in,out);}
        virtual void init() {discont_init(q);}
        virtual void print_dense(double t_,double *in) {}
};

class discont_rk4 : public rk4, public discont {
    public:
        discont_rk4(bool step_) : rk4(2), discont(step_) {}
        virtual void ff(double t_,double *in,double *out) {discont_ff(t_,in,out);}
        virtual void init() {discont_init(q);}
};

#endif
