#ifndef SWARM_HH
#define SWARM_HH

#include <cstdlib>

#include "sol_fsal.hh"
#include "sol_ckr.hh"

class swarm {
    public:
        /** The number of actors. */
        const int N;
        const double J;
        const double K;
        const char* filename;
        swarm(int N_,double J_,double K_,const char* filename_);
        void swarm_ff(double t_,double *in,double *out);
        void swarm_init(double *q);
        void swarm_print_dense(double t_,double *in);
        inline double rnd(double l,double u) {
            return l+(u-l)/RAND_MAX*static_cast<double>(rand());
        }
    private:
        int d_o_frame;
};

class swarm_fsal : public sol_fsal, public swarm {
    public:
        swarm_fsal(int N_,double J_,double K_,const char* filename_) :
            sol_fsal(3*N_), swarm(N_,J_,K_,filename_) {}
        virtual void ff(double t_,double *in,double *out) {swarm_ff(t_,in,out);}
        virtual void init() {swarm_init(q);}
        virtual void print_dense(double t_,double *in) {swarm_print_dense(t_,in);}
};

class swarm_ckr : public sol_ckr, public swarm {
    public:
        swarm_ckr(int N_,double J_,double K_,const char* filename_) :
            sol_ckr(3*N_), swarm(N_,J_,K_,filename_) {}
        virtual void ff(double t_,double *in,double *out) {swarm_ff(t_,in,out);}
        virtual void init() {swarm_init(q);}
        virtual void print_dense(double t_,double *in) {swarm_print_dense(t_,in);}
};

#endif
