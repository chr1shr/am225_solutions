#ifndef SOL_CKR_HH
#define SOL_CKR_HH

#include "sol_adapt.hh"

class sol_ckr : public sol_adapt {
    public:
        sol_ckr(int dof_);
        virtual ~sol_ckr();
        virtual double step_and_error(double dt,double atol,double rtol);
        virtual void dense_output(double theta,double dt);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    protected:
        virtual void pre_accept();
    private:
        void cash_karp_step(double t_,double dt,double *in,double *out,double *k1_);
        double* const k5;
        double* const k1b;
        double* const q1;
        double* const q2;
};

#endif
