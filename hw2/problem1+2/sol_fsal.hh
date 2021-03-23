#ifndef SOL_FSAL_HH
#define SOL_FSAL_HH

#include "sol_adapt.hh"

class sol_fsal : public sol_adapt {
    public:
        sol_fsal(int dof_) : sol_adapt(4,dof_) {}
        virtual ~sol_fsal() {}
        virtual double step_and_error(double dt,double atol,double rtol);
        virtual void dense_output(double theta,double dt);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    protected:
        void pre_accept() {}
};

#endif
