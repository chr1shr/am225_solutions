#ifndef SOL_ADAPT_HH
#define SOL_ADAPT_HH

/** The minimum amount that the timestep can be reduced by in one step. */
const double facmin=1/3.;

/** The maximum amount that the timestep can be enlarged by in one step. */
const double facmax=3.;

/** A safety factor to scale the optimal timestep choice by, so that the next
 * step will be accepted with high probability. */
const double safe_fac=0.9;

class sol_adapt {
    public:
        /** The order of the method. */
        const int order;
        /** The total number of degrees of freedom in the ODE system. */
        int dof;
        /** A counter for the number of function evaluations. */
        int fcount;
        /** The current time. */
        double t;
        /** The solution vector. */
        double* const q;
        sol_adapt(int order_,int dof_);
        virtual ~sol_adapt();
        void solve(double t_end,double atol,double rtol,bool output=false,int d_steps=0);
        double initial_step_size(double atol,double rtol);
        virtual double step_and_error(double dt,double atol,double rtol) = 0;
        virtual void dense_output(double theta,double dt) = 0;
        virtual void print_step();
        virtual void print_dense(double t_,double *in);
        virtual void init() = 0;
        virtual void ff(double t_,double *in,double *out) = 0;
    protected:
        virtual void pre_accept() {}
        inline double min(double a,double b) {return a<b?a:b;}
        inline double max(double a,double b) {return a>b?a:b;}
        double* const dq;
        double* const k1;
        double* const k2;
        double* const k3;
        double* const k4;
    private:
        double scaled_norm(double *in,double atol,double rtol);
};

#endif
