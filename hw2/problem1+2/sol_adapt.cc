#include "sol_adapt.hh"

#include <cstdio>
#include <cmath>
#include <cstring>

/** The class constructor when the number of degrees of freedom are known.
 * \param[in] dof_ the number of degrees of freedom. */
sol_adapt::sol_adapt(int order_,int dof_) : order(order_), dof(dof_),
    fcount(0), t(0.), q(new double[dof]),
    dq(new double[dof]), k1(new double[dof]), k2(new double[dof]),
    k3(new double[dof]), k4(new double[dof]) {}

/** The class destructor frees the dynamically allocated memory. */
sol_adapt::~sol_adapt() {
    delete [] k4;
    delete [] k3;
    delete [] k2;
    delete [] k1;
    delete [] dq;
    delete [] q;
}

/** Prints the current state of the solution. */
void sol_adapt::print_step() {
    printf("%g",t);
    for(int i=0;i<dof;i++) printf(" %g",q[i]);
    puts("");
}

/** Prints the current state of the solution. */
void sol_adapt::print_dense(double t_,double *in) {
    printf("%g",t_);
    for(int i=0;i<dof;i++) printf(" %g",in[i]);
    puts("");
}

/** Performs a time integration of the ODE problem, using adaptive step size
 * selection.
 * \param[in] t_end the end point of the integration.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance.
 * \param[in] output whether to output the integration steps.
 * \param[in] d_steps the number of dense output intervals, set to zero if
 *                    dense output is not required. */
void sol_adapt::solve(double t_end,double atol,double rtol,bool output,int d_steps) {

    // Set up initial condition and compute timestep
    init();
    double dt=initial_step_size(atol,rtol),err,t_den=0,dt_den;
    bool last=false;
    if(t+dt>t_end) {last=true;dt=t_end-t;}

    // Do any required output of the initial step
    if(output) print_step();
    if(d_steps>0) {
        dt_den=t_end/d_steps;
        print_dense(t,q);
    }

    while(true) {

        // Perform integration step and check if the scaled error term is
        // acceptable
        if((err=step_and_error(dt,atol,rtol))<1.0) {
            t+=dt;

            // Some methods (e.g. Cash--Karp Richardson) may need to do extra
            // evaluations here, to prepare for the next step
            pre_accept();

            // Do any dense output interpolation
            if(d_steps>0) {
                while(t_den+dt_den<t) {
                    t_den+=dt_den;
                    dense_output(1.+(t_den-t)/dt,dt);
                    print_dense(t_den,k3);
                }
            }

            // Copy the solution and first RK step into the correct arrays
            memcpy(k1,k4,dof*sizeof(double));
            memcpy(q,dq,dof*sizeof(double));

            // Print solution and check for the termination condition
            if(output) print_step();
            if(last) return;
        }

        // Compute the new timestep. If it exceeds the end of the integration
        // interval, then truncate the timestep and mark this as the last step.
        dt*=min(facmax,max(facmin,safe_fac*pow(err,-1./order)));
        if(t+dt>t_end) {last=true;dt=t_end-t;} else last=false;
    }
}

/** Estimates an initial step size choice.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance. */
double sol_adapt::initial_step_size(double atol,double rtol) {

    // Compute d_0 and d_1
    ff(t,q,k1);
    double d0=scaled_norm(q,atol,rtol),
           d1=scaled_norm(k1,atol,rtol),d2=0,h0,h1,md12,o;

    // Make initial step size guess
    h0=fabs(d0)<1e-5||fabs(d1)<1e-5?1e-6:0.01*(d0/d1);

    // Perform one explicit step with initial step size
    for(int i=0;i<dof;i++) k2[i]=q[i]+h0*k1[i];
    ff(t+h0,k2,k3);fcount+=2;

    // Estimate second derivative of the solution
    for(int i=0;i<dof;i++) {
        o=(k1[i]-k3[i])/(atol+rtol*max(fabs(q[i]),fabs(k2[i])));
        d2+=o*o;
    }
    d2=sqrt(d2/dof)/h0;

    // Compute second step size estimate
    md12=max(d1,d2);
    h1=md12<=1e-15?max(1e-6,h0*1e-3)
                  :pow(0.01/md12,1/(order+1));

    // Choose the minimum of the two possibilities
    return min(100*h0,h1);
}

/** Computes the scaled norm used in initial step size selection.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance. */
double sol_adapt::scaled_norm(double *in,double atol,double rtol) {
    double no=0.,o;
    for(int i=0;i<dof;i++) {
        o=in[i]/(atol+rtol*fabs(q[i]));
        no+=o*o;
    }
    return sqrt(no/dof);
}
