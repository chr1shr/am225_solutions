#include "sol_ckr.hh"

#include <cmath>

sol_ckr::sol_ckr(int dof_) : sol_adapt(6,dof_), k5(new double[dof]),
        k1b(new double[dof]), q1(new double[dof]), q2(new double[dof]) {}

sol_ckr::~sol_ckr() {
    delete [] q2;
    delete [] q1;
    delete [] k1b;
    delete [] k5;
}

/** Performs a fifth-order Cash--Karp step.
 * \param[in] t the current time.
 * \param[in] dt the timestep.
 * \param[in] in the initial solution.
 * \param[in] out the computed solution at the end of the step.
 * \param[in] k1_ a place to look for the first RK step (which is assumed to be
 *                already computed). */
void sol_ckr::cash_karp_step(double t_,double dt,double *in,double *out,double *k1_) {

    // Second CK step
    for(int i=0;i<dof;i++) dq[i]=in[i]+(1/5.)*dt*k1_[i];
    ff(t_+(1/5.)*dt,dq,k2);

    // Third CK step
    for(int i=0;i<dof;i++) dq[i]=in[i]+dt*((3/40.)*k1_[i]+(9/40.)*k2[i]);
    ff(t_+(3/10.)*dt,dq,k3);

    // Fourth RK step
    for(int i=0;i<dof;i++) dq[i]=in[i]+dt*((3/10.)*k1_[i]+(-9/10.)*k2[i]+(6/5.)*k3[i]);
    ff(t_+(3/5.)*dt,dq,k4);

    // Fifth RK step
    for(int i=0;i<dof;i++) dq[i]=in[i]+dt*((-11/54.)*k1_[i]+(5/2.)*k2[i]+(-70/27.)*k3[i]+(35/27.)*k4[i]);
    ff(t_+dt,dq,k5);

    // Sixth RK step (overwriting k5 since it won't be needed any more)
    for(int i=0;i<dof;i++) dq[i]=in[i]+dt*((1631/55296.)*k1_[i]+(175/512.)*k2[i]+(575/13824.)*k3[i]
                                          +(44275/110592.)*k4[i]+(253/4096.)*k5[i]);
    ff(t_+(7/8.)*dt,dq,k5);

    // Complete solution that is fifth-order accurate
    for(int i=0;i<dof;i++) out[i]=in[i]+dt*((37/378.)*k1_[i]+(250/621.)*k3[i]
                                           +(125/594.)*k4[i]+(512/1771.)*k5[i]);
    fcount+=5;
}

/** Performs a Runge-Kutta step. It assumes k1 is already available.
 * \param[in] dt the step size.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance.
 * \return A scaled estimate of error over the current integration step. */
double sol_ckr::step_and_error(double dt,double atol,double rtol) {
    const double rfac=1/31.;
    double hdt=0.5*dt;
    double err=0.,o;

    // Take two half-steps
    cash_karp_step(t,hdt,q,q1,k1);
    ff(t+hdt,q1,k1b);
    cash_karp_step(t+hdt,hdt,q1,q2,k1b);

    // Take a big step
    cash_karp_step(t,dt,q,dq,k1);
    fcount+=16;

    for(int i=0;i<dof;i++) {

        // Compute error estimate
        o=(q2[i]-dq[i])*rfac;

        // Correct q1 and dq to be sixth-order
        q1[i]+=0.5*o;
        dq[i]=q2[i]+o;

        // Store scaled error contribution
        o/=atol+rtol*max(fabs(dq[i]),fabs(q[i]));

        // Compute scaled error contribution
        err+=o*o;
    }
    err=sqrt(err/dof);

    return std::isnan(err)?100:sqrt(err/dof);
}

/** Computes a quintic interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] th the fraction of the timestep at which to evaluate the quintic
 *               interpolant.
 * \param[in] dt the length of the current timestep. */
void sol_ckr::dense_output(double th,double dt) {
    double al=2*th-1,be=th-1;
    for(int i=0;i<dof;i++) {
        k3[i]=th*th*(be*be)*(16*q1[i]+dt*8*al*k1b[i])
             +al*al*(be*be*(dt*th*k1[i]+(1+6*th)*q[i])
                    +th*th*(dt*be*k4[i]+(1-6*be)*dq[i]));
    }
}

void sol_ckr::pre_accept() {
    ff(t,dq,k4);fcount++;
}
