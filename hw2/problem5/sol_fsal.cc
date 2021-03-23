#include "sol_fsal.hh"

#include <cmath>

/** Computes a Hermite interpolation of the solution, for dense output. The
 * result is stored into the k3 array.
 * \param[in] th the fraction of the timestep at which to evaluate the Hermite
 *               interpolant.
 * \param[in] dt the length of the current timestep. */
void sol_fsal::dense_output(double th,double dt) {
    double mth=1-th;

    // The function assumes that the current solution is in q, the new solution
    // is in dq, the current derivative is in k1, and the new derivative is in
    // k2
    for(int i=0;i<dof;i++)
        k3[i]=mth*q[i]+th*dq[i]
             -th*mth*((1-2*th)*(dq[i]-q[i])+dt*(th*k4[i]-mth*k1[i]));
}

/** Performs a Runge-Kutta step. It assumes k1 is already available.
 * \param[in] dt the step size.
 * \param[in] atol the absolute tolerance.
 * \param[in] rtol the relative tolerance.
 * \return A scaled estimate of error over the current integration step. */
double sol_fsal::step_and_error(double dt,double atol,double rtol) {

    // Second RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+(1/3.)*dt*k1[i];
    ff(t+(1/3.)*dt,dq,k2);

    // Third RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(-(1/3.)*k1[i]+k2[i]);
    ff(t+(2/3.)*dt,dq,k3);

    // Fourth RK step
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*(k1[i]-k2[i]+k3[i]);
    ff(t+dt,dq,k4);

    // Complete solution that is fourth-order accurate
    for(int i=0;i<dof;i++) dq[i]=q[i]+dt*0.125*(k1[i]+3*(k2[i]+k3[i])+k4[i]);

    // Perform FSAL step needed for error estimation, reusing k4 since it is no
    // longer needed
    ff(t+dt,dq,k4);fcount+=4;

    // Compute normalized error estimate
    double err=0.,qhat,o;
    for(int i=0;i<dof;i++) {

        // Compute third-order solution for error
        qhat=q[i]+dt*((1/12.)*k1[i]+0.5*k2[i]+0.25*k3[i]+(1/6.)*k4[i]);

        // Compute scaled error contribution
        o=(dq[i]-qhat)/(atol+rtol*max(fabs(dq[i]),fabs(qhat)));
        err+=o*o;
    }
    return sqrt(err/dof);
}
