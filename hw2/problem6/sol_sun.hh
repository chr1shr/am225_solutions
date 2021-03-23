#ifndef SOL_SUN_HH
#define SOL_SUN_HH

#include "sol_base.hh"

/** Class for solving an ODE IVP using Geng Sun's symplectic fifth-order
 * method. */
class sun : public sol_base {
    public:
        sun(int dof_);
        virtual ~sun();
        virtual void step(double dt);
    private:
        double *dq;
        double *k1;
        double *k2;
        double *k3;
        double *k1b;
        double *k2b;
        double *k3b;
};

#endif
