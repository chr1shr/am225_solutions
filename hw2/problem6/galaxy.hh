#ifndef GALAXY_HH
#define GALAXY_HH

#include "sol_rk4.hh"
#include "sol_sun.hh"

/** This class has functions to specify the test Brusselator problem. */
class galaxy {
    public:
        const double a,b,c;
        const double aai,bbi,cci;
        const double A,C;
        const double Omega;
        galaxy(double a_,double b_,double c_,double A_,double C_,double Omega_) :
            a(a_), b(b_), c(c_), aai(1/(a*a)), bbi(1/(b*b)), cci(1/(c*c)), 
            A(A_), C(C_), Omega(Omega_) {}
        /** Evaluates the function f(x,y) on the RHS of the ODE.
         * \param[in] t the dependent x variable in Hairer et al.'s notation.
         * \param[in] in the array containing y variable.
         * \param[in] out the function f(x,y). */
        inline void galaxy_ff(double t_,double *in,double *out) {
            double fac=-2*A/(C+*in*(*in)*aai+in[1]*in[1]*bbi+in[2]*in[2]*cci);
            *out=in[3]+Omega*(in[1]);
            out[1]=in[4]-Omega*(*in);
            out[2]=in[5];
            out[3]=aai*fac*(*in)+Omega*in[4];
            out[4]=bbi*fac*in[1]-Omega*in[3];
            out[5]=cci*fac*in[2];
        }
        /** Sets up the initial conditions for the ODE.
         * \param[in] q the array to write to. */
        inline void galaxy_init(double *q) {
            *q=2.5;
            q[1]=0;
            q[2]=0;
            q[3]=0;
            q[4]=1.68883700590448;
            q[5]=0.2;
        }
};

/** Class to solve the Brusselator problem with the Euler method. */
class galaxy_sun : public sun, public galaxy {
    public:
        galaxy_sun(double a_,double b_,double c_,double A_,double C_,double Omega_)
            : sun(6), galaxy(a_,b_,c_,A_,C_,Omega_) {}
        virtual void ff(double t_,double *in,double *out) {
            galaxy_ff(t_,in,out);
        }
        virtual void init() {galaxy_init(q);}
};

/** Class to solve the Brusselator problem with the fourth-order Runge-Kutta
 * method. */
class galaxy_rk4 : public rk4, public galaxy {
    public:
        galaxy_rk4(double a_,double b_,double c_,double A_,double C_,double Omega_)
            : rk4(6), galaxy(a_,b_,c_,A_,C_,Omega_) {}
        virtual void ff(double t_,double *in,double *out) {
            galaxy_ff(t_,in,out);
        }
        virtual void init() {galaxy_init(q);}
};

#endif
