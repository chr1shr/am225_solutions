#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <cmath>
#include <complex>
#include <fftw3.h>

typedef std::complex<double> dcomplex;

class Helmholtz
{
    public:
        Helmholtz(int n_, double f_, dcomplex* er_ = nullptr);
        ~Helmholtz();
        void set_frequency(double f_);
        virtual void solve() = 0;
        void print_solution();
        /** The number of subdivisions in one dimension */
        const int n;
        /** The total number of unknowns in the system */
        const int nn;
        /** The inverse grid spacing squared */
        const double ih2;
        /** The free space magnetic permeability */
        const double mu0 = 4e-7*M_PI;
        /** The free space electric permittivity */
        const double e0 = 8.8542e-12;
        /** The source term */
        dcomplex* const v;
    protected:
        /** The frequency */
        double f;
        /** The angular frequency */
        double omega;
        /** The wavenumber squared */
        double k02;
        /** The spatially-varying relative electric permittivity */
        dcomplex* er;
};

class SlowHelmholtz : public Helmholtz
{
    public:
        SlowHelmholtz(int n_, double f_, dcomplex* er_ = nullptr);
        virtual void solve();
};

class FastHelmholtz : public Helmholtz
{
    public:
        FastHelmholtz(int n_, double f_, dcomplex* er_ = nullptr);
        ~FastHelmholtz();
        virtual void solve();
        /** True solution to compare against */
        const dcomplex* vtrue;
    private:
        void dst1_2d(dcomplex* x);
        void direct_solve();
        void iter_solve(int niter = 20);
        /** Eigenvalues of the discrete second derivative */
        double* lam;
        /** Storage for DST */
        double* temp;
        /** DST plan */
        fftw_plan plan;
};

#endif
