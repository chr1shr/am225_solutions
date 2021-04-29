#ifndef POISSON_FFT
#define POISSON_FFT

#include <functional>
#include <fftw3.h>

class poisson_fft
{
    public:
        poisson_fft(int n_);
        ~poisson_fft();
        void solve();
        void solve(const std::function<double(double,double)>& fun);
        void print_solution();
        /** The number of interior gridpoints in one dimension */
        const int n;
        /** The total number of interior gridpoints */
        const int nn;
        /** The grid spacing */
        const double h;
        /** The inverse grid spacing squared */
        const double ih2;
        /** The source term / solution */
        double* const v;
    protected:
        /** Eigenvalues of the discrete second derivative */
        double* lam;
        /** DST plan */
        fftw_plan plan;
};

#endif
