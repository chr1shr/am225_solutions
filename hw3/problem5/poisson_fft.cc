#include <cmath>
#include "poisson_fft.hh"

poisson_fft::poisson_fft(int n_) :
    n(n_-1), nn(n*n), h(1./n_), ih2(n_*n_),
    v(fftw_alloc_real(nn)), lam(new double[n]),
    plan(fftw_plan_r2r_2d(n,n,v,v,FFTW_RODFT00,FFTW_RODFT00,FFTW_MEASURE))
{
    for (int i=0; i<n; ++i) {
        double s = std::sin((i+1)*M_PI/2.*h);
        lam[i] = 4*ih2*s*s;
    }
}

poisson_fft::~poisson_fft()
{
    fftw_destroy_plan(plan);
    delete [] lam;
    fftw_free(v);
}

/** Solve Poisson's equation using the DST. */
void poisson_fft::solve(const std::function<double(double,double)>& fun)
{
    // Evaluate the source term at the gridpoints
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            v[n*i+j] = fun((j+1)*h,(i+1)*h);
        }
    }
    solve();
}

/** Solve Poisson's equation using the DST. */
void poisson_fft::solve()
{
    fftw_execute(plan);

    double fac = 1./(2*(n+1))*1./(2*(n+1));
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            v[n*i+j] *= fac/(lam[i]+lam[j]);
        }
    }

    fftw_execute(plan);
}

/** Print the solution, padding with zeros for the Dirichlet boundary
 *  conditions. */
void poisson_fft::print_solution()
{
    for (int j=0; j<n+2; ++j) printf("0 ");
    puts("");
    for (int i=0; i<n; ++i) {
        printf("0 ");
        for (int j=0; j<n; ++j) printf("%g ", v[n*i+j]);
        puts("0");
    }
    for (int j=0; j<n+2; ++j) printf("0 ");
    puts("");
}