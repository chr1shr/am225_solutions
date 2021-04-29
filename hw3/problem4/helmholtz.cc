#include <cstdio>
#include <vector>
#include <algorithm>
#include "helmholtz.hh"

extern "C" {
    int zgbsv_(const int *n, const int *kl, const int *ku, const int *nrhs,
               dcomplex *ab, const int *ldab, int *ipiv, dcomplex *b,
               const int *ldb, int *info);
}

Helmholtz::Helmholtz(int n_, double f_, dcomplex* er_) :
    n(n_),
    nn((n-1)*(n-1)),
    ih2(n*n),
    v(reinterpret_cast<dcomplex*>(fftw_alloc_complex(nn))),
    f(f_),
    omega(2*M_PI*f),
    k02(omega*omega*mu0*e0),
    er(er_)
{
    for (int i=0; i<nn; ++i) v[i]=0;
}

Helmholtz::~Helmholtz()
{
    fftw_free(v);
}

/** Set the frequency.
 * \param[in] f_ : The given frequency. */
void Helmholtz::set_frequency(double f_)
{
    f = f_;
    omega = 2*M_PI*f;
    k02 = omega*omega*mu0*e0;
}

/** Print the absolute value of the solution, padding with zeros for the
 *  Dirichlet boundary conditions. */
void Helmholtz::print_solution()
{
    for (int j=0; j<n+1; ++j) printf("0 ");
    puts("");
    for (int i=0; i<n-1; ++i) {
        printf("0 ");
        for (int j=0; j<n-1; ++j) {
            dcomplex u = v[(n-1)*i+j];
            if (er!=nullptr && er[(n+1)*(i+1)+(j+1)].real()<2) u=0;
            printf("%g ", std::abs(u));
        }
        puts("0");
    }
    for (int j=0; j<n+1; ++j) printf("0 ");
}

SlowHelmholtz::SlowHelmholtz(int n_, double f_, dcomplex* er_) :
    Helmholtz(n_, f_, er_)
{}

/** Solve the Helmholtz equation by inverting the banded matrix directly. */
void SlowHelmholtz::solve()
{
    int nsub = n-1;
    int nrhs = 1;
    int ldab = 3*nsub+1;
    std::vector<dcomplex> diags(ldab*nn, 0);
    dcomplex* ab = new dcomplex[ldab*nn];
    int* ipiv = new int[nn];
    int info;

    // The first nsub rows of diags are for internal use by LAPACK
    auto start = diags.begin() + nn*nsub;

    // Construct the diagonals
    if (er == nullptr) {
        std::fill_n(start + nn*nsub, nn, -4*ih2+k02);
    } else {
        auto it = start + nn*nsub;
        for (int i=1; i<n; ++i) {
            for (int j=1; j<n; ++j) {
                *it = dcomplex(-4*ih2,0) + k02*er[(n+1)*i+j];
                ++it;
            }
        }
    }
    for (int d : {n-1, 1, -1, -(n-1)})
        std::fill_n(start + nn*(nsub-d) + d*(d>0), nn-std::abs(d), ih2);
    for (int i=1; i<nn; ++i) {
        if (i%(n-1)==0) {
            *(start+nn*(nsub-1)+i) = 0;
            *(start+nn*(nsub+1)+i-1) = 0;
        }
    }

    // Transpose to column-major for LAPACK
    for (int i=0; i<ldab; ++i)
        for (int j=0; j<nn; ++j)
            ab[ldab*j+i] = diags[nn*i+j];

    // Call LAPACK
    zgbsv_(&nn, &nsub, &nsub, &nrhs, ab, &ldab, ipiv, v, &nn, &info);

    delete [] ab;
    delete [] ipiv;

    if (info < 0) {
        fprintf(stderr, "Invalid argument at position %d\n", -info);
        exit(1);
    } else if (info > 0) {
        fprintf(stderr, "Matrix is singular at position %d\n", info);
        exit(1);
    }
}

FastHelmholtz::FastHelmholtz(int n_, double f_, dcomplex* er_) :
    Helmholtz(n_, f_, er_),
    vtrue(nullptr),
    lam(new double[n-1]),
    temp(fftw_alloc_real(nn)),
    plan(fftw_plan_r2r_2d(n-1, n-1, temp, temp, FFTW_RODFT00,  FFTW_RODFT00, FFTW_MEASURE))
{
    for (int i=0; i<n-1; ++i) {
        double s = std::sin((i+1)*M_PI/2./n);
        lam[i] = -4*ih2*s*s + k02/2.;
    }
}

FastHelmholtz::~FastHelmholtz()
{
    fftw_destroy_plan(plan);
    fftw_free(temp);
    delete [] lam;
}

/** Perform a 2D DST-I. */
void FastHelmholtz::dst1_2d(dcomplex* x)
{
    // The data is complex, so perform two real DSTs using the fact that
    // dst(a+bi) = dst(a) + dst(b)i.
    for(int i=0; i<nn; ++i) temp[i] = x[i].real();
    fftw_execute(plan);
    for(int i=0; i<nn; ++i) x[i].real(temp[i]);
    for(int i=0; i<nn; ++i) temp[i] = x[i].imag();
    fftw_execute(plan);
    for(int i=0; i<nn; ++i) x[i].imag(temp[i]);
}

/** Solve the Helmholtz equation, deciding which method to use. */
void FastHelmholtz::solve()
{
    if (er != nullptr) {
        iter_solve();
    } else {
        direct_solve();
    }
}

/** Solve the Helmholtz equation using the DST-based direct method. */
void FastHelmholtz::direct_solve()
{
    dst1_2d(v);

    double fac = 1./(2*n)*1./(2*n);
    for (int i=0; i<n-1; i++) {
        for (int j=0; j<n-1; j++) {
            v[(n-1)*i+j] *= fac/(lam[i]+lam[j]);
        }
    }

    dst1_2d(v);
}

/** Solve the Helmholtz equation using the iterative method.
 * \param[in] niter : The number of iterations to perform. */
void FastHelmholtz::iter_solve(int niter)
{
    dcomplex* rhs = new dcomplex[(n-1)*(n-1)];
    for (int i=0; i<nn; ++i) rhs[i] = v[i];

    // Calculate the norm of the true solution
    double rel = 0;
    if (vtrue != nullptr) {
        for (int j=0; j<nn; ++j) {
            rel = std::max(rel, std::abs(vtrue[j]));
        }
    }

    for (int i=0; i<niter; ++i) {
        for (int j=1, l=0; j<n; ++j) {
            for (int k=1; k<n; ++k) {
                v[l] = rhs[l] - (k02*er[(n+1)*j+k]-k02)*v[l];
                l++;
            }
        }
        direct_solve();

        // Calculate the relative error
        if (vtrue != nullptr) {
            double norm = 0;
            for (int j=0; j<nn; ++j) {
                norm = std::max(norm, std::abs(v[j]-vtrue[j]));
            }
            printf("%i %g\n", i, norm/rel);
        }
    }

    delete [] rhs;
}
