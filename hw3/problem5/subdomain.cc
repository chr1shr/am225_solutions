#include "subdomain.hh"

subdomain::subdomain(poisson_schur* schur_, int x_, int y_, int s_) :
    poisson_fft(s_), x(x_), y(y_), s(s_), scl(s*s*schur_->h*schur_->h), schur(schur_)
{
    // Build the glue-to-interior lookup table
    int sn = schur->n;
    if (x>0)    for (int i=0; i<n; ++i) glue[sn*(y+i)+(x-1)] = n*i;       // Left
    if (x+n<sn) for (int i=0; i<n; ++i) glue[sn*(y+i)+(x+n)] = n*i+(n-1); // Right
    if (y>0)    for (int j=0; j<n; ++j) glue[sn*(y-1)+(x+j)] = j;         // Bottom
    if (y+n<sn) for (int j=0; j<n; ++j) glue[sn*(y+n)+(x+j)] = n*(n-1)+j; // Top
}

void subdomain::solve()
{
    this->poisson_fft::solve();
    for (int k=0; k<nn; ++k) v[k]*=scl;
}

void subdomain::solve(const std::function<double(double,double)>& fun)
{
    this->poisson_fft::solve(fun);
    for (int k=0; k<nn; ++k) v[k]*=scl;
}
