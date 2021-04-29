#ifndef SUBDOMAIN_HH
#define SUBDOMAIN_HH

#include <vector>
#include <functional>
#include <unordered_map>
#include "poisson_fft.hh"
#include "poisson_schur.hh"

// Forward declaration
class poisson_schur;

class subdomain : public poisson_fft
{
    public:
        subdomain(poisson_schur* schur_, int x_, int y_, int s_);
        void solve();
        void solve(const std::function<double(double,double)>& fun);
        const int x;
        const int y;
        const int s;
        const double scl;
        std::unordered_map<int,int> glue;
    private:
        poisson_schur* const schur;
};

#endif
