#ifndef POISSON_SCHUR_HH
#define POISSON_SCHUR_HH

#include <vector>
#include <functional>
#include <unordered_map>
#include "conj_grad.hh"
#include "subdomain.hh"

class subdomain;

class poisson_schur : public conj_grad
{
    public:
        poisson_schur(int n_, const char* filename);
        ~poisson_schur();
        void solve();
        void solve(const std::function<double(double,double)>& fun);
        virtual void mul_A(double* in, double* out);
        double v(int i, int j, bool glue = false);
        void print_solution(bool glue = false);
        /** The number of interior gridpoints in one dimension */
        const int n;
        /** The total number of interior gridpoints */
        const int nn;
        /** The grid spacing */
        const double h;
        /** The inverse grid spacing squared */
        const double ih2;
        /** The source term */
        double* const f;
    private:
        void init_subdomains(const char* filename);
        /** The number of gridpoints on the glue */
        int ng;
        /** Map */
        int* c;
        /** Glue indices */
        std::unordered_map<int,int> g2l;
        /** The subdomains */
        std::vector<subdomain*> subdomains;
};

#endif
