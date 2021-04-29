#include <cstdio>
#include <cmath>
#include <string>
#include <algorithm>
#include "poisson_fft.hh"
#include "poisson_schur.hh"

int main(int argc, char* argv[])
{
    int n = 112;
    auto f = [] (double x, double y) { return std::exp(x-y); };
    poisson_fft p1(n);
    p1.solve(f);
    poisson_schur p2(n, "pss.txt");
    p2.solve(f);

    p1.print_solution();
    puts("\n");
    p2.print_solution();
    puts("\n");
    p2.print_solution(true);
    puts("\n");

    // Compute error
    double err, norm = 0;
    for (int j=0; j<n+1; ++j) printf("0 ");
    puts("");
    for (int i=0,k=0; i<p1.n; ++i) {
        printf("0 ");
        for (int j=0; j<p1.n; ++j,++k) {
            err = std::abs(p1.v[k]-p2.v(i,j));
            norm = std::max(norm, err);
            printf("%g ", err);
        }
        puts("0");
    }
    for (int j=0; j<n+1; ++j) printf("0 ");
    puts("");
    printf("# Maximum norm error: %g\n", norm);
}
