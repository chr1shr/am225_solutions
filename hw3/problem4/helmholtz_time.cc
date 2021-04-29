#include <cstdio>
#include <cmath>
#include <chrono>
#include "helmholtz.hh"

int main(int argc, char* argv[])
{
    double f = 21.3e6;
    for (int n=16; n<=256; n*=2) {
        int x = n/2-1;
        auto start = std::chrono::system_clock::now();
        SlowHelmholtz hh(n, f);
        hh.v[(n-1)*x+x] = hh.ih2;
        hh.solve();
        std::chrono::duration<double> diff = std::chrono::system_clock::now()-start;
        printf("%i %g\n", n, diff.count());
    }
    puts("\n");
    for (int n=16; n<=1024; n*=2) {
        int x = n/2-1;
        auto start = std::chrono::system_clock::now();
        FastHelmholtz hh(n, f);
        hh.v[(n-1)*x+x] = hh.ih2;
        hh.solve();
        std::chrono::duration<double> diff = std::chrono::system_clock::now()-start;
        printf("%i %g\n", n, diff.count());
    }
}
