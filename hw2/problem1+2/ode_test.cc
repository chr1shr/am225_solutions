#include "sol_fsal.hh"
#include "brusselator.hh"
#include "osc.hh"
#include <cstdio>

int main() {

    // Solve the oscillator problem using the sixth-order
    // Richardson-extrapolated Cash-Karp method
//    osc_ckr bh;
brus_ckr bh;
    double t=20, lamda=1e-3;
    int d_steps=1200;

    // Output the integration time points
    bh.solve(t,lamda,lamda,true);
    puts("\n");

    // Do dense output
    bh.t=0;
    bh.solve(t,lamda,lamda,false,d_steps);
}
