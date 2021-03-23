#include "sol_sun.hh"
#include "brusselator.hh"

int main() {

    // Solve up to x=20 using the HH method with 500 steps
    brus_sun bs;
    bs.solve_fixed(20.,500,true);
}
