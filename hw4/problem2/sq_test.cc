#include <cstdio>
#include <cmath>

#include "squircle.hh"

int main() {
    squircle sq(40,8);
    sq.set_source(1);
    sq.solve();
    sq.output("x.sol");
    printf("%g\n",sq.l2_norm_mms(1));
}
