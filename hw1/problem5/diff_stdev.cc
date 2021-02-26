#include "diffuse.hh"

#include <cstdio>

int main() {
    //diffuse di(400,1.,0,1.);
    diffuse di(400,1,-1.,1.);
    di.num_t=2;
    di.init();
    di.integrate(0,true);
}
