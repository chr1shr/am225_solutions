#include <cstdio>
#include <cmath>
#include "helmholtz.hh"

int main(int argc, char* argv[])
{
    if (argc != 2) {
        puts("Usage: ./helmholtz_test <frequency>");
        exit(0);
    }

    int n = 256;
    int nn = (n+1)*(n+1);
    double f = atof(argv[1]);
    int x = std::floor(n*0.6)-1, y = std::floor(n*0.7)-1;

    // Read the MRI data
    dcomplex* er = new dcomplex[nn];
    FILE* file = fopen("mri_data_256.dat", "rb");
    fread(er, sizeof(dcomplex), nn, file);
    fclose(file);

    SlowHelmholtz slow(n, f, er);
    slow.v[(n-1)*y+x] = slow.ih2;
    slow.solve();
    slow.print_solution();

    puts("\n\n");
    FastHelmholtz fast(n, f, er);
    fast.v[(n-1)*y+x] = fast.ih2;
    //fast.vtrue = slow.v;
    fast.solve();
    fast.print_solution();

    delete [] er;
}
