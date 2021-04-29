#include <cmath>

#include "lax_wend.hh"

int main(int argc,char **argv) {

	if (argc<2) {
		fputs("Usage: ./lw_solve <grid size> [opt: filename]\n",stderr);
		exit(-1);
	}

    // Period of the solution
    const double T=3*M_PI/sqrt(5);

    // Number of snapshots to output, and iterations between snapshot
    const int snaps=4;

    // Number of gridpoints
    int m=atoi(argv[1]);

    // Integration timestep safety factor
    const double sf=1./3;

    // Create the diffusion simulation. Initialize the solution and the nu
    // table.
    lax_wend lw(m);
    lw.init_exp_sine();
    lw.init_triangle();

    // Integrate and save the solution snapshots to file
    lw.solve(argc==2?"lw3.out":argv[2],snaps,T,sf);
}
