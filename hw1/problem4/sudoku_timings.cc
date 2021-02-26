#include <cstdlib>
#include <string>
#include <sstream>
#include "sudoku_parallel.hh"

#ifdef _OPENMP
#include "omp.h"
inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
inline double wtime() {return double(std::clock())*(1./CLOCKS_PER_SEC);}
#endif

int main(int argc, char* argv[])
{
    double t[3][3];
    int N[3] = {1000, 100, 3};
    for (int i=0; i<3; ++i) {
        std::stringstream ss;
        ss << "test" << i+1 << ".sdk";
        std::printf("Puzzle %d:\n", i+1);
        int l=0;
        for (int j : {1,2,4}) {
            omp_set_num_threads(j);
            double t0 = wtime();
            for (int k=0; k<N[i]; ++k) {
                SudokuParallel s(ss.str());
                s.solve();
            }
            t[i][l] = (wtime()-t0)/N[i];
            std::printf("%g / %g\n", t[i][l], t[i][0]/(j*t[i][l]));
            l++;
        }
        std::puts("\n");
    }
}
