#include <iostream>
#include <cstdlib>
#include <string>
#include "sudoku_serial.hh"
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
    if (argc != 3) {
        std::cerr << "Syntax: ./sudoku_run <serial/parallel> <input_file>" << std::endl;
        return 1;
    }

    std::string execution(argv[1]);
    std::string filename(argv[2]);
    int total;
    double t0 = wtime();
    if (execution == "serial") {
        SudokuSerial s(filename);
        total = s.solve();
    } else if (execution == "parallel") {
        SudokuParallel s(filename);
        total = s.solve();
    } else {
        std::cerr << "Unknown argument: " << execution << std::endl;
        return 1;
    }

    std::printf("Found %d solution%s in %g seconds\n", total, total==1?"":"s", wtime()-t0);
}
