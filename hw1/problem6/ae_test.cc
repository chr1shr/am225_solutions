#include "ae_puzzle.hh"

int main() {
    ae_puzzle aep("ae_2x8x10.txt");
    aep.solve_parallel(0);
}
