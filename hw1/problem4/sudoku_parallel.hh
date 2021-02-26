#ifndef SUDOKU_PARALLEL_HH
#define SUDOKU_PARALLEL_HH

#include <deque>
#include "sudoku.hh"
#include "sudoku_serial.hh"

/** A parallel Sudoku solver. */
class SudokuParallel : public Sudoku {
    public:
        SudokuParallel() : Sudoku() {}
        SudokuParallel(const Board& board) : Sudoku(board) { analyze(); }
        SudokuParallel(const std::string& filename) : Sudoku(filename) { analyze(); }
        virtual int solve(int i=0);
    protected:
        void analyze();
        void generate(std::deque<SudokuSerial>& starts);
        /** Number of generations of bootstrapping to perform.
         * \note Optimal value determined experimentally. */
        int ngen;
        /** Sequence of boards used to bootstrap parallelization. */
        std::deque<SudokuSerial> starts;
};

#endif
