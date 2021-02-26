#ifndef SUDOKU_SERIAL_HH
#define SUDOKU_SERIAL_HH

#include "sudoku.hh"

/** A serial Sudoku solver. */
class SudokuSerial : public Sudoku {
    public:
        SudokuSerial() : Sudoku() {}
        SudokuSerial(const Board& board) : Sudoku(board) {}
        SudokuSerial(const std::string& filename) : Sudoku(filename) {}
        virtual int solve(int i=0);
        void print();
        friend class SudokuParallel;
};

#endif
