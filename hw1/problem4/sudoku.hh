#ifndef SUDOKU_HH
#define SUDOKU_HH

#include <string>
#include <array>

/** A Sudoku solver. */
class Sudoku {
    public:
        typedef std::array<int,81> Board;
        Sudoku();
        Sudoku(const Board& board);
        Sudoku(const std::string& filename);
        void load(const std::string& filename);
        virtual int solve(int i=0) = 0;
    protected:
        int look(int i) const;
        void init_indices();
        int total;
        Board board;
        int index[243];
};

#endif
