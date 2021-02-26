#include <iostream>
#include <fstream>
#include <sstream>
#include "sudoku.hh"

Sudoku::Sudoku() :
    total(0)
{
    init_indices();
}

Sudoku::Sudoku(const Board& board) :
    total(0),
    board(board)
{
    init_indices();
}

Sudoku::Sudoku(const std::string& filename) :
    total(0)
{
    init_indices();
    load(filename);
}

/** Load a Sudoku board from a file.
 * \param[in] filename : The file to load. */
void Sudoku::load(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file " << filename << std::endl;
        exit(1);
    } else {
        std::string line;
        int i=0;
        while (std::getline(file, line)) {
            std::istringstream ss(line);
            int j=0;
            while (ss >> board[9*i+j]) j++;
            if (++i>9 || j>9) {
                printf("%d %d\n", i, j);
                std::cerr << "Error: malformed file" << std::endl;
                exit(1);
            }
        }
    }
}

/** Accumulate the entries seen in this row, column, and box.
 * \param[in] i : The square defining the row, column, and box.
 * \return A bitmask whose n-th bit is 1 if n was seen. */
int Sudoku::look(int i) const
{
    int seen=0, row=index[3*i], col=index[3*i+1], box=index[3*i+2];
    for (int c=row; c<row+9; ++c)   seen |= 1<<board[c];
    for (int r=col; r<col+81; r+=9) seen |= 1<<board[r];
    for (int r=box; r<box+27; r+=9) {
        for (int c=r; c<r+3; ++c) {
            seen |= 1<<board[c];
        }
    }
    return seen;
}

/** Precompute indices for rows, columns, and boxes. */
void Sudoku::init_indices()
{
    for(int i=0; i<81; i++) {
        index[3*i]   = i/9*9;               // Row index
        index[3*i+1] = i%9;                 // Column index
        index[3*i+2] = i/27*27 + (i%9)/3*3; // Box row index + box column index
    }
}
