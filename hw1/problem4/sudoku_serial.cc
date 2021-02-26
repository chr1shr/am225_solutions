#include <iostream>
#include "sudoku_serial.hh"

/** Solve the Sudoku.
 * \param[in] i : The starting position.
 * \return The number of solutions found. */
int SudokuSerial::solve(int i)
{
    if (i==81) {
        // The grid is complete
        total++;
    } else if (board[i]==0) {
        // Insert each possible digit here
        int seen = look(i);
        for (int j=1; j<=9; ++j) {
            if ((seen&(1<<j))==0) {
                board[i] = j;
                solve(i+1);
            }
        }
        board[i] = 0;
    } else {
        // Move on to the next square
        solve(i+1);
    }
    return total;
}

/** Print the board. */
void SudokuSerial::print()
{
    for (int i=0; i<9; ++i) {
        for (int j=0; j<9; ++j) {
            std::cout << board[9*i+j] << " ";
        }
        std::cout << std::endl;
    }
}