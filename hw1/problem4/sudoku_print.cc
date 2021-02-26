#include "sudoku_serial.hh"

int main()
{
    SudokuSerial s("test1.sdk");
    s.solve();
    s.print();
}
