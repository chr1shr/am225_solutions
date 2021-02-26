#include "sudoku_parallel.hh"

/** Solve the Sudoku.
 * \param[in] i : The starting position.
 * \return The number of solutions found. */
int SudokuParallel::solve(int i)
{
    // Generate starting boards via breadth-first search
    starts.emplace_back(board);
    for (int i=0; i<ngen; ++i) generate(starts);

    // Solve starting from each of these boards in parallel
    int size = starts.size();
    #pragma omp parallel for reduction(+:total)
    for (int i=0; i<size; ++i) total += starts[i].solve();

    return total;
}

/** Analyze the board to determine the amount of bootstrapping to do. */
void SudokuParallel::analyze()
{
    // Count the number of zeros in the board
    int nz = 0;
    for (int i=0; i<81; ++i) if (board[i]==0) nz++;
    // Interpolate between some magic numbers
    ngen = 4+(4.-13.)/(53.-64.)*(nz-53);
}

/** Generate a new family of boards by performing one step of breadth-first
 *  search on each of the given boards.
 * \param[in]  starts : A list of boards to fill in.
 * \param[out] starts : The new boards generated. */
void SudokuParallel::generate(std::deque<SudokuSerial>& starts)
{
    int size = starts.size();
    for (int j=0; j<size; ++j) {
        const SudokuSerial& s = starts.front();
        int gen = 0;
        for (int i=0; i<81; ++i) {
            if (s.board[i]==0) {
                // Insert each possible digit here
                int seen = s.look(i);
                for (int j=1; j<=9; ++j) {
                    if ((seen&(1<<j))==0) {
                        starts.push_back(s);
                        starts.back().board[i] = j;
                        gen++;
                    }
                }
                if (gen>0) break;
            }
        }
        if (gen!=0) starts.pop_front();
    }
}
