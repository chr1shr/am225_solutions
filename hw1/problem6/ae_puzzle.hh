#ifndef AE_PUZZLE_HH
#define AE_PUZZLE_HH

#include <cstdio>
#include <vector>

#include "letter.hh"
#include "adj_table.hh"

const int ae_file_buf_size=1024;

class ae_puzzle {
    public:
        ae_puzzle* master;
        /** The size of the domain in the x direction. */
        int m;
        /** The size of the domain in the y direction. */
        int n;
        /** The size of the domain in the z direction. */
        int o;
        /** The size of the domain in the x direction, with padding added. */
        int mn;
        /** The size of the domain in the y direction, with padding added. */
        int mno;
        /** The total number of letter types. */
        int nletter;
        /** The grid memory structure, containing the boundary padding layers. */
        char *grid;
        /** The shortcut table, only used for pattern solves. */
        int *sc;
        /** The number of each letter type. */
        int *num;
        /** Total number of solutions. */
        int sol;
        ae_puzzle(const char *filename);
        ae_puzzle(ae_puzzle &a);
        ~ae_puzzle();
        void solve(int w);
        void solve_parallel(int w);
        void print_solution();
    private:
        /** The array with the letter information. */
        letter *le;
        adj_table **at;
        std::vector<adj_table*> atmem;
        bool put(char *p,letter *lp,int s,int q);
        void remove(char *p,letter *lp,int s,int q);
        void build_adjacency_table();
        bool same(adj_table *nat,adj_table *&pre);
        int read_pattern(FILE *fp,char *buf);
        inline unsigned int find_mask(int w);
        bool compute_mask(unsigned int &mask,letter *lp,int s,int q,int i,int j,int k);
};

#endif
