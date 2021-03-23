#include <iostream>
#include <cstdlib>
#include <cstring>

/** \brief Compute the next unique tree using the Beyer-Hedetniemi algorithm.
 *  \param[in,out] level : Levels of each node
 *  \param[in,out] parent : Parents of each node
 *  \param[in] p : Position of rightmost element with level > 1
 *  \return Flag indicating there are more trees remaining */
bool next_tree(int* level, int* parent, int p, int n)
{
    while (level[p]==1) p--;
    if (p==0) {
        return false;
    } else if (level[p]==2 && level[p-1]==1) {
        level[p] = 1;
        parent[p] = parent[p-1];
    } else {
        int q = p-parent[p]+1;
        for (int i=p; i<n; ++i) {
            level[i] = level[i-q];
            parent[i] = (q+parent[i-q]<p+1) ? parent[i-q] : q+parent[i-q];
        }
        p = n-1;
    }
    return true;
}

/** \brief Print the parent list of all trees of order n. */
void print_trees(int n)
{
    int* level  = new int[n];
    int* parent = new int[n];

    // Initialize the tree
    int p = n-1;
    for (int i=0; i<n; ++i) {
        level[i] = i;
        parent[i] = i;
    }

    do {
        for (int i=0; i<n; ++i) std::cout << parent[i] << " ";
        std::cout << std::endl;
    } while (next_tree(level, parent, p, n));

    delete [] level;
    delete [] parent;
}

/** \brief Count the number of trees up to order n.
 *  \note This is based on the mathematical formula in 10.1137/0209055. */
void count_trees(int n)
{
    int* counts = new int[n];

    // Count from the bottom up
    counts[0] = 1;
    std::cout << counts[0] << std::endl;
    for (int i=1; i<n; ++i) {
        counts[i] = 0;
        for (int j=0; j<i; ++j) {
            int s=0;
            for (int k=1; k<=i/(j+1); ++k) s += counts[i-k*(j+1)];
            counts[i] += (j+1)*counts[j]*s;
        }
        counts[i] /= i;
        std::cout << counts[i] << std::endl;
    }

    delete [] counts;
}

int main(int argc, char* argv[])
{
    if (argc == 2) {
        int n = atoi(argv[1]);
        count_trees(n);
    } else if (argc == 3 && strcmp(argv[1], "--print") == 0) {
        int n = atoi(argv[2]);
        print_trees(n);
    } else {
        std::cerr << "Usage: ./bh_trees [--print] <order>" << std::endl;
        exit(1);
    }
}