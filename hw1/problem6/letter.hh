#ifndef LETTER_HH
#define LETTER_HH

#include <cstdio>
#include <vector>

class letter {
    public:
        /** The array displacements that describe the letter shape. All
         * orientations of the letter are packed contiguously into this array.
         * */
        int *d;
        /** The current number of available pieces of this type. */
        int num;
        /** The number of distinct orientation types for this letter. */
        int subtypes;
        /** The number of blocks that make up this letter. */
        int len;
        /** The character that this wooden letter represents. */
        char character;
        letter() : d(NULL) {}
        ~letter();
        void read_data(char *data);
        void rebase(int m,int mn);
        /** Returns a pointer to the array displacements for a given
         * orientation type.
         * \param[in] s the orientation type. */
        int *disp(int s) {
            return d+len*s;
        }
    private:
        void gen_sym_disp(int *&dp,int xmul,int ymul,std::vector<int> &v,int sym);
        void gen_disp(int *&dp,int xmul,int ymul,std::vector<int> &v);
};

#endif
