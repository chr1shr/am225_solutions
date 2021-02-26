#include "letter.hh"
#include "common.hh"

#include <cstdlib>
#include <cstring>

/** The class destructor removes the dynamically allocated memory for the
 * letter displacements (assuming that it had been allocated). */
letter::~letter() {
    if(d!=NULL) delete [] d;
}

/** Sets up the data structure based on a line of information from an input data file.
 * \param[in] buf the line of information.
 * \param[in] m the memory step in the y direction.
 * \param[in] mn the memory step in the z direction. */
void letter::read_data(char *buf) {

    // Read in the alphabetic character that this letter represents
    char *bp;
    bp=strtok(buf," \t\n");
    if(strlen(bp)>1) fatal_error("Error reading character");
    character=*bp;

    // Read in the total number of letters of this type
    bp=strtok(NULL," \t\n");
    num=atoi(bp);
    if(num<=0||num>128) fatal_error("Error reading number of pieces");

    // Read in the symmetry type for this letter
    bp=strtok(NULL," \t\n");
    int sym=atoi(bp);
    if(sym<0||sym>3) fatal_error("Error reading symmetry type");

    // Set the total number of orientations based on the symmetry type
    const int s_table[4]={3,12,12,24};
    subtypes=s_table[sym];

    // Read in the (x,y) coordinates of the blocks making up this letter
    std::vector<int> v;
    bp=strtok(NULL," \t\n");
    if(bp==NULL) fatal_error("Error reading locations");
    do {
        v.push_back(atoi(bp));
        bp=strtok(NULL," \t\n");
        if(bp==NULL) fatal_error("Error reading locations");
        v.push_back(atoi(bp));
        bp=strtok(NULL," \t\n");
    } while(bp!=NULL);

    // Allocate memory for the array displacements for all of the different
    // orientations of this letter
    len=v.size()/2;
    d=new int[len*subtypes];

    // Generate the different array displacements for when the letter is
    // arranged in the xy plane, the yz plane, or the xz plane
    int *dp=d;
    gen_sym_disp(dp,1,256,v,sym);
    gen_sym_disp(dp,256,65536,v,sym);
    gen_sym_disp(dp,65536,1,v,sym);
}

void letter::gen_sym_disp(int *&dp,int xmul,int ymul,std::vector<int> &v,int sym) {
    gen_disp(dp,xmul,ymul,v);
    switch(sym) {
        case 1:
            gen_disp(dp,ymul,-xmul,v);
            gen_disp(dp,-xmul,-ymul,v);
            gen_disp(dp,-ymul,xmul,v);break;
        case 3:
            gen_disp(dp,-xmul,-ymul,v);
            gen_disp(dp,-ymul,-xmul,v);
            gen_disp(dp,xmul,-ymul,v);
            gen_disp(dp,ymul,-xmul,v);
        case 2:
            gen_disp(dp,ymul,xmul,v);
            gen_disp(dp,-xmul,ymul,v);
            gen_disp(dp,-ymul,xmul,v);
        case 0:
            break;
    }
}

void letter::rebase(int m,int mn) {
    for(int *dp=d;dp<d+len*subtypes;dp++) {
        int dd=*dp+(1+256+65536)*128,
            k=((dd>>16)&255)-128,j=((dd>>8)&255)-128,i=(dd&255)-128;
        *dp=i+m*j+mn*k;
    }
}

void letter::gen_disp(int *&dp,int xmul,int ymul,std::vector<int> &v) {
    for(int i=0;i<2*len;i+=2) *(dp++)=xmul*v[i]+ymul*v[i+1];
}
