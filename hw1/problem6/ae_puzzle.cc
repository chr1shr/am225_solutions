#include "ae_puzzle.hh"
#include "letter.hh"
#include "common.hh"

#include <cstdlib>
#include <cstring>

#ifdef _OPENMP
#include "omp.h"
#endif

/** Sets up the puzzle solver class, reading the puzzle configuration from a
 * file.
 * \param[in] filename the name of the configuration file. */
ae_puzzle::ae_puzzle(const char *filename) : master(NULL), sol(0) {

    // Open the configuration, read in the puzzle dimensions, and the total
    // number of puzzle pieces
    FILE *fp=safe_fopen(filename,"r");
    fscanf(fp,"%d %d %d %d",&m,&n,&o,&nletter);
    if(m<=0||m>128||n<=0||n>128||o<=0||o>128) fatal_error("Grid dimensions out of range");
    if(nletter<=0||nletter>128) fatal_error("Number of letters out of range");
    mn=m*n;mno=mn*o;

    // Read in the puzzle piece information
    le=new letter[nletter];
    int pvol=0;
    char buf[ae_file_buf_size];
    fgets(buf,ae_file_buf_size,fp);
    for(int l=0;l<nletter;l++) {
        if(fgets(buf,ae_file_buf_size,fp)==NULL) fatal_error("Error reading file");
        le[l].read_data(buf);
        pvol+=le[l].num*le[l].len;
    }

    // Set up the puzzle grid. If the 'pattern' keyword is provided, then scan
    // the configuration file for the grid pattern
    int dvol=0;
    grid=new char[mno];
    fgets(buf,ae_file_buf_size,fp);
    if(strcmp(buf,"block\n")==0) {
        for(char *gp=grid;gp<grid+mno;gp++) *gp=0;
        dvol=mno;
        sc=NULL;
    } else if(strcmp(buf,"pattern\n")==0) {
        dvol=read_pattern(fp,buf);
    } else fatal_error("Can't read 'block' or 'pattern' keyword");
    fclose(fp);

    // Print diagnostic information
    printf("# Total blocks in letters: %d\n"
           "# Total blocks in domain: %d\n",pvol,dvol);
    if(pvol!=dvol) fatal_error("Block/domain mismatch");

    // Set up the adjacency table
    at=new adj_table*[mno];
    build_adjacency_table();

    // Allocate memory to store the number of available pieces. Reorganize the
    // letter information to be appropriate for the specific problem grid size.
    num=new int[nletter];
    for(int l=0;l<nletter;l++) {
        le[l].rebase(m,mn);
        num[l]=le[l].num;
    }
}

/** Sets up a slave puzzle solver that has its own solution grid, but makes
 * references to all of the puzzle configuration information in a master
 * solver.
 * \param[in] a a master puzzle solver. */
ae_puzzle::ae_puzzle(ae_puzzle &a) : master(&a), m(a.m), n(a.n), o(a.o),
    mn(a.mn), mno(a.mno), nletter(a.nletter), grid(new char[mno]), sc(a.sc),
    num(new int[nletter]), sol(0), le(a.le), at(a.at) {
    memcpy(grid,a.grid,mno*sizeof(char));
    memcpy(num,a.num,nletter*sizeof(int));
}

/** The class destructor frees the dynamically allocated memory. */
ae_puzzle::~ae_puzzle() {
    delete [] num;
    delete [] grid;
    if(master==NULL) {
        if(sc!=NULL) delete [] sc;
        for(unsigned int k=0;k<atmem.size();k++) delete atmem[k];
        delete [] at;
        delete [] le;
    }
}

/** Reads a grid pattern from a configuration file.
 * \param[in] fp an open file handle to read from.
 * \param[in] buf temporary space for reading in lines from the file. */
int ae_puzzle::read_pattern(FILE *fp,char *buf) {
    int dvol=0;
    char *gp=grid;
    for(int k=0;k<o;k++) {
        for(int j=0;j<n;j++) {

            // Process a line from the file, using '.' entries for empty slots
            // and '+' entries for occupied slots
            fgets(buf,ae_file_buf_size,fp);
            if(strlen(buf)<size_t(m+1)) fatal_error("Pattern line too short");
            for(int i=0;i<m;i++) {
                if(buf[i]=='.') {*(gp++)=0;dvol++;}
                else if(buf[i]=='+') *(gp++)=1;
                else fatal_error("Corrupt pattern information");
            }
        }
        fgets(buf,ae_file_buf_size,fp);
    }

    // Allocate the array for rapid jumping
    sc=new int[mno];
    int p=mno;
    for(int q=mno-1;q>=0;q--) {
        sc[q]=p-q;
        if(grid[q]==0) p=q;
    }
    return dvol;
}

/** Builds the tables of possible pieces that can be placed at each grid
 * location, creating separate lists depending on whether pieces occupy
 * adjacent squares. */
void ae_puzzle::build_adjacency_table() {
    adj_table **ap=at;

    // Loop over the grid locations
    for(int k=0;k<o;k++) for(int j=0;j<n;j++) for(int i=0;i<m;i++,ap++) {
        adj_table *nat=new adj_table;

        //
        unsigned int mask,bmask=(k==o-1?4:0)|(j==n-1?2:0)|(i==m-1?1:0),key;bmask=0;
        int w;
        for(int l=0;l<nletter;l++) for(int s=0;s<le[l].subtypes;s++) for(int q=0;q<le[l].len;q++) {
            mask=bmask;
            if(compute_mask(mask,le+l,s,q,i,j,k)) continue;
            for(w=0;w<8;w++) if((w&mask)==0) nat->c[w]++;
        }
        nat->allocate();

        for(int l=0;l<nletter;l++) for(int s=0;s<le[l].subtypes;s++) for(int q=0;q<le[l].len;q++) {
            mask=bmask;
            if(compute_mask(mask,le+l,s,q,i,j,k)) continue;
            key=(l<<16)|(s<<8)|q;
            for(w=0;w<8;w++) if((w&mask)==0) nat->add(key,w);
        }

        adj_table *pre;
        if(same(nat,pre)) {
            delete nat;
            *ap=pre;
        } else {
            atmem.push_back(nat);
            *ap=nat;
        }
    }
}

bool ae_puzzle::same(adj_table *nat,adj_table *&pre) {
    for(unsigned int z=0;z<atmem.size();z++)
        if(nat->same(atmem[z])) {pre=atmem[z];return true;}
    return false;
}

bool ae_puzzle::compute_mask(unsigned int &mask,letter *lp,int s,int q,int i,int j,int k) {
    int *ds=lp->disp(s),*de=ds+lp->len;
    for(int *dp=ds;dp<de;dp++) {
        int d=*dp-ds[q];
        if(d<0) return true;
        int dd=d+(1+256+65536)*128;
        int ck=((dd>>16)&255)+k-128,cj=((dd>>8)&255)+j-128,ci=(dd&255)+i-128;
        if(ci<0||ci>=m||cj<0||cj>=n||ck<0||ck>=o||grid[ci+cj*m+ck*mn]!=0) return true;
        if(d==1) mask|=1;
        if(d==256) mask|=2;
        if(d==65536) mask|=4;
    }
    return false;
}

inline unsigned int ae_puzzle::find_mask(int w) {
    unsigned int mask=at[w]->rec[0]==at[w]->rec[1]||grid[w+1]!=0?1:0;
    if(at[w]->rec[0]==at[w]->rec[2]||grid[w+m]!=0) mask|=2;
    if(at[w]->rec[0]==at[w]->rec[4]||grid[w+mn]!=0) mask|=4;
    return mask;
}

void ae_puzzle::solve_parallel(int w) {
    if(w==mno) {
        print_solution();
    } else if(grid[w]==0) {
        unsigned int mask=find_mask(w),
                     *rs=at[w]->rec[mask],
                     *re=at[w]->rec[mask+1];
#pragma omp parallel
        {
            ae_puzzle ae(*this);
#pragma omp for schedule(dynamic)
            for(unsigned int *rp=rs;rp<re;rp++) {
                unsigned int key=*rp;
                int l=(key>>16)&255;
                if(ae.num[l]>0) {
                    int s=(key>>8)&255,q=key&255;
                    if(ae.put(ae.grid+w,le+l,s,q)) {
                        ae.num[l]--;
                        ae.solve(w+1);
                        ae.remove(ae.grid+w,le+l,s,q);
                        ae.num[l]++;
                    }
                }
            }
        }
    } else if(grid[w]==1) {
        solve(w+sc[w]);
    } else solve(w+1);
}

void ae_puzzle::solve(int w) {
    if(w==mno) {
#pragma omp critical
        {
            print_solution();
        }
    } else if(grid[w]==0) {
        unsigned int mask=find_mask(w),
                     *rp=at[w]->rec[mask],
                     *re=at[w]->rec[mask+1];
        for(;rp<re;rp++) {
            unsigned int key=*rp;
            int l=(key>>16)&255;
            if(num[l]==0) continue;
            int s=(key>>8)&255,q=key&255;
            if(put(grid+w,le+l,s,q)) {
                num[l]--;
                solve(w+1);
                remove(grid+w,le+l,s,q);
                num[l]++;
            }
        }
    } else if(grid[w]==1) {
        solve(w+sc[w]);
    } else solve(w+1);
}

/** Tries to add a piece to the puzzle grid; if successful, then it updates the
 * grid with the letter.
 * \param[in] p the location in the grid to add the piece.
 * \param[in] lp a pointer to the letter to add.
 * \param[in] s the block index in the letter to place at the location p.
 * \param[in] q the letter orientation.
 * \return True if the piece was successfully added; false otherwise. */
bool ae_puzzle::put(char *p,letter *lp,int s,int q) {
    int *ds=lp->disp(s),*de=ds+lp->len;
    p-=ds[q];

    for(int *dp=ds;dp<de;dp++) if(p[*dp]>0) return false;

    for(int *dp=ds;dp<de;dp++) p[*dp]=lp->character;
    return true;
}

/** Removes a piece to the puzzle grid.
 * \param[in] p the location in the grid to add the piece.
 * \param[in] lp a pointer to the letter to add.
 * \param[in] s the block index in the letter to place at the location p.
 * \param[in] q the letter orientation. */
void ae_puzzle::remove(char *p,letter *lp,int s,int q) {
    int *ds=lp->disp(s),*de=ds+lp->len;p-=ds[q];
    for(int *dp=ds;dp<de;dp++) p[*dp]=0;
}

void ae_puzzle::print_solution() {
#ifdef _OPENMP
    printf("Solution %d (T%d)\n",(master==NULL?sol:master->sol)++,omp_get_thread_num());
#else
    printf("Solution %d\n",(master==NULL?sol:master->sol)++);
#endif
    for(int j=0;j<n;j++) {
        for(int k=0;k<o;k++) {
            for(int i=0;i<m;i++) {
                char g=grid[i+m*j+mn*k];
                putchar(g<32?(g==0?'.':'+'):g);
            }
            putchar(' ');
        }
        putchar('\n');
    }
    putchar('\n');
}
