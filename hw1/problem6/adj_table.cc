#include "adj_table.hh"

#include <cstdio>

adj_table::adj_table() {
    rec[0]=NULL;
    for(int l=0;l<8;l++) c[l]=0;
}

adj_table::~adj_table() {
    if(rec[0]!=NULL) delete [] rec[0];
}

void adj_table::allocate() {
    int tot=total();
    rec[0]=new unsigned int[tot];
    for(int l=0;l<8;l++) {
        rec[l+1]=rec[l]+c[l];
        c[l]=0;
    }
}

int adj_table::total() {
    int tot=*c;
    for(int l=1;l<8;l++) tot+=c[l];
    return tot;
}

bool adj_table::same(adj_table* ap) {
    for(int k=0;k<8;k++) if(c[k]!=ap->c[k]) return false;
    for(int k=0;k<*c;k++) if(rec[0][k]!=ap->rec[0][k]) return false;
    return true;
}
