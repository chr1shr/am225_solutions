#include "rbf_hilbert.hh"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

void rbf_hilbert::matrix_analysis(int bls,int &tot,int &pre) {
    pre=0;tot=0;
	printf("BLS %d\n",bls);
    double dx,dy,rsq;

    // Count the number of non-zero terms
    for(int j=0;j<n;j++) {
        int li=(j/bls)*bls;
        for(int i=0;i<j;i++) {
            dx=px[i]-px[j];
            dy=py[i]-py[j];
            rsq=dx*dx+dy*dy;
            if(rsq<lsq) (i>=li?pre:tot)++;
        }
    }

    pre=n+2*pre;
    tot=pre+2*tot;
}

void rbf_hilbert::init_hilbert() {
    double *p=new double[2*n],*pp=p;
    hil_coord *h=new hil_coord[n];

    for(int i=0;i<n;i++,pp+=2) {
        *pp=urand();
        pp[1]=urand();
        h[i].ind=i;
        h[i].q=pos(*pp,pp[1]);
    }
    std::sort(h,h+n);

    // Create random positions
    pp=p;
    for(int i=0;i<n;i++) {
        pp=p+2*h[i].ind;
        px[i]=*pp;
        py[i]=pp[1];
        rs[i]=pf[i]=exp(-2*(px[i]*px[i]+py[i]*py[i]));
    }
}

unsigned int rbf_hilbert::pos(double x,double y) {
    unsigned int q=0;
    double d=0.5,t;

    for(int i=0;i<12;i++,d*=0.5) {
        q<<=2;
        if(x>0) {x-=d;q|=3;} else x+=d;
        if(y>0) {y-=d;q^=1;} else {
            y+=d;
            t=y;y=q&1?-x:x;x=q&1?-t:t;
        }
    }
    return q;
}
