#include <cstdio>
#include <cstdlib>

#include "c_automaton.hh"
#include "custom_rng.hh"

c_automaton::c_automaton(int m_,int n_,unsigned int survive_,unsigned int birth_) :
    m(m_), n(n_), mn(m*n), survive(survive_), birth(birth_),
    r(new bool[mn]), s(new bool[mn]) {}

c_automaton::~c_automaton() {
    delete [] s;
    delete [] r;
}

void c_automaton::clear() {
    for(bool *rp=r;rp<r+mn;rp++) *rp=false;
}

void c_automaton::set_random(int li,int ui,int lj,int uj,float prop) {
    int thresh=int(RAND_MAX*prop);
    for(int j=lj;j<uj;j++) for(int i=li;i<ui;i++)
        r[i+m*j]=rand()<thresh;
}

void c_automaton::step() {

#pragma omp parallel for
    for(int j=0;j<n;j++) {
        bool *rp=r+m*j,*sp=s+m*j,
             *rd=j>0?rp-m:r+(mn-m),
             *ru=j<n-1?rp+m:r;
        for(int i=0;i<m;i++) {
            int il=i<m-1?i+1:0,ir=i>0?i-1:m-1,c=rd[il]?1:0;
            if(rd[i]) c++;
            if(rd[ir]) c++;
            if(rp[il]) c++;
            if(rp[ir]) c++;
            if(ru[il]) c++;
            if(ru[i]) c++;
            if(ru[ir]) c++;

            sp[i]=(rp[i]?survive:birth)&(1<<c);
        }
    }
    bool *q=r;r=s;s=q;
}

int c_automaton::count() {
    int sum=0;

#pragma omp parallel for reduction(+:sum)
    for(bool *rp=r;rp<r+mn;rp++) sum+=*rp;
    return sum;
}

void c_automaton::ascii_art() {
    char *buf=new char[m+1],*bp;
    buf[m]=0;

    bool *rp=r;
    for(int j=0;j<n;j++) {
        for(bp=buf;bp<buf+m;bp++,rp++)
            *bp=*rp?'#':' ';
        puts(buf);
    }
}

void c_automaton::time_steps(int steps) {
    double t0=-wtime();
    for(int k=0;k<steps;k++) step();
    t0+=wtime();

    printf("# Time for %d step%s : %g s\n# Time per step : %g ms\n",
            steps,steps!=1?"s":"",t0,1e3*t0/steps);
}
