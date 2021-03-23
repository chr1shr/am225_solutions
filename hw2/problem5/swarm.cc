#include "swarm.hh"

#include <cstdio>
#include <cmath>
#include <sys/stat.h>
#include <sys/types.h>

swarm::swarm(int N_,double J_,double K_,const char* filename_) :
    N(N_), J(J_), K(K_), filename(filename_), d_o_frame(0) {

    // Create directory for output frames
    mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

void swarm::swarm_ff(double t_,double *in,double *out) {
    const double iN=1./N,KiN=K*iN;

    for(int i=0;i<N;i++) {
        double *ini=in+3*i,*outi=out+3*i;

        *outi=0;
        outi[1]=0;
        outi[2]=0.1;

        for(int j=0;j<i;j++) {
            double *inj=in+3*j,*outj=out+3*j,
                   dx=*inj-*ini,
                   dy=inj[1]-ini[1],
                   dth=inj[2]-ini[2],
                   irsq=1/(dx*dx+dy*dy),
                   ir=sqrt(irsq),
                   a=iN*((1+J*cos(dth))*ir-irsq),
                   vx=dx*a,vy=dy*a,
                   vth=KiN*sin(dth)*ir;

            *outi+=vx;*outj-=vx;
            outi[1]+=vy;outj[1]-=vy;
            outi[2]+=vth;outj[2]-=vth;
        }
    }
}

void swarm::swarm_init(double *q) {
    double x,y;

    for(int i=0;i<N;i++) {
        do {
            x=rnd(-1,1);
            y=rnd(-1,1);
        } while(x*x+y*y>1);

        *(q++)=x;*(q++)=y;
        *(q++)=rnd(0,2*M_PI);
    }
}

void swarm::swarm_print_dense(double t_,double *in) {
    char buf[256];

    sprintf(buf,"%s/fr.%d",filename,d_o_frame++);
    FILE *f=fopen(buf,"w");
    if(f==NULL) {
        fputs("Error opening output file\n",stderr);
        exit(1);
    }

    fprintf(f,"# Time: %.10g\n",t_);
    for(double *p=in;p<in+3*N;p+=3)
        fprintf(f,"%.8g %.8g %.8g\n",*p,p[1],p[2]);

    fclose(f);
}
