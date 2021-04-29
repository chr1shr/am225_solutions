#include "strassen.hh"

void set_full(int m,double *q,double *c) {
    for(double *qe=q+m*m;q<qe;c+=m) for(double *qr=q+m;q<qr;q++,c++)
        *c=*q;
}

void add_full(int m,double *q,double *c) {
    for(double *qe=q+m*m;q<qe;c+=m) for(double *qr=q+m;q<qr;q++,c++)
        *c+=*q;
}

void sub_full(int m,double *q,double *c) {
    for(double *qe=q+m*m;q<qe;c+=m) for(double *qr=q+m;q<qr;q++,c++)
        *c-=*q;
}

void set_half(int m,double *s,double *x) {
    for(double *se=s+m*m;s<se;x+=m) for(double *sr=s+m;s<sr;s++,x++)
        *s=*x;
}

void add_half(int m,double *s,double *x,double *y) {
    for(double *se=s+m*m;s<se;x+=m,y+=m) for(double *sr=s+m;s<sr;s++,x++,y++)
        *s=*x+*y;
}

void sub_half(int m,double *s,double *x,double *y) {
    for(double *se=s+m*m;s<se;x+=m,y+=m) for(double *sr=s+m;s<sr;s++,x++,y++)
        *s=*x-*y;
}

void strassen(int n,double *a,double *b,double *c) {
    if(n<=2) {
        if(n==1) *c=*a*(*b);
        else {
            *c=*a*(*b)+a[2]*b[1];
            c[1]=a[1]*(*b)+a[3]*b[1];
            c[2]=*a*b[2]+a[2]*b[3];
            c[3]=a[1]*b[2]+a[3]*b[3];
        }
        return;
    }
    int m=n>>1;
    double *s=new double[3*m*m],*t=s+m*m,*q=t+m*m;
    double *a10=a+m,*a01=a+2*m*m,*a11=a01+m;
    double *b10=b+m,*b01=b+2*m*m,*b11=b01+m;
    double *c10=c+m,*c01=c+2*m*m,*c11=c01+m;

    // Q0 contribution
    add_half(m,s,a,a11);add_half(m,t,b,b11);
    strassen(m,s,t,q);set_full(m,q,c);set_full(m,q,c11);

    // Q1 contribution
    add_half(m,s,a10,a11);set_half(m,t,b);
    strassen(m,s,t,q);set_full(m,q,c10);sub_full(m,q,c11);

    // Q2 contribution
    set_half(m,s,a);sub_half(m,t,b01,b11);
    strassen(m,s,t,q);set_full(m,q,c01);add_full(m,q,c11);

    // Q3 contribution
    set_half(m,s,a11);sub_half(m,t,b10,b);
    strassen(m,s,t,q);add_full(m,q,c);add_full(m,q,c10);

    // Q4 contribution
    add_half(m,s,a,a01);set_half(m,t,b11);
    strassen(m,s,t,q);sub_full(m,q,c);add_full(m,q,c01);

    // Q5 contribution
    sub_half(m,s,a10,a);add_half(m,t,b,b01);
    strassen(m,s,t,q);add_full(m,q,c11);

    // Q6 contribution
    sub_half(m,s,a01,a11);add_half(m,t,b10,b11);
    strassen(m,s,t,q);add_full(m,q,c);

    delete [] s;
}
