#ifndef STRASSEN_HH
#define STRASSEN_HH

void set_full(int m,double *q,double *c);
void add_full(int m,double *q,double *c);
void sub_full(int m,double *q,double *c);
void set_half(int m,double *s,double *x);
void add_half(int m,double *s,double *x,double *y);
void sub_half(int m,double *s,double *x,double *y);
void strassen(int n,double *a,double *b,double *c);

#endif
