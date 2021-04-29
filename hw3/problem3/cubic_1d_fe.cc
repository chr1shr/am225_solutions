#include <cstdlib>

#include "cubic_1d_fe.hh"
#include "blas.h"

/** Initializes the source function to be a constant. */
void cubic_1d_fe::init_const() {
    for(int i=0;i<2*n+2;i++) f[i]=i&1?0:1;
    assemble_b();
}

/** Initializes the source function to be a linear slope, f(x)=1.5-x. */
void cubic_1d_fe::init_slope() {
    for(int i=0;i<2*n+2;i++) f[i]=i&1?-h:0.5-i*h;
    assemble_b();
}

/** Initializes the source function so that the solution will match
 * a manufactured solution, u(x)=exp(1-x)*sin(5*pi*x). */
void cubic_1d_fe::init_mms() {
    const double o=5*M_PI;
    int i=0;
    double xx=1,co,so;
    for(i=0;i<2*n+2;i+=2,xx+=h) {
        co=cos(o*xx);so=sin(o*xx);
        f[i]=exp(1-xx)*(o*(1-2*xx)*co+((1-o*o)*xx-1)*so);
        f[i+1]=-h*exp(1-xx)*(o*(4+(o*o-3)*xx)*co+(-2+o*o*(2-3*xx)+xx)*so);
    }
    assemble_b();
}

/** Prints the stiffness matrix as a text array. */
void cubic_1d_fe::print_matrix() {
    int i,j;

    // Allocate zero vector, and workspace to compute a matrix product
    double *r=new double[2*dof],*s=r+dof;
    for(i=0;i<dof;i++) r[i]=0.;

    for(int i=0;i<dof;i++) {

        // Apply the black box matrix multiplication routine to a unit vector,
        // in order to extract a column of matrix entries
        r[i]=1.;mul_A(r,s);r[i]=0.;

        // Print a row of matrix entries. This assumes the matrix is symmetric
        // (as required for conjugate gradient) so that the row<->column switch
        // is permissible.
        for(j=0;j<dof-1;j++) printf("%g ",s[j]);
        printf("%g\n",s[dof-1]);
    }
    delete [] r;
}

/** Performs multiplication on a vector by the stiffness matrix.
 * \param[in] in the input vector.
 * \param[in] out the matrix-vector product. */
void cubic_1d_fe::mul_A(double *in,double *out) {
    int i,j,k;

    // Pre-computed integrals of derivatives of Lagrange polynomials, which are
    // required to construct the stiffness matrix
    const double B[16]={6/5.,1/10.,-6/5.,1/10.,
                        1/10.,2/15.,-1/10.,-1/30.,
                        -6/5.,-1/10.,6/5.,-1/10.,
                        1/10.,-1/30.,-1/10.,2/15.},
                 C[16]={3/5.,1/10.,-3/5.,0.,
                        1/10.,1/30.,-1/10.,-1/60.,
                        -3/5.,-1/10.,3/5.,0.,
                        0.,-1/60.,0.,1/10.};

    // Set the output vector to initially be zero
    for(i=0;i<dof;i++) out[i]=0.;

    // Loop over each interval, and compute the contribution from
    // each
    for(k=0;k<2*n;k+=2) for(i=(k==0?1:0);i<4;i++)
        for(j=(k==0?1:0);j<4;j++)
            out[-1+k+i]+=((0.5*k+1./h)*B[i+4*j]+C[i+4*j])*in[-1+k+j];
}

/** Computes the source vector in the linear system, which is based
 * on multiplying the source function by the mass matrix. */
void cubic_1d_fe::assemble_b() {
    int i,j,k;

    // Pre-computed integrals of Lagrange polynomial products
    const double D[16]={13/35.,11/210.,9/70.,-13/420.,
                        11/210.,1/105.,13/420.,-1/140.,
                        9/70.,13/420.,13/35.,-11/210.,
                        -13/420.,-1/140.,-11/210.,1/105.};

    // Clear the source function
    for(i=0;i<dof;i++) b[i]=0.;

    // Loop over each interval, and compute the contributions
    // from considering each Lagrange polynomial pair
    for(k=0;k<2*n;k+=2) for(i=(k==0?1:0);i<4;i++)
        for(j=0;j<4;j++) b[-1+k+i]-=D[i+4*j]*f[k+j];

    // Normalize the results, and add in the Neumann condition to the last
    // entry
    for(i=0;i<dof;i++) b[i]*=h;
    b[2*n-1]+=2*g;
}

/** Prints the solution.
 * \param[in] fp a file handle to write to. */
void cubic_1d_fe::print(FILE *fp) {
    double xx=1+h;
    fprintf(fp,"1 0 0 %g\n",*f);
    for(int i=1;i<=2*n;i+=2,xx+=h) fprintf(fp,"%g %g %g %g\n",xx,x[i],b[i],f[i+1]);
}

/** Prints the finite element solution function, sampled at equal intervals.
 * \param[in] fp a file handle to write to.
 * \param[in] I the number of intervals to use. */
void cubic_1d_fe::print_fe_function(FILE *fp,int I) {
    double spc=1./I,xx;
    for(int i=0;i<=I;i++) {
        xx=1+spc*i;
        fprintf(fp,"%g %g\n",xx,sol(xx));
    }
}

/** Prints the solution.
 * \param[in] filename the name of the file to write to. */
void cubic_1d_fe::print(const char* filename) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    print(fp);
    fclose(fp);
}

/** Prints the finite element solution function, sampled at equal intervals.
 * \param[in] filename the name of the file to write to. */
void cubic_1d_fe::print_fe_function(const char* filename,int I) {
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Error opening file\n",stderr);
        exit(1);
    }
    print_fe_function(fp,I);
    fclose(fp);
}

/** Computes the L2 norm between the numerical solution and the true
 * manufactured solution. The routine uses the trapezoid rule to evaluate the
 * integral.
 * \return The L2 norm. */
double cubic_1d_fe::l2_norm_mms() {
    double l2=0.,xx=1.;
    for(int i=1;i<2*n-1;i+=2) {
        xx+=h;
        l2+=mms_dsq(xx,x[i]);
    }

    // Add the contribution at the last point, including a factor of 1/2 due to
    // the trapezoid rule. Note that there is no contribution at x=1, since the
    // numerical solution is zero there, and hence matches the manufactured
    // solution perfectly.
    l2+=0.5*mms_dsq(2.,x[2*n-1]);
    return sqrt(h*l2);
}

/** Evaluates the solution at given x position in terms of a sum of finite element
 * basis functions.
 * \param[in] xhe position to consider. */
double cubic_1d_fe::sol(double xx) {
    xx=(1./h)*(xx-1);

    int k=static_cast<int>(xx);
    if(k<0) k=0;
    else if(k>n-1) k=n-1;

    xx-=k;
    double *p=x+2*k,al=k==0?0:p[-1];

    return al+xx*((p[1]-al)*xx*(3-2*xx)+(xx-1)*(*p*(xx-1)+p[2]*xx));
}
