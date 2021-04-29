#include "squircle.hh"
#include "file_output.hh"

#include <cstdio>
#include <cmath>

/** Initializes the class for solving a Poisson problem on a circle using a
 * mapping from a square.
 * \param[in] n_ the number of squares to divide the domain into.
 * \param[in] q_ the number of quadrature points to use for the finite element
 *               calculations. */
squircle::squircle(int n_,int q_) : conj_grad((n_-1)*(n_-1)), n(n_), m(n-1),
    q(q_), h(2./n), hh(0.5*h), gq(q), a(new double[9*dof]),
    phi_base(new double[2*q]), phi(phi_base+q) {
    int i,j,ci,cj,di,dj,li,lj,ui,uj,qi,qj;
    double phix,phiy,phiu,phiv,phix2,phiy2,phiu2,phiv2,
           d11,d12,d21,d22,det,xx,yy;

    // Compute tables of the finite element functions
    for(i=0;i<q;i++) phi_base[i]=0.5*(1-gq.x[i]);
    for(i=0;i<q;i++) phi[i]=0.5*(1+gq.x[i]);

    // Scale the quadrature weights for integrating on squares of side length h
    gq.scale(hh);

    // Clear the solution array
    for(i=0;i<dof;i++) x[i]=0.;

    // Clear stiffness matrix entries
    for(i=0;i<9*dof;i++) a[i]=0;

    // Assemble the stiffness matrix. Loop over the squares in the grid.
    for(j=0;j<n;j++) {
        lj=j==0?0:-1;uj=j==n-1?0:1;
        yy=-1+hh+j*h;
        for(i=0;i<n;i++) {
            li=i==0?0:-1;ui=i==n-1?0:1;
            xx=-1+hh+i*h;

            // Loop over the quadrature points in the square
            for(qj=0;qj<q;qj++) for(qi=0;qi<q;qi++) {

                // Compute the Jacobian integration factors at the
                // quadrature point
                jacobian(xx+gq.x[qi],yy+gq.x[qj],d11,d12,d21,d22,det);

                // Loop over all interior corners of the square, to
                // consider all basis function patches
                for(cj=lj;cj<uj;cj++) for(ci=li;ci<ui;ci++) {
                    phix=phi[cj*q+qj]*(ci==-1?-1:1)/h;
                    phiy=phi[ci*q+qi]*(cj==-1?-1:1)/h;
                    phiu=d11*phix+d12*phiy;
                    phiv=d21*phix+d22*phiy;

                    // Loop again over all interior corners of the square,
                    // to consider the integral of any pair of basis
                    // function patches
                    for(dj=lj;dj<uj;dj++) for(di=li;di<ui;di++) {
                        phix2=phi[dj*q+qj]*(di==-1?-1:1)/h;
                        phiy2=phi[di*q+qi]*(dj==-1?-1:1)/h;
                        phiu2=d11*phix2+d12*phiy2;
                        phiv2=d21*phix2+d22*phiy2;
                        a[9*(m*(j+cj)+(i+ci))+3*(dj-cj)+(di-ci)+4]
                            +=gq.w[qj]*gq.w[qi]*(phiu*phiu2+phiv*phiv2)*det;
                    }
                }
            }
        }
    }
}

/** The class destructor frees the dynamically allocated memory. */
squircle::~squircle() {
    delete [] phi_base;
    delete [] a;
}

/** Performs multiplication on a vector by the stiffness matrix.
 * \param[in] in the input vector.
 * \param[in] out the matrix-vector product. */
void squircle::mul_A(double *in,double *out) {
    int i,j,ci,cj,li,lj,ui,uj;
    double *ap=a+4;

    // Loop over the interior nodes of the grid, and find the range of entries
    // of the stiffness matrix to use
    for(j=0;j<m;j++) {
        lj=j==0?0:-1;uj=j==m-1?1:2;
        for(i=0;i<m;i++,in++,out++,ap+=9) {
            li=i==0?0:-1;ui=i==m-1?1:2;

            // Perform the multiplication to evaluate one term of the
            // matrix-vector product
            *out=0;
            for(cj=lj;cj<uj;cj++) for(ci=li;ci<ui;ci++)
                *out+=ap[3*cj+ci]*in[m*cj+ci];
        }
    }
}

/** Prints the stiffness matrix as a text array. */
void squircle::print_matrix() {
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

/** Prints diagnostic information about the matrix coefficients. */
void squircle::diagnostic_a() {
    double *ap=a;
    int i,j;
    for(j=0;j<m;j++) for(i=0;i<m;i++,ap+=9)
        printf("%d %d %g %g %g %g %g %g %g %g %g\n",i,j,*ap,ap[1],
               ap[2],ap[3],ap[4],ap[5],ap[6],ap[7],ap[8]);
}

/** Sets the source term (right-hand side) for the finite element calculation.
 * \param[in] type the type of exact solution. */
void squircle::set_source(int type) {
    int i,j,ci,cj,li,lj,ui,uj,qi,qj;
    double v,w,det,xx,yy;
    for(j=0;j<n;j++) {
        lj=j==0?0:-1;uj=j==n-1?0:1;
        yy=-1+hh+j*h;
        for(i=0;i<n;i++) {
            li=i==0?0:-1;ui=i==n-1?0:1;
            xx=-1+hh+i*h;
            for(qj=0;qj<q;qj++) for(qi=0;qi<q;qi++) {
                mapping(xx+gq.x[qi],yy+gq.x[qj],v,w,det);
                for(cj=lj;cj<uj;cj++) for(ci=li;ci<ui;ci++)
                    b[m*(j+cj)+(i+ci)]-=phi[cj*q+qj]*gq.w[qj]*phi[ci*q+qi]*gq.w[qi]
                        *det*(type==0?1:-exp(-v)*(3+(v-4)*v+w*w));
            }
        }
    }
}

/** Computes the L2 norm between the numerical solution and an exact solution.
 * \param[in] type the type of exact solution. */
double squircle::l2_norm_mms(int type) {
    int i,j;
    double v,w,xx,yy,det,sum=0.;
    for(j=0;j<m;j++) {
        yy=-1+(j+1)*h;
        for(i=0;i<m;i++) {
            xx=-1+(i+1)*h;

            // Perform trapezoid rule calculation, applying determinant factor
            // due to change of variable from circle to square
            mapping(xx,yy,v,w,det);
            sum+=det*(x[j*m+i]-(v*v+w*w-1)*(type==0?0.25:-exp(-v)));
        }
    }
    return h*h*sum;
}

/** Computes the inverse and determinant of the Jacobian of the
 * square-to-circle mapping, for use in the integration.
 * \param[in] (x,y) the position in the square.
 * \param[out] (d11,d12,d21,d22) the components of the inverse Jacobian.
 * \param[out] det the determinant of the Jacobian. */
inline void squircle::jacobian(double x,double y,double &d11,double &d12,double &d21,double &d22,double &det) {
    double fac=2-x*x-y*y,ifac=1./fac,sx=sqrt(1-0.5*x*x),sy=sqrt(1-0.5*y*y);
    d11=ifac*(2-x*x)*sy;
    d12=ifac*x*y*sy;
    d21=ifac*x*y*sx;
    d22=ifac*(2-y*y)*sx;
    det=0.5*fac/(sx*sy);
}

/** Computes the mapped variables and the determinant of the mapping.
 * \param[in] (x,y) the position in the square.
 * \param[out] (v,w) the position in the circle.
 * \param[out] det the determinant of the mapping. */
inline void squircle::mapping(double x,double y,double &v,double &w,double &det) {
    double sx=sqrt(1-0.5*x*x),sy=sqrt(1-0.5*y*y);
    v=x*sy;w=y*sx;
    det=(1-0.5*(x*x+y*y))/(sx*sy);
}

/** Saves the solution on the square in the Gnuplot matrix binary format.
 * \param[in] filename the name of the file to write to. */
void squircle::output(const char *filename) {
    int i,j;

    // Open file and print error if there is a problem
    FILE *outf=fopen(filename,"wb");
    if(outf==NULL) {
        fputs("Can't open file\n",stderr);
        exit(1);
    }

    // Allocate memory and write the header file
    float *fbuf=new float[n+2],*fp=fbuf;
    double *pp;
    *(fp++)=n+1;for(i=0;i<=n;i++) *(fp++)=-1+i*h;
    sfwrite(fbuf,sizeof(float),n+2,outf);

    // Write field entries line-by-line
    fbuf[1]=fbuf[n+1]=0;
    for(j=0;j<=n;j++) {

        // Write header entry
        *fbuf=-1+j*h;fp=fbuf+2;

        // Write a horizontal line to the buffer
        if(j==0||j==n) {while(fp<fbuf+n+1) *(fp++)=0;}
        else {
            pp=x+(j-1)*m;
            for(i=0;i<m;i++) *(fp++)=static_cast<float>(*(pp++));
        }
        sfwrite(fbuf,sizeof(float),n+2,outf);
    }

    // Remove temporary memory and close file
    delete [] fbuf;
    fclose(outf);
}
