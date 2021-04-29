#include <cstdlib>
#include <cstring>

#include "lax_wend.hh"

/** Initializes the class for solving the diffusion equation on the periodic
 * unit interval, using a spatially dependent diffusion constant.
 * \param[in] m_ the number of gridpoints to use. */
lax_wend::lax_wend(int m_) : m(m_), dx(2*M_PI/m), a(new double[m]),
    b(new double[m]), A_grid(new double[m]), A_off(new double[m+1]) {

    // Initialize the A arrays
    for (int i=0;i<m;i++) {
        A_grid[i]=A_fun((i+0.5)*dx);
        A_off[i]=A_fun(i*dx);
    }
    A_off[m]=*A_off;
}

/** The class destructor frees the dynamically allocated memory. */
lax_wend::~lax_wend() {
    delete [] A_off;
    delete [] A_grid;
    delete [] b;
    delete [] a;
}

/** Performs one explicit timestep of the solution.
 * \param[in] dt the timestep to use. */
void lax_wend::step_forward(double dt) {
    double pre=dt/dx;
    for(int j=0;j<m;j++) {

        // Compute indices on left and right, taking into account periodicity
        int jl=j==0?m-1+j:(j-1),
            jr=j==m-1?1-m+j:(j+1);

        // Compute the fluxes
        double Fl=0.5*((A_grid[jl]*a[jl]+A_grid[j]*a[j])
                      -A_off[j]*pre*(A_grid[j]*a[j]-A_grid[jl]*a[jl])),
               Fr=0.5*((A_grid[j]*a[j]+A_grid[jr]*a[jr])
                      -A_off[j+1]*pre*(A_grid[jr]*a[jr]-A_grid[j]*a[j]));

        // Compute the update in the b array
        b[j]=a[j]-pre*(Fr-Fl);
    }

    // Swap pointers so that b becomes the primary array
    double *c=a;a=b;b=c;
}

/** Solves the advection equation with the spatially dependent velocity using the
 * generalized Lax-Wendroff method. This solving routine produces no output and is
 * used for convergence tests.
 * \param[in] snaps the number of snapshots to save (not including the initial
 *                  snapshot).
 * \param[in] duration the number of iterations to step the solution forward by
 *                     between snapshots.
 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction. */
void lax_wend::solve(double duration,double safe_fac) {
    int iters;
    double dt=timestep_select(duration,safe_fac,iters);
    for(int k=0;k<iters;k++) step_forward(dt);
}

/** Solves the diffusion equation, storing snapshots of the solution to a file.
 * \param[in] filename the name of the file to write to.
 * \param[in] snaps the number of snapshots to save (not including the initial snapshot).
 * \param[in] iters the number of iterations to step the solution forward by
 *                  between snapshots.
 * \param[in] type the integration type to use. 0: a finite-difference method,
 *                 1: a finite-volume method. */
void lax_wend::solve(const char* filename,int snaps,double duration,double safe_fac) {
    int iters;
    double dt=timestep_select(duration/snaps,safe_fac,iters);

    // Allocate memory to store solution snapshots. Integrate the system and
    // store snapshots after fixed time intervals
    double *z=new double[m*(snaps+1)];
    memcpy(z,a,m*sizeof(double));
    for(int i=1;i<=snaps;i++) {
        for(int k=0;k<iters;k++) step_forward(dt);

        // Store the snapshot
        memcpy(z+i*m,a,m*sizeof(double));
    }

    // Open the output file to store the snapshots
    FILE *fp=fopen(filename,"w");
    if(fp==NULL) {
        fputs("Can't open output file\n",stderr);
        exit(1);
    }

    // Print the snapshots, including periodic copies
    // at either end to get a full line over the interval
    // from 0 to 1
    print_line(fp,-0.5*dx,z+(m-1),snaps);
    for(int j=0;j<m;j++) print_line(fp,(j+0.5)*dx,z+j,snaps);
    print_line(fp,2*M_PI+0.5*dx,z,snaps);

    // Delete snapshots and close file
    fclose(fp);
    delete [] z;
}

/** Initializes the solution to be a step function. */
void lax_wend::init_step_function() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=x>0.25&&x<0.75?1:0;
    }
}

/** Initializes the solution to be the exponential of a combination of sine
 * waves. */
void lax_wend::init_exp_sine() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5);
        a[i]=exp(sin(x)+0.5*sin(4*x));
    }
}

/** Initializes the solution to be the triangle hat function in part (e). */
void lax_wend::init_triangle() {
    for(int i=0;i<m;i++) {
        double x=dx*(i+0.5),v=1+0.5*M_PI-fabs(x-M_PI);
        a[i]=v<0?0:v;
    }
}

/** Computes the integral of the solution over the domain.
 * \return The integral. */
double lax_wend::integral() {
    double sum=0;
    for(double *ap=a;ap<a+m;ap++) sum+=*ap;
    return dx*sum;
}

/** Selects the timestep to use, based on the largest acceptable value where a multiple
 * fits into the given time interval.
 * \param[in] interval the time interval for the timesteps to exactly fit into.
 * \param[in] safe_fac a safety factor to apply to the CFL timestep restriction.
 * \param[out] iters the number of timesteps required.
 * \return The timestep. */
inline double lax_wend::timestep_select(double interval,double safe_fac,int &iters) {
    double dt=dx/A_max()*safe_fac;
    iters=static_cast<int>(interval/dt)+1;
    return interval/iters;
}

/** Prints a line of stored snapshots to a file.
 * \param[in] fp a pointer to the file to write to.
 * \param[in] x the position in the domain corresponding to this line.
 * \param[in] zp a pointer to the first snapshot data point to print.
 * \param[in] snaps the number of snapshots (not including the starting
 *                  snapshot). */
void lax_wend::print_line(FILE *fp,double x,double *zp,int snaps) {
    fprintf(fp,"%g",x);
    for(int i=0;i<=snaps;i++) fprintf(fp," %g",zp[i*m]);
    fputc('\n',fp);
}
