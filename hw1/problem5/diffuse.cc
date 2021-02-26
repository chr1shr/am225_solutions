#include "diffuse.hh"

#include <cstdio>
#include <cmath>

/** The class constructor sets constants and dynamically allocate memory
 * for the solution vector.
 * \param[in] N_ the number of gridpoints in one dimension of the grid.
 * \param[in] (al_,be_,ga_) the three parameters in the 2D diffusion equation.
 */
diffuse::diffuse(int N_,double al_,double be_,double ga_) : N(N_),
    NN(N*N), al(al_), be(be_), ga(ga_), dx(1./N), xsp(double(N_)),
    u(new double[NN]), v(new double[NN]) {}

/* The class destructor frees the dynamically allocated memory. */
diffuse::~diffuse() {
    delete [] v;
    delete [] u;
}

/** Initializes the solution to be a periodic Gaussian-like function. */
void diffuse::init() {
    for(int j=0;j<N;j++) {
        double y=j*dx;
        for(int i=0;i<N;i++) {
            double x=i*dx;
            u[i+N*j]=exp(-cos(2*M_PI*x)-cos(2*M_PI*y));
        }
    }
}

/** Computes the standard deviation of the current solution. */
double diffuse::std_deviation() {
    const double iNN=1./NN;
    double su=0,suu=0,de;

    // Compute the mean of the solution field
#pragma omp parallel for reduction(+:su)
    for(int i=0;i<NN;i++) su+=u[i];
    su*=iNN;

    // Compute the variance of the solution field
#pragma omp parallel for reduction(+:suu)
    for(int i=0;i<NN;i++) {de=u[i]-su;suu+=de*de;}

    // Normalize and take the square root to obtain the standard deviation
    return sqrt(suu*iNN);
}

/** Integrates the solution forward.
 * \param[in] steps the number of timesteps to take. If this is set to zero,
 *                  then the code steps forward by a duration of 0.1.
 * \param[in] std_dev whether to output status messages about the standard
 *                    deviation of the solution.
 * \return The number of steps taken. */
int diffuse::integrate(int steps,bool std_dev) {
    const double duration=0.1;

    // Set timestep
    double adt=0.08*dx*dx/al,dt;
    if(steps==0) {
        steps=static_cast<int>(duration/adt)+1;dt=duration/steps;
    } else dt=adt;

    // Perform timesteps
    for(int k=0;k<steps;k++) {
        if(std_dev&&k%1000==0) printf("%g %g\n",dt*k,std_deviation());
        step_forward(dt);
    }
    if(std_dev) printf("%g %g\n",dt*steps,std_deviation());
    return steps;
}

/** Steps the solution forward in time according to the 2D diffusion equation.
 * \param[in] dt the timestep to take. */
void diffuse::step_forward(double dt) {
    const double f1=0.25*be*xsp*xsp*dt,
                 f2=al*xsp*xsp*dt,
                 f3=ga*xsp*xsp*dt,
                 f4=1-2*f2-2*f3;

    // Perform the finite difference update, computing the new solution in the
    // v array
#pragma omp parallel for num_threads(num_t)
    for(int j=0;j<N;j++) {
        double *uc=u+j*N,*um=uc+(j==0?NN-N:-N),*up=uc+(j==N-1?N-NN:N),
               *vc=v+j*N;
        for(int i=0;i<N;i++) {
            int im=i==0?N-1:i-1,ip=i==N-1?0:i+1;
            vc[i]=(um[im]-um[ip]-up[im]+up[ip])*f1+(uc[im]+uc[ip])*f2
                 +(um[i]+up[i])*f3+uc[i]*f4;
        }
    }

    // Swap the pointers to u and v, so that the new solution is now in u
    double *w=v;v=u;u=w;
}
