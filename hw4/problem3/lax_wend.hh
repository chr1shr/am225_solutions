#ifndef DIFFUSE_HH
#define DIFFUSE_HH

#include <cstdio>
#include <cmath>

class lax_wend {
    public:
        /** The grid size. */
        const int m;
        /** The grid spacing. */
        const double dx;
        /** The primary grid for storing the solution. */
        double *a;
        /** The secondary grid for storing the solution. */
        double *b;
        /** The advection velocities at the interval centers. */
        double* const A_grid;
        /** The advection velocities at the interval boundaries. */
        double* const A_off;
        lax_wend(int m_);
        ~lax_wend();
        void step_forward(double dt);
        void init_step_function();
        void init_exp_sine();
        void init_triangle();
        void solve(double duration,double safe_fac);
        void solve(const char* filename,int snaps,double duration,double safe_fac);
        double integral();
    private:
        inline double timestep_select(double interval,double safe_fac,int &iters);
        void print_line(FILE *fp,double x,double *zp,int snaps);
        /** Calculates the diffusion constant at a position.
         * \param[in] x the position to consider.
         * \return The diffusion constant. */
        inline double A_fun(double x) {
            return 2+(4./3)*sin(x);
        }
        /** The maximum of the beta function. */
        inline double A_max() {return 10./3;}
};

#endif
