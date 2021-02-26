#ifndef DIFFUSE_HH
#define DIFFUSE_HH



class diffuse {
    public:
        /** The number of gridpoints in one dimension of the grid. */
        const int N;
        /** The total number of gridpoints. */
        const int NN;
        /** The alpha parameter in the 2D diffusion equation to solve. */
        const double al;
        /** The beta parameter in the 2D diffusion equation to solve. */
        const double be;
        /** The gamma parameter in the 2D diffusion equation to solve. */
        const double ga;
        /** The grid spacing. */
        const double dx;
        /** The inverse grid spacing. */
        const double xsp;
        /** The solution vector. */
        double* u;
        /** A copy of the solution vector for stepping forward. */
        double* v;
        /** The number of threads to use. */
        int num_t;
        diffuse(int N_,double al_,double be_,double ga_);
        ~diffuse();
        void step_forward(double dt);
        void init();
        double std_deviation();
        int integrate(int steps,bool std_dev);
};

#endif
