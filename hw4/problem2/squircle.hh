#ifndef SQUIRCLE_HH
#define SQUIRCLE_HH

#include "quadrat.hh"
#include "conj_grad.hh"

class squircle : public conj_grad {
    public:
        /** The number of grid squares along one side of the [-1,1]^2 domain.
         */
        const int n;
        /** The number of interior nodes on the grid in one direction. */
        const int m;
        /** The number of quadrature points for evaluating the stiffness
         * matrix. */
        const int q;
        /** The grid spacing. */
        const double h;
        /** Half the grid spacing. */
        const double hh;
        squircle(int n_,int q_);
        ~squircle();
        void print_matrix();
        void diagnostic_a();
        void set_source(int type);
        double l2_norm_mms(int type);
        void output(const char *filename);
        virtual void mul_A(double *in,double *out);
    private:
        inline void jacobian(double x,double y,double &d11,double &d12,double &d21,double &d22,double &det);
        inline void mapping(double x,double y,double &v,double &w,double &det);
        /** The quadrature weights. */
        quadrat gq;
        /** The stiffness matrix terms. */
        double* const a;
        /** The table to the phi values at quadrature points. */
        double* const phi_base;
        /** A pointer to an element in the phi table, used for easy
         * referencing. */
        double* const phi;
};

#endif
