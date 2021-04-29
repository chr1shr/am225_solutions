#include <cstdio>
#include "poisson_schur.hh"

poisson_schur::poisson_schur(int n_, const char* filename) :
    n(n_-1), nn(n*n), h(1./n_), ih2(n_*n_), f(new double[nn]), c(new int[nn])
{
    for (int i=0; i<nn; ++i) c[i]=0;
    init_subdomains(filename);
    this->conj_grad::init(ng);
}

poisson_schur::~poisson_schur()
{
    for (int i=subdomains.size()-1; i>=0; --i) delete subdomains[i];
    delete [] c;
    delete [] f;
}

/** Solve Poisson's equation with the given source term using the Schur
 *  complement method. */
void poisson_schur::solve(const std::function<double(double,double)>& fun)
{
    // Evaluate the source term at the gridpoints
    for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
            f[n*i+j]=fun((j+1)*h,(i+1)*h);
        }
    }
    solve();
}

/** Solve Poisson's equation using the Schur complement method. */
void poisson_schur::solve()
{
    // First, solve the Schur complement system for the unknowns on the glue
    //
    //    S*ug = fg - Ag1*A11^-1*f1 - Ag2*A22^-1*f2 - ... - AgK*AKK^-1*fK
    //
    // where S = Agg - Ag1*A11^-1*A1g - Ag2*A22^-1*A2g - ... - AgK*AKK^-1*AKg.
    // Here Aii^-1 means "call the Poisson FFT solver on subdomain i."

    // Compute b = fg - Ag1*(A11^-1*f1) - ... - AgK*(AKK^-1*fK)
    for (const auto& g : g2l) b[g.second] = f[g.first];
    for (subdomain* dp : subdomains) {
        subdomain& dom = *dp;
        for (int i=0; i<dom.n; ++i) {
            for (int j=0; j<dom.n; ++j) {
                dom.v[dom.n*i+j] = f[n*(dom.y+i)+(dom.x+j)];
            }
        }
        dom.solve();
        for (const auto& g : dom.glue) b[g2l[g.first]] += dom.v[g.second]*ih2;
        for (int k=0; k<dom.nn; ++k) dom.v[k]=0;
    }

    // Solve for ug using conjugate gradient
    this->conj_grad::solve();

    // Now use ug to solve on the subdomains
    for (subdomain* dp : subdomains) {
        subdomain& dom = *dp;
        for (int i=0; i<dom.n; ++i) {
            for (int j=0; j<dom.n; ++j) {
                dom.v[dom.n*i+j] = f[n*(dom.y+i)+(dom.x+j)];
            }
        }
        for (const auto& g : dom.glue) dom.v[g.second] += x[g2l[g.first]]*ih2;
        dom.solve();
    }
}

void poisson_schur::mul_A(double* in, double* out)
{
    // Compute the glue-to-glue contributions
    for (const auto& g : g2l) {
        int global = g.first;
        int local  = g.second;
        out[local] = 4*in[local]*ih2;
        int i=global/n, j=global%n;
        if (j>0   && c[n*i+(j-1)]==0) out[local] -= in[g2l[n*i+(j-1)]]*ih2; // Left
        if (j<n-1 && c[n*i+(j+1)]==0) out[local] -= in[g2l[n*i+(j+1)]]*ih2; // Right
        if (i>0   && c[n*(i-1)+j]==0) out[local] -= in[g2l[n*(i-1)+j]]*ih2; // Bottom
        if (i<n-1 && c[n*(i+1)+j]==0) out[local] -= in[g2l[n*(i+1)+j]]*ih2; // Top
    }

    // Use "in" as glue to solve on the subdomains
    for (subdomain* dp : subdomains) {
        subdomain& dom = *dp;
        for (int k=0; k<dom.nn; ++k) dom.v[k]=0;
        for (const auto& g : dom.glue) dom.v[g.second] -= in[g2l[g.first]]*ih2;
        dom.solve();
        for (const auto& g : dom.glue) out[g2l[g.first]] += dom.v[g.second]*ih2;
    }
}

/** Initialize the subdomains from a file. */
void poisson_schur::init_subdomains(const char* filename)
{
    FILE* file = fopen(filename, "r");
    if (file == nullptr) {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }

    // Loop over the squares in the file
    int x, y, s, q, id=1;
    while ((q=fscanf(file,"%d %d %d",&x,&y,&s)) == 3) {
        subdomains.push_back(new subdomain(this,x,y,s));
        for(int i=y; i<y+s-1; ++i) for(int j=x; j<x+s-1; ++j) c[n*i+j]=id;
        id++;
    }

    if (q != EOF) {
        fprintf(stderr, "Error reading file %s\n", filename);
        exit(1);
    }
    fclose(file);

    // Count the number of glue points and build the global-to-local table
    ng = 0;
    for (int i=0; i<nn; ++i) if (c[i]==0) g2l[i]=ng++;
}

/** Get the solution at position (i,j).
 * \param[in] glue : Flag to evaluate just the glue. */
double poisson_schur::v(int i, int j, bool glue)
{
    int k = c[n*i+j];
    if (k==0) {
        return x[g2l[n*i+j]];
    } else if (!glue) {
        subdomain* dom = subdomains[k-1];
        return dom->v[dom->n*(i-dom->y)+(j-dom->x)];
    } else {
        return 0;
    }
}

/** Print the solution, padding with zeros for the Dirichlet boundary conditions.
 * \param[in] glue : Flag to print just the glue. */
void poisson_schur::print_solution(bool glue)
{
    for (int j=0; j<n+2; ++j) printf("0 ");
    puts("");
    for (int i=0; i<n; ++i) {
        printf("0 ");
        for (int j=0; j<n; ++j) printf("%g ", v(i,j,glue));
        puts("0");
    }
    for (int j=0; j<n+2; ++j) printf("0 ");
    puts("");
}
