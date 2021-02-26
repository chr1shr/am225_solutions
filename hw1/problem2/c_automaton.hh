#ifndef C_AUTOMATON_HH
#define C_AUTOMATON_HH

#ifdef _OPENMP
#include "omp.h"
#endif

class c_automaton {
    public:
        const int m;
        const int n;
        const int mn;
        const unsigned int survive;
        const unsigned int birth;
        bool* r;
        bool* s;
        c_automaton(int m_,int n_,unsigned int survive_,unsigned int birth_);
        ~c_automaton();
        void time_steps(int steps);
        void clear();
        void set_random(int li,int ui,int lj,int uj,float prop);
        void step();
        void ascii_art();
        int count();
    private:
#ifdef _OPENMP
        inline double wtime() {return omp_get_wtime();}
#else
#include <ctime>
        inline double wtime() {return double(clock())*(1./CLOCKS_PER_SEC);}
#endif
};

#endif
