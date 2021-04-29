#include <cstdio>
#include <cstdlib>
#include <cmath>

#include <cstdio>
#include <cmath>
#include <algorithm>

#include "rbf_hilbert.hh"
#include "omp.h"

int main() {

	int max = 25000;
	for(int k=10;k<=max;k+=std::max(1,std::min(k>>2,max-k))) {
		fprintf(stderr,"%d\n",k);
		rbf r(k,2);
		r.set_length_scale(5./sqrt(k));
		r.init_random();
		r.assemble_matrix();

		double t0=omp_get_wtime(),t1,t_cg=0,tc1;
		int l=0;
		/*
		do {
			r.solve_weights_lapack();
			l++;t1=omp_get_wtime();
		} while(t1<t0+0.5);
		t_lapack=(t1-t0)/l;
		fprintf(stderr,"lapack done...");

		r.make_table();
		t0=omp_get_wtime();l=0;
		*/
		int bls=sqrt(k);
		do {
			r.solve_weights_conj_grad(bls);
			l++;t1=omp_get_wtime();
		} while(t1<t0+1);
		t_cg=(t1-t0)/l;
		tc1=t_cg;

		rbf_hilbert rh(k,2);
		rh.set_length_scale(5./sqrt(k));
		rh.init_hilbert();
		rh.assemble_matrix();

		t0=omp_get_wtime();t_cg=0;
		l=0;
		/*
		do {
			rh.solve_weights_lapack();
			l++;t1=omp_get_wtime();
		} while(t1<t0+0.5);
		t_lapack=(t1-t0)/l;
		fprintf(stderr,"lapack done...");

		rh.make_table();
		t0=omp_get_wtime();l=0;
		bls=sqrt(k);
		*/
		do {
			rh.solve_weights_conj_grad(bls);
			l++;t1=omp_get_wtime();
		} while(t1<t0+1);
		t_cg=(t1-t0)/l;

		printf("%d %g %g\n",k,tc1,t_cg);
	}
}
