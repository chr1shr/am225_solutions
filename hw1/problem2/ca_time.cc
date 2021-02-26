#include <cstdio>
#include <cstdlib>

#include "omp.h"
#include "c_automaton.hh"
#include "stopwatch.hh"
#include <vector>

int main() {

	static const int n_nts = 4;
	static const int n_gss = 7;
	static double duration=2;

	int gs0=16;

	stopwatch sw;

	for (int nt=1;nt<=n_nts;nt++) {
		omp_set_num_threads(nt);
		for (int g=0,gs=gs0;g<n_gss;g++,gs<<=1) {

			c_automaton ca(gs,gs,0b111110,0b1000);

			int iters=0;

			sw.reset();
			
			do {

				ca.clear();
				ca.set_random(gs/2-6,gs/2+6,gs/2-6,gs/2+6,0.75);

				sw.tic();
				ca.step();
				sw.toc();

				iters++;
			} while (sw.elapsed() < duration);

			printf("%12.6g ",sw.elapsed()/iters);
			fflush(stdout);
		}
		putchar('\n');
	}
}
