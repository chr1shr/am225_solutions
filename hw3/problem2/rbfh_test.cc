#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "rbf_hilbert.hh"
#include <algorithm>

int main() {

	for (int type=0; type<2; type++) {
		for (int k=10; k<=10000; k+=std::max(1,std::min(k>>2,10000-k))) {
			rbf_hilbert r(k,2);
			r.set_length_scale(5./sqrt(k));

			if (type==0) {
				r.init_random();
			} else {
				r.init_hilbert();
		    }

			int tot,pre;
			r.matrix_analysis(sqrt(k),tot,pre);
			printf("%d %d %d %.2f%%\n",k,tot,pre,float(100.*pre/tot));
		}
		printf("\n\n");
	}

    //r.solve_weights_lapack();
//    r.output_points("rbf.pts");
    //r.output_interpolant("rbf.fld",101,1.1);
}
