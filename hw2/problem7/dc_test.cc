#include <cstdio>

#include "discont.hh"

int main() {

	sol_adapt* sa = new discont_fsal(true);
	sa->solve(48+exp(-1),1e-16,0,true,1000);

}
