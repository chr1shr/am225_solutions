#include <cstdio>
#include <cstdlib>

#include "c_automaton.hh"

int main() {
    c_automaton ca(80,40,0b111110,0b1000);
//    c_automaton ca(1000,1000,0b111100000,0b111101000);
//    ca.clear();
srand(10);
    ca.set_random(34,46,14,26,0.75);
    ca.ascii_art();

	printf("generation %d:\n",0);
	for(int i=0;i<80;i++)putchar('-');
	putchar('\n');
	ca.ascii_art();

    for(int q=0;q<6;q++) {
        for(int p=0;p<25;p++) ca.step();
        putchar('\n');
	
		printf("generation %d:\n",q+1);
	for(int i=0;i<80;i++)putchar('-');
	putchar('\n');
	ca.ascii_art();
    }
}
