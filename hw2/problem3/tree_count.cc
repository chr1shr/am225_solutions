#include "tree.hh"
#include <cstdio>

int main() {

	const int max_order=15;

	tree t;
	std::vector<unsigned*> pattern_sets[max_order+1];
	pattern_sets[1].push_back(new unsigned[1]);
	*pattern_sets[1].back() = 0;

	puts("order  1:     1 trees");

	for (int order=2;order<=max_order;order++) {


		unsigned *curr = new unsigned[order];
		for(unsigned j=0;j<pattern_sets[order-1].size();j++) {

			for (int node=0;node<order-1;node++) {

				// copy tree
				tree t_new(pattern_sets[order-1][j]);
//				printf("o=%d, j=%d, node=%d | stem: ",
//						order,j,node); t_new.print_list();

				// add at node
				t_new.insert(node,pattern_sets);
//				printf("o=%d, j=%d, node=%d | add : ",
//						order,j,node); t_new.print_list();
				t_new.list(curr);

				bool seen = false;
				for(unsigned i=0;i<pattern_sets[order].size()&&!seen;i++) {
					bool seen_p=true;
					for(int j=0;j<order&&seen_p;j++) {
						if(curr[j]!=pattern_sets[order][i][j])seen_p=false;
					}
					if(seen_p) seen = true;
				}

				if (!seen) {
					pattern_sets[order].push_back(curr);
					curr = new unsigned[order];
				}

			}
		}
		delete[] curr;

		printf("order %2d: %5d trees\n",order,
				static_cast<int>(pattern_sets[order].size()));

		char fn[256];
		sprintf(fn,"trees%d",order);
		FILE *fh = fopen(fn,"w");
		//*
		for(unsigned i=0;i<pattern_sets[order].size();i++) {
			for(int j=0;j<order;j++)
				fprintf(fh,"%2d ",static_cast<int>(pattern_sets[order][i][j]));
			fprintf(fh,"\n");
		}
		fclose(fh);
		// */
	}
}


