#include "tree.hh"
#include <cstring>
#include <cstdio>
#include <cstdlib>

tree::tree(unsigned *c) {

	std::vector<unsigned> nodes;
	std::vector<unsigned> label;

	nodes.push_back(*c++);
	label.push_back(0);

	int val,target=0;
	while(nodes.size()>0 && nodes.front()>0) {
		val = (*c++);	
		insert(target);
		if (val>0) {

			nodes.push_back(val);
			label.push_back(this->count()-1);
			target=label.back();

		} else {

			while (nodes.size()>0 && (--nodes.back())==0) {
				nodes.pop_back();
				label.pop_back();
			}

			if (nodes.size()>0) target=label.back();
		}
	}
}

tree::tree(tree *t) {
	for(unsigned i=0;i<t->kids.size();i++)
		kids.push_back(new tree(t->kids[i]));
}

tree::~tree() {
	for(unsigned i=0;i<kids.size();i++) delete kids[i];
}

int tree::pattern(std::vector<unsigned *> *patterns) {
	if (patterns==NULL) return 0;
	int o = count();
	for (unsigned j=0;j<patterns[o].size();j++)
		if (is_list(patterns[o][j])) return j;
	fputs("error: pattern not found\n",stderr);
	exit(-1);
}

void tree::sort_kids(std::vector<unsigned*> *patterns,int start) {

	// first see if we have to move to the left
	bool good=false;
	unsigned i=static_cast<unsigned>(start);
	while (i>0 && !good) {

		if (kids[i]->count()<kids[i-1]->count()) {
			good=true;
		} else if (kids[i]->count()==kids[i-1]->count() &&
			kids[i]->pattern(patterns) < kids[i-1]->pattern(patterns)) {
			good=true;
		}

		if (!good) {
			tree *tmp=kids[i];
			kids[i]=kids[i-1];
			kids[--i]=tmp;
		}
	}

	// now see if we have to move to the right
	good=false;
	while (i<kids.size()-1 && !good) {

		if (kids[i]->count()>kids[i+1]->count()) {
			good=true;
		} else if (kids[i]->count()==kids[i+1]->count() &&
			kids[i]->pattern(patterns) > kids[i+1]->pattern(patterns)) {
			good=true;
		}

		if (!good) {
			tree *tmp=kids[i];
			kids[i]=kids[i+1];
			kids[++i]=tmp;
		}
	}
}

void tree::add() {
	kids.push_back(new tree());
}

int tree::add(int node,std::vector<unsigned*> *patterns) {
	if (node==0) {
		add();
		return -1;
	} else {
		for (unsigned i=0;i<kids.size();i++) {
			node = kids[i]->add(node-1,patterns);
			if (node<0) {
				sort_kids(patterns,i);
				return -1;
			}
		}
	}
	return node;
}

int tree::count(int start) {
	for(unsigned i=0;i<kids.size();i++)
		start = kids[i]->count(start);
	return start+1;
}

void tree::get_counts(unsigned *(&c)) {
	(*c++)=kids.size();
	for(unsigned i=0;i<kids.size();i++)kids[i]->get_counts(c);
}

bool tree::equals(unsigned *(&c)) {
	bool same = (*c++)==kids.size();
	for(unsigned i=0;i<kids.size();i++) {
		same = same && kids[i]->equals(c);
		// maybe short-circuit
	}
	return same;
}

void tree::print_list() {
	int n = count();
	unsigned *c = new unsigned[n];
	list(c);
	printf("[ ");
	for(int i=0;i<n;i++)
		printf("%2d ",static_cast<int>(c[i]));
	printf("]\n");
	delete[] c;
}
