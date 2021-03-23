#ifndef TREE_HH
#define TREE_HH

#include <vector>
#include <cstddef>

class tree {

	public:
	tree() {};
	tree(unsigned *c);
	tree(tree *t);
	~tree();

	bool is_list(unsigned *c) {return equals(c);};
	void print_list();
	void list(unsigned *c) {get_counts(c);};
	bool insert(int node,std::vector<unsigned*> *patterns=NULL) {return add(node,patterns)<1;};
	inline int count() {return count(0);};

	private:
	std::vector<tree*> kids;

	int pattern(std::vector<unsigned*> *patterns);
	bool equals(unsigned *(&c));
	void get_counts(unsigned *(&c));
	void sort_kids(std::vector<unsigned*> *patterns,int start);
	void add();
	int add(int node,std::vector<unsigned*> *patterns);
	int count(int start);

};

#endif
