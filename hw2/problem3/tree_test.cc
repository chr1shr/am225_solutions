#include "tree.hh"
#include <cstdio>

int main() {

	puts("allocation");
	tree t;
	printf("count %d should be 1\n\n",t.count());

	puts("failed insert");
	printf("%s should be false\n",t.insert(1)?"true":"false");
	printf("count: %d should be 1\n\n",t.count());

	puts("successful adds");
	printf("%s should be true\n",t.insert(0)?"true":"false");
	printf("%s should be true\n",t.insert(1)?"true":"false");
	printf("%s should be true\n",t.insert(0)?"true":"false");
	printf("%s should be true\n\n",t.insert(2)?"true":"false");
	printf("count: %d should be 5\n\n",t.count());
	
	puts("failed insert");
	printf("%s should be false\n",t.insert(5)?"true":"false");
	printf("count: %d should be 5\n\n",t.count());

	unsigned cs[5];
	t.list(cs);
	printf("counts [ %d %d %d %d %d ] should be [ 2 1 1 0 0 ]\n",
			cs[0],cs[1],cs[2],cs[3],cs[4]);

	unsigned counts[5] = {2,1,1,0,0};
	unsigned wrong[5] = {2,2,0,0,0};
	printf("%s should be true\n",t.is_list(counts)?"true":"false");
	printf("%s should be false\n\n",t.is_list(wrong)?"true":"false");

	puts("add one at end so that leaf needs to be moved");
	printf("%s should be true\n",t.insert(0)?"true":"false");
	printf("%s should be true\n",t.insert(5)?"true":"false");
	unsigned cnew[7] = {3,1,1,0,1,0,0};
	printf("%s should be true\n",t.is_list(cnew)?"true":"false");
	t.list(cnew);
	printf("counts [ %d %d %d %d %d %d %d ] should be [ 3 1 1 0 1 0 0 ]\n\n",
			cnew[0],cnew[1],cnew[2],cnew[3],cnew[4],cnew[5],cnew[6]);

	puts("try copying and adding to one");
	tree t2(&t);
	t.insert(0);
	unsigned c2[8];

	t2.list(c2);
	printf("counts [ %d %d %d %d %d %d %d ] should be [ 3 1 1 0 1 0 0 ]\n",
			c2[0],c2[1],c2[2],c2[3],c2[4],c2[5],c2[6]);
	t.list(c2);
	printf("counts [ %d %d %d %d %d %d %d %d ] should be [ 4 1 1 0 1 0 0 0 ]\n\n",
			c2[0],c2[1],c2[2],c2[3],c2[4],c2[5],c2[6],c2[7]);

	puts("make a new tree from one of these lists");
	tree t3(c2);
	t3.list(c2);
	printf("counts [ %d %d %d %d %d %d %d %d ] should be [ 4 1 1 0 1 0 0 0 ]\n\n",
			c2[0],c2[1],c2[2],c2[3],c2[4],c2[5],c2[6],c2[7]);

}
