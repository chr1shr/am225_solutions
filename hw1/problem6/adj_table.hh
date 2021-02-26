#ifndef ADJ_TABLE_HH
#define ADJ_TABLE_HH

struct adj_table {
    int c[8];
    unsigned int* rec[9];
    adj_table();
    ~adj_table();
    void allocate();
    int total();
    bool same(adj_table* ap);
    inline void add(unsigned int key,int l) {
        rec[l][c[l]++]=key;
    }
};

#endif
