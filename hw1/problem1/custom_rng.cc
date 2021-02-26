#include "custom_rng.hh"

custom_rng::custom_rng(unsigned long seed) : b(4101842887655102017L), c(1) {
    a=seed^b;int64();
    b=a;int64();
    c=b;int64();
}

unsigned long custom_rng::int64() {
    a=a*2862933555777941757L+7046029254386353087L;
    b^=b>>17;b^=b<<31;b^=b>>8;
    c=4294957665U*(c&0xffffffff)+(c>>32);
    unsigned long d=a^(a<<21);
    d^=d>>35;d^=d<<4;
    return (d+b)^c;
}
