#ifndef CUSTOM_RNG_HH
#define CUSTOM_RNG_HH

/** Custom random number generator based of "Ran" routine in Numerical Recipes
 * by Press et al. */
class custom_rng {
    public:
        unsigned long a,b,c;
        custom_rng(unsigned long seed);
        unsigned long int64();
        inline double doub() {
            return 5.42101086242752217E-20*int64();
        }
};

#endif
