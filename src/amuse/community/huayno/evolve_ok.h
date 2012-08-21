// Optimal Kick header
struct force {
    struct particle *parti;
    struct particle *partj;
    FLOAT timestep;
};
struct forces {
    UINT n;
    struct force *forc;
    struct force *last;
};
extern struct forces zeroforces;

void evolve_ok_init(struct sys s);
void evolve_ok_stop();
void evolve_ok2(int clevel,struct sys s, struct forces f, DOUBLE stime, DOUBLE etime, DOUBLE dt, int calc_timestep);
